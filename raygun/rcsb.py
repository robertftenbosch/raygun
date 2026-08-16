# Copyright 2024  Kapil Devkota, Rohit Singh
# All rights reserved
# This code is available under the terms of the license available at https://github.com/rohitsinghlab/raygun
"""Fetch template sequences straight from the RCSB PDB.

Two public endpoints are used, both without authentication:

* the Search API (``search.rcsb.org``) to select entries by resolution,
  experimental method and so on;
* the download service (``www.rcsb.org/fasta/entry/...``) for the sequences.

PDB FASTA headers look like ``>4R8P_1|Chains A, E|Histone H3.2|Xenopus laevis``.
Biopython truncates that id at the first whitespace (``4R8P_1|Chains``), which
silently breaks the key matching in ``raygun-sample-multiple``; the records
produced here therefore carry a clean ``<entry>_<entity>`` id instead. Nucleotide
entities are dropped, since ACGT are valid amino-acid codes and would otherwise
be processed as if they were protein.
"""
import json
import logging
import random
import time
import urllib.error
import urllib.request
from collections import OrderedDict

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from raygun.validation import is_nucleotide_sequence

logger = logging.getLogger(__name__)

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
FASTA_URL = "https://www.rcsb.org/fasta/entry/{ids}"
USER_AGENT = "raygun-rcsb/1.0 (https://github.com/rohitsinghlab/raygun)"

# The download service takes comma-separated ids; keep batches modest so a
# single failure costs little and the URL stays within sane limits.
FASTA_BATCH_SIZE = 50


def _request(url, data=None, retries=3, backoff=2.0, timeout=60):
    """GET/POST with a couple of retries; RCSB occasionally rate-limits bursts."""
    headers = {"User-Agent": USER_AGENT}
    if data is not None:
        headers["Content-Type"] = "application/json"
        data = json.dumps(data).encode()
    lasterr = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=data, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode()
        except (urllib.error.URLError, TimeoutError) as err:
            lasterr = err
            if attempt < retries - 1:
                wait = backoff ** attempt
                logger.warning(f"RCSB request failed ({err}); retrying in {wait:.0f}s")
                time.sleep(wait)
    raise RuntimeError(f"RCSB request to {url} failed after {retries} attempts: {lasterr}")


def build_query(resolution_max=None, method=None, min_entities=1):
    """Assemble a Search API query for protein-containing entries."""
    terms = [{
        "type": "terminal", "service": "text",
        "parameters": {"attribute": "rcsb_entry_info.polymer_entity_count_protein",
                       "operator": "greater_or_equal", "value": min_entities},
    }]
    if resolution_max is not None:
        terms.append({
            "type": "terminal", "service": "text",
            "parameters": {"attribute": "rcsb_entry_info.resolution_combined",
                           "operator": "less_or_equal", "value": resolution_max},
        })
    if method is not None:
        terms.append({
            "type": "terminal", "service": "text",
            "parameters": {"attribute": "exptl.method",
                           "operator": "exact_match", "value": method},
        })
    return {"type": "group", "logical_operator": "and", "nodes": terms}


def _search_page(query, start, rows):
    payload = {"query": query, "return_type": "entry",
               "request_options": {"paginate": {"start": start, "rows": rows}}}
    result = json.loads(_request(SEARCH_URL, data=payload))
    return ([hit["identifier"] for hit in result.get("result_set", [])],
            result.get("total_count", 0))


def search_entries(limit=100, resolution_max=None, method=None, start=0,
                   scatter=True, chunks=20, seed=0):
    """Return up to `limit` PDB entry ids matching the filters.

    The Search API returns hits in a stable order, which for PDB ids is
    essentially alphabetical and therefore chronological: reading from offset 0
    yields 101M, 102L, 102M, 103L ... - long runs of myoglobin and T4 lysozyme
    point mutants that share one sequence. Taking a contiguous page gives a
    sample that looks large but collapses to a handful of distinct proteins.

    With `scatter` (the default) the ids are instead drawn from `chunks` random
    offsets spread across the whole result set. Pass `scatter=False` with
    `start` to read a specific contiguous slice.
    """
    query = build_query(resolution_max=resolution_max, method=method)

    if not scatter:
        ids, total = _search_page(query, start, limit)
        logger.info(f"RCSB search matched {total} entries; took {len(ids)} "
                    f"starting at offset {start}")
        return ids

    _, total = _search_page(query, 0, 1)
    if total == 0:
        return []
    nchunks = max(1, min(chunks, limit))
    perchunk = -(-limit // nchunks)          # ceil, so the chunks cover `limit`
    rng = random.Random(seed)
    highest = max(0, total - perchunk)
    offsets = sorted(rng.sample(range(highest + 1), min(nchunks, highest + 1)))

    ids, seen = [], set()
    for offset in offsets:
        page, _ = _search_page(query, offset, perchunk)
        for pdbid in page:
            if pdbid not in seen:
                seen.add(pdbid)
                ids.append(pdbid)
    logger.info(f"RCSB search matched {total} entries; sampled {len(ids)} from "
                f"{len(offsets)} offsets spread across the result set")
    return ids[:limit]


def parse_header(header):
    """Parse a PDB FASTA header into its components.

    ``4R8P_1|Chains A, E|Histone H3.2|Xenopus laevis (8355)`` becomes a dict with
    entry_id, entity_id, chains, description and organism.
    """
    parts = header.split("|")
    ident = parts[0].strip()
    entry_id, _, entity_id = ident.partition("_")
    chains = parts[1].replace("Chains", "").replace("Chain", "").strip() if len(parts) > 1 else ""
    return {
        "id": ident,
        "entry_id": entry_id,
        "entity_id": entity_id,
        "chains": chains,
        "description": parts[2].strip() if len(parts) > 2 else "",
        "organism": parts[3].strip() if len(parts) > 3 else "",
    }


def _parse_fasta_text(text):
    """Yield (header, sequence) pairs; the download service returns plain FASTA."""
    header, chunks = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
            header, chunks = line[1:], []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def fetch_entry_records(entry_ids, batch_size=FASTA_BATCH_SIZE, pause=0.2):
    """Download the FASTA for each entry and return annotated SeqRecords."""
    records = []
    for i in range(0, len(entry_ids), batch_size):
        batch = entry_ids[i:i + batch_size]
        text = _request(FASTA_URL.format(ids=",".join(batch)))
        for header, seq in _parse_fasta_text(text):
            meta = parse_header(header)
            # `<entry>_<entity>` survives Biopython's whitespace truncation.
            rec = SeqRecord(Seq(seq), id=f"{meta['entry_id']}_{meta['entity_id']}",
                            description=f"{meta['description']} | {meta['organism']}")
            rec.annotations = dict(meta)
            records.append(rec)
        logger.info(f"fetched {min(i + batch_size, len(entry_ids))}/{len(entry_ids)} entries "
                    f"({len(records)} entities so far)")
        if pause:
            time.sleep(pause)
    return records


def filter_records(records, minlength=50, maxlength=1000, dedupe=True,
                   drop_nucleotides=True):
    """Apply the same constraints Raygun imposes, but report what was dropped.

    Returns (kept, stats). `stats` counts each rejection reason so a caller can
    show why 100 entries turned into far fewer usable templates.
    """
    stats = OrderedDict([("nucleotide", 0), ("too_short", 0), ("too_long", 0),
                         ("duplicate_sequence", 0)])
    seen, kept = set(), []
    for rec in records:
        seq = str(rec.seq)
        if drop_nucleotides and is_nucleotide_sequence(seq):
            stats["nucleotide"] += 1
        elif len(seq) < minlength:
            stats["too_short"] += 1
        elif len(seq) > maxlength:
            stats["too_long"] += 1
        elif dedupe and seq in seen:
            stats["duplicate_sequence"] += 1
        else:
            seen.add(seq)
            kept.append(rec)
    return kept, stats


def lengthinfo_for(records, shrink=0.9, spread=0.03):
    """Build the JSON length ranges `raygun-sample-multiple` expects.

    `shrink` is the target fraction of the template length, `spread` the
    half-width of the range around it.
    """
    info = {}
    for rec in records:
        centre = len(rec.seq) * shrink
        low = max(1, round(centre - len(rec.seq) * spread))
        high = max(low, round(centre + len(rec.seq) * spread))
        info[rec.id] = [low, high]
    return info
