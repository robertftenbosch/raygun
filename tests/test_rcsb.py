"""Tests for the RCSB helpers. Pure parsing/filtering only, no network access."""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import raygun.rcsb as rcsb
from raygun.rcsb import (_parse_fasta_text, build_query, filter_records,
                         lengthinfo_for, parse_header)

# Verbatim from https://www.rcsb.org/fasta/entry/4R8P
FASTA_4R8P = """>4R8P_1|Chains A, E|Histone H3.2|Xenopus laevis (8355)
ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKT
DLRFQSSAVMALQEASEAYLVALFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA
>4R8P_5|Chain I|DNA (147-mer)|Synthetic DNA (32630)
ATCGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCC
CCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCGAT
"""


class TestParseHeader:
    def test_full_header(self):
        meta = parse_header("4R8P_1|Chains A, E|Histone H3.2|Xenopus laevis (8355)")
        assert meta["entry_id"] == "4R8P"
        assert meta["entity_id"] == "1"
        assert meta["chains"] == "A, E"
        assert meta["description"] == "Histone H3.2"
        assert meta["organism"] == "Xenopus laevis (8355)"

    def test_single_chain_header(self):
        assert parse_header("4R8P_5|Chain I|DNA (147-mer)|Synthetic DNA")["chains"] == "I"

    def test_sparse_header(self):
        meta = parse_header("1ABC_1")
        assert meta["entry_id"] == "1ABC" and meta["description"] == ""


class TestParseFastaText:
    def test_joins_wrapped_lines(self):
        parsed = list(_parse_fasta_text(FASTA_4R8P))
        assert len(parsed) == 2
        assert len(parsed[0][1]) == 135   # H3.2, wrapped over two lines
        assert len(parsed[1][1]) == 147   # nucleosomal DNA

    def test_handles_empty_input(self):
        assert list(_parse_fasta_text("")) == []


class TestFilterRecords:
    def records(self):
        seqs = list(_parse_fasta_text(FASTA_4R8P))
        return [SeqRecord(Seq(s), id=parse_header(h)["id"]) for h, s in seqs]

    def test_drops_dna_by_default(self):
        kept, stats = filter_records(self.records())
        assert [r.id for r in kept] == ["4R8P_1"]
        assert stats["nucleotide"] == 1

    def test_keeps_dna_when_asked(self):
        kept, stats = filter_records(self.records(), drop_nucleotides=False)
        assert len(kept) == 2 and stats["nucleotide"] == 0

    def test_deduplicates_identical_sequences(self):
        recs = self.records()
        twin = SeqRecord(Seq(str(recs[0].seq)), id="9XXX_1")
        kept, stats = filter_records(recs + [twin])
        assert [r.id for r in kept] == ["4R8P_1"]
        assert stats["duplicate_sequence"] == 1

    def test_length_bounds_counted_separately(self):
        recs = [SeqRecord(Seq("ACDEFGHIKL"), id="short"),
                SeqRecord(Seq("ACDEFGHIKL" * 200), id="long")]
        kept, stats = filter_records(recs)
        assert kept == []
        assert stats["too_short"] == 1 and stats["too_long"] == 1


class TestLengthinfo:
    def test_range_brackets_the_target(self):
        rec = SeqRecord(Seq("A" * 100), id="x")
        assert lengthinfo_for([rec], shrink=0.9, spread=0.03) == {"x": [87, 93]}

    def test_never_returns_an_inverted_range(self):
        rec = SeqRecord(Seq("A" * 60), id="x")
        low, high = lengthinfo_for([rec], shrink=0.5, spread=0.0)["x"]
        assert low <= high

    def test_ids_match_fasta_record_ids(self):
        # The whole point of the clean id: these keys must line up with what
        # raygun-sample-multiple reads back out of the FASTA.
        rec = SeqRecord(Seq("A" * 100), id="4R8P_1")
        assert "4R8P_1" in lengthinfo_for([rec])


class TestSearchEntries:
    """The PDB orders hits by id, which groups point mutants of one protein
    together; scattering the offsets is what keeps a sample diverse."""

    def fake_index(self, monkeypatch, total=1000):
        calls = []

        def fake_page(query, start, rows):
            calls.append((start, rows))
            return [f"E{i:04d}" for i in range(start, min(start + rows, total))], total

        monkeypatch.setattr(rcsb, "_search_page", fake_page)
        return calls

    def test_contiguous_reads_one_page(self, monkeypatch):
        calls = self.fake_index(monkeypatch)
        ids = rcsb.search_entries(limit=10, scatter=False, start=40)
        assert ids == [f"E{i:04d}" for i in range(40, 50)]
        assert calls == [(40, 10)]

    def test_scatter_spreads_across_the_result_set(self, monkeypatch):
        calls = self.fake_index(monkeypatch)
        ids = rcsb.search_entries(limit=20, chunks=10, seed=1)
        offsets = [start for start, _ in calls[1:]]     # first call probes the count
        assert len(offsets) == 10
        assert max(offsets) - min(offsets) > 100        # genuinely spread out
        assert len(ids) == 20

    def test_scatter_is_reproducible_for_a_seed(self, monkeypatch):
        self.fake_index(monkeypatch)
        first = rcsb.search_entries(limit=20, chunks=5, seed=7)
        self.fake_index(monkeypatch)
        assert rcsb.search_entries(limit=20, chunks=5, seed=7) == first

    def test_never_exceeds_limit(self, monkeypatch):
        self.fake_index(monkeypatch)
        assert len(rcsb.search_entries(limit=7, chunks=5, seed=0)) == 7

    def test_empty_result_set(self, monkeypatch):
        monkeypatch.setattr(rcsb, "_search_page", lambda q, s, r: ([], 0))
        assert rcsb.search_entries(limit=10) == []

    def test_handles_result_set_smaller_than_limit(self, monkeypatch):
        self.fake_index(monkeypatch, total=3)
        ids = rcsb.search_entries(limit=10, chunks=5, seed=0)
        assert len(ids) == len(set(ids)) <= 3


class TestBuildQuery:
    def test_protein_filter_always_present(self):
        nodes = build_query()["nodes"]
        assert len(nodes) == 1
        assert nodes[0]["parameters"]["attribute"].endswith("polymer_entity_count_protein")

    def test_optional_filters_are_added(self):
        nodes = build_query(resolution_max=2.0, method="X-RAY DIFFRACTION")["nodes"]
        attrs = [n["parameters"]["attribute"] for n in nodes]
        assert "rcsb_entry_info.resolution_combined" in attrs
        assert "exptl.method" in attrs
