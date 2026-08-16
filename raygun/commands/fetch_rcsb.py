# Copyright 2024  Kapil Devkota, Rohit Singh
# All rights reserved
# This code is available under the terms of the license available at https://github.com/rohitsinghlab/raygun
"""Fetch template sequences from the RCSB PDB, ready for raygun-sample-multiple."""
import argparse
import json
import logging
import sys

from Bio import SeqIO

from raygun.rcsb import (fetch_entry_records, filter_records, lengthinfo_for,
                         search_entries)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
_handler = logging.StreamHandler()
_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(_handler)


def get_params():
    parser = argparse.ArgumentParser(
        description="Download protein sequences from the RCSB PDB as Raygun templates.")
    parser.add_argument("output_fasta", help="Where to write the template FASTA")
    parser.add_argument("--limit", type=int, default=100,
                        help="Number of PDB entries to fetch (default: 100)")
    parser.add_argument("--entries", default=None,
                        help="Comma-separated PDB ids to fetch instead of searching (e.g. 4R8P,1LGH)")
    parser.add_argument("--resolution-max", type=float, default=None,
                        help="Only entries at or below this resolution, in angstrom")
    parser.add_argument("--method", default=None,
                        help="Experimental method, e.g. 'X-RAY DIFFRACTION'")
    parser.add_argument("--start", type=int, default=0,
                        help="Offset into the search results; only used with --contiguous")
    parser.add_argument("--contiguous", action="store_true", default=False,
                        help="Read one contiguous page from --start instead of scattering. "
                             "PDB ids are ordered roughly chronologically, so a contiguous "
                             "page returns long runs of point mutants of the same protein")
    parser.add_argument("--chunks", type=int, default=20,
                        help="Number of random offsets to draw the sample from (default: 20)")
    parser.add_argument("--seed", type=int, default=0, help="Seed for that sampling")
    parser.add_argument("--minlength", type=int, default=50,
                        help="Drop entities shorter than this (default: 50, Raygun's own floor)")
    parser.add_argument("--maxlength", type=int, default=1000,
                        help="Drop entities longer than this (default: 1000)")
    parser.add_argument("--keep-duplicates", action="store_true", default=False,
                        help="Keep entities whose sequence is identical to one already seen")
    parser.add_argument("--keep-nucleotides", action="store_true", default=False,
                        help="Keep DNA/RNA entities. They are dropped by default because "
                             "ACGT are valid amino-acid codes and would be silently "
                             "processed as protein")
    parser.add_argument("--lengthinfo", default=None,
                        help="Also write a length-range JSON for raygun-sample-multiple")
    parser.add_argument("--shrink", type=float, default=0.9,
                        help="Target length as a fraction of the template (default: 0.9)")
    parser.add_argument("--spread", type=float, default=0.03,
                        help="Half-width of the target range, as a fraction (default: 0.03)")
    return parser.parse_args().__dict__


def main():
    config = get_params()

    if config["entries"]:
        entries = [e.strip().upper() for e in config["entries"].split(",") if e.strip()]
        logger.info(f"using {len(entries)} explicitly requested entries")
    else:
        entries = search_entries(limit=config["limit"],
                                 resolution_max=config["resolution_max"],
                                 method=config["method"],
                                 start=config["start"],
                                 scatter=not config["contiguous"],
                                 chunks=config["chunks"],
                                 seed=config["seed"])
        if not entries:
            logger.error("RCSB search returned no entries; loosen the filters")
            return 1
        logger.info(f"selected {len(entries)} entries")

    records = fetch_entry_records(entries)
    logger.info(f"{len(records)} polymer entities downloaded")

    kept, stats = filter_records(records,
                                 minlength=config["minlength"],
                                 maxlength=config["maxlength"],
                                 dedupe=not config["keep_duplicates"],
                                 drop_nucleotides=not config["keep_nucleotides"])
    dropped = ", ".join(f"{reason}: {count}" for reason, count in stats.items() if count)
    logger.info(f"kept {len(kept)} of {len(records)} entities" +
                (f" (dropped -> {dropped})" if dropped else ""))
    if not kept:
        logger.error("nothing left after filtering; loosen --minlength/--maxlength")
        return 1

    SeqIO.write(kept, config["output_fasta"], "fasta")
    logger.info(f"wrote {config['output_fasta']}")

    if config["lengthinfo"]:
        info = lengthinfo_for(kept, shrink=config["shrink"], spread=config["spread"])
        with open(config["lengthinfo"], "w") as fh:
            json.dump(info, fh, indent=4)
        logger.info(f"wrote {config['lengthinfo']} "
                    f"(targets at {config['shrink']:.0%} of template length)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
