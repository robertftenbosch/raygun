# Copyright 2024  Kapil Devkota, Rohit Singh
# All rights reserved
# This code is available under the terms of the license available at https://github.com/rohitsinghlab/raygun
"""Input validation helpers.

Raygun silently ignores input it cannot handle: records outside the length
bounds are dropped by `RaygunData`, and nucleotide sequences are accepted as
if they were protein (ACGT are all valid amino-acid codes). Both produce a
successful-looking run with missing or meaningless results. The helpers here
make those cases detectable.
"""

# ACGT/U plus N; every one of these is also a valid amino-acid code, which is
# exactly why a nucleotide FASTA passes through ESM-2 without complaint.
NUCLEOTIDE_CHARS = set("ACGTUN")

# Below this length the noise injected by the encoder is structurally zero:
# `Reduction` averages over windows of `len(seq) // reduce_size` residues and
# uses the within-window standard deviation as sigma, so a window holding a
# single residue has sigma exactly 0. See `noise_effective_fraction`.
NOISE_RELIABLE_LENGTH = 100


def is_nucleotide_sequence(seq, threshold=0.95, minlength=20):
    """True if `seq` looks like DNA/RNA rather than protein.

    Uses a fraction rather than a strict subset check so that a sequence with a
    few ambiguity codes is still caught. Short sequences are never flagged:
    a peptide like `GATTACA` is genuinely ambiguous.
    """
    if len(seq) < minlength:
        return False
    seq = seq.upper()
    return sum(c in NUCLEOTIDE_CHARS for c in seq) / len(seq) >= threshold


def noise_effective_fraction(length, reduce_size=50):
    """Fraction of the fixed-length representation that `noise` can actually perturb.

    `Reduction` splits the sequence into `reduce_size` windows: `gap` windows of
    size `w + 1` and `reduce_size - gap` windows of size `w`, where
    `w = length // reduce_size`. Sigma is the standard deviation within a window,
    so any window covering a single residue contributes exactly zero noise.

    Returns a value in [0, 1]; 0.0 means `--noiseratio` is completely inert.
    """
    if length <= 0:
        return 0.0
    window = length // reduce_size
    if window == 0:
        # Sequence shorter than reduce_size: every window holds one residue and
        # the representation is shorter than reduce_size as well.
        return 0.0
    if window == 1:
        # Only the `gap` windows of size 2 carry a non-zero sigma.
        gap = length - window * reduce_size
        return gap / reduce_size
    return 1.0


def describe_noise_support(length, reduce_size=50):
    """Human-readable warning if `noise` is degraded at this length, else None."""
    frac = noise_effective_fraction(length, reduce_size=reduce_size)
    if frac >= 1.0:
        return None
    if frac <= 0.0:
        return (f"sequence length {length} is below the reduction size "
                f"({reduce_size}): the noise parameter has no effect at all and "
                f"generation is deterministic. Vary the target length instead.")
    return (f"sequence length {length} only supports noise on {frac:.0%} of the "
            f"fixed-length representation (positions near the termini); expect "
            f"few unique samples. Sigma is non-zero everywhere from "
            f"{NOISE_RELIABLE_LENGTH} residues upward, though windows of two "
            f"residues still give little variation; on a 82-template PDB "
            f"benchmark, sequences under 200 residues yielded noticeably fewer "
            f"distinct samples than longer ones.")


# The README's own criterion: below this zero-shot reconstruction identity it
# recommends fine-tuning before trusting any generated candidate.
RECONSTRUCTION_THRESHOLD = 0.90


def reconstruction_identity(template, reconstruction):
    """Positional identity between a template and its zero-noise reconstruction.

    Both are the same length by construction, so no alignment is needed. A short
    reconstruction (which should not happen) is scored against the template
    length rather than silently flattering the result.
    """
    if not template:
        return 0.0
    matches = sum(a == b for a, b in zip(template, reconstruction))
    return matches / max(len(template), len(reconstruction))


def describe_reconstruction_quality(name, identity,
                                    threshold=RECONSTRUCTION_THRESHOLD):
    """Warning text if the model cannot reproduce this template, else None."""
    if identity >= threshold:
        return None
    return (f"{name}: zero-shot reconstruction identity is only {identity:.3f}, "
            f"below the {threshold:.2f} at which fine-tuning is recommended. "
            f"Generated candidates start from a representation that does not "
            f"reproduce this template, so they will drift further still. "
            f"Consider --finetune, or a shorter template: reconstruction "
            f"accuracy falls off with sequence length.")


def partition_records(records, minlength=50, maxlength=1000,
                      allow_nucleotides=False):
    """Split SeqIO records into (accepted, rejected).

    `rejected` is a list of (record, reason) so callers can report exactly what
    was dropped rather than silently shrinking the dataset.
    """
    accepted, rejected = [], []
    for rec in records:
        seq = str(rec.seq)
        if not allow_nucleotides and is_nucleotide_sequence(seq):
            rejected.append((rec, "looks like a nucleotide sequence, not protein"))
        elif len(seq) < minlength:
            rejected.append((rec, f"length {len(seq)} < minlength {minlength}"))
        elif len(seq) > maxlength:
            rejected.append((rec, f"length {len(seq)} > maxlength {maxlength}"))
        else:
            accepted.append(rec)
    return accepted, rejected
