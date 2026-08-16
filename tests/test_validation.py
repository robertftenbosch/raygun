"""Tests for the input-validation helpers."""
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from raygun.validation import (describe_noise_support, describe_reconstruction_quality,
                               is_nucleotide_sequence, noise_effective_fraction,
                               partition_records, reconstruction_identity)

# Chain I of PDB 4R8P: 147-mer nucleosomal DNA. Every character is also a valid
# amino-acid code, so nothing downstream complains about it.
DNA_4R8P = ("ATCGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCT"
            "GTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCGAT")
HISTONE_H4 = ("SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEH"
              "AKRKTVTAMDVVYALKRQGRTLYGFGG")
# Chain B of PDB 1LGH: 45 aa, below the reduction size of 50.
LGH_BETA = "AERSLSGLTEEEAIAVHDQFKTTFSAFIILAAVAHVLVWVWKPWF"


def rec(name, seq):
    return SeqRecord(Seq(seq), id=name, description="")


class TestNucleotideDetection:
    def test_detects_real_nucleosomal_dna(self):
        assert is_nucleotide_sequence(DNA_4R8P)

    def test_accepts_real_protein(self):
        assert not is_nucleotide_sequence(HISTONE_H4)
        assert not is_nucleotide_sequence(LGH_BETA)

    def test_short_ambiguous_peptides_are_not_flagged(self):
        # Genuinely ambiguous; flagging these would reject valid short peptides.
        assert not is_nucleotide_sequence("GATTACA")

    def test_tolerates_ambiguity_codes(self):
        assert is_nucleotide_sequence("ACGT" * 10 + "NNNN")

    def test_rna_is_detected(self):
        assert is_nucleotide_sequence("ACGU" * 10)


class TestNoiseEffectiveFraction:
    @pytest.mark.parametrize("length,expected", [
        (45, 0.0),    # below reduce_size: every window holds one residue
        (56, 0.12),   # measured: 88% of sigma is exactly zero
        (95, 0.90),   # measured: 10% zero
        (102, 1.0),
        (236, 1.0),
    ])
    def test_matches_measured_behaviour(self, length, expected):
        assert noise_effective_fraction(length) == pytest.approx(expected, abs=0.01)

    def test_full_support_at_and_above_two_residue_windows(self):
        assert noise_effective_fraction(100) == 1.0
        assert noise_effective_fraction(1000) == 1.0

    def test_degenerate_input(self):
        assert noise_effective_fraction(0) == 0.0


class TestDescribeNoiseSupport:
    def test_silent_when_noise_works(self):
        assert describe_noise_support(236) is None

    def test_flags_complete_inertness(self):
        msg = describe_noise_support(45)
        assert msg is not None and "no effect at all" in msg

    def test_flags_partial_support(self):
        msg = describe_noise_support(56)
        assert msg is not None and "12%" in msg


class TestReconstructionQuality:
    def test_perfect_reconstruction(self):
        assert reconstruction_identity(HISTONE_H4, HISTONE_H4) == 1.0

    def test_counts_mismatches(self):
        mutated = "C" + HISTONE_H4[1:]
        expected = (len(HISTONE_H4) - 1) / len(HISTONE_H4)
        assert reconstruction_identity(HISTONE_H4, mutated) == pytest.approx(expected)

    def test_length_mismatch_is_not_flattered(self):
        # Scoring only the overlap would report 1.0 for a truncated result.
        assert reconstruction_identity("ACDEFGHIKL", "ACDEF") == pytest.approx(0.5)

    def test_empty_template(self):
        assert reconstruction_identity("", "ACDE") == 0.0

    def test_silent_above_threshold(self):
        assert describe_reconstruction_quality("1ABC_1", 0.95) is None
        assert describe_reconstruction_quality("1ABC_1", 0.90) is None

    def test_warns_below_threshold(self):
        # 1Y2M_1 in the benchmark: 716 residues, reconstructed at 0.855.
        msg = describe_reconstruction_quality("1Y2M_1", 0.855)
        assert msg is not None
        assert "1Y2M_1" in msg and "0.855" in msg and "finetune" in msg.lower()

    def test_threshold_is_configurable(self):
        assert describe_reconstruction_quality("x", 0.93, threshold=0.95) is not None


class TestModelRegistry:
    """Both sampling commands must resolve model names the same way; they used
    to disagree, with raygun-sample-multiple hardcoded to 4.4M."""

    def test_default_is_the_recommended_model(self):
        from raygun.pretrained import DEFAULT_MODEL
        assert DEFAULT_MODEL == "8.8M"

    def test_registry_covers_the_released_models(self):
        from raygun.pretrained import PRETRAINED_MODELS
        assert set(PRETRAINED_MODELS) == {"2.2M", "4.4M", "8.8M"}

    def test_unknown_name_is_rejected(self):
        from raygun.pretrained import load_pretrained
        with pytest.raises(ValueError, match="Unknown Raygun model"):
            load_pretrained("16M")


class TestPartitionRecords:
    def test_rejects_dna_and_keeps_protein(self):
        accepted, rejected = partition_records(
            [rec("dna", DNA_4R8P), rec("h4", HISTONE_H4)])
        assert [r.id for r in accepted] == ["h4"]
        assert rejected[0][0].id == "dna"
        assert "nucleotide" in rejected[0][1]

    def test_reports_length_reason(self):
        accepted, rejected = partition_records([rec("beta", LGH_BETA)])
        assert accepted == []
        assert "length 45 < minlength 50" in rejected[0][1]

    def test_minlength_is_configurable(self):
        accepted, _ = partition_records([rec("beta", LGH_BETA)], minlength=40)
        assert [r.id for r in accepted] == ["beta"]

    def test_allow_nucleotides_opt_out(self):
        accepted, rejected = partition_records([rec("dna", DNA_4R8P)],
                                               allow_nucleotides=True)
        assert [r.id for r in accepted] == ["dna"] and rejected == []
