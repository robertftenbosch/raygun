"""Benchmark Raygun zero-shot behaviour across a set of RCSB templates.

For every template this measures the three things that determine whether Raygun
is usable on it: how faithfully it reconstructs at the original length, whether
the noise parameter can actually produce variation, and what the generated
candidates look like at a shortened target length.
"""
import argparse
import numpy as np
import pandas as pd
import torch
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from esm.pretrained import esm2_t33_650M_UR50D
from raygun.pll import get_PLL
from raygun.pretrained import raygun_8_8mil_800M
from raygun.validation import noise_effective_fraction
from tqdm import tqdm

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score, aligner.extend_gap_score, aligner.mode = -11, -1, "global"


def aligned_identity(a, b):
    alignment = aligner.align(a, b)[0]
    top, bottom = (str(s) for s in alignment)
    return sum(x == y for x, y in zip(top, bottom)) / max(len(a), len(b))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta")
    ap.add_argument("out_tsv")
    ap.add_argument("--samples", type=int, default=5, help="generations per template")
    ap.add_argument("--noise", type=float, default=0.2)
    ap.add_argument("--shrink", type=float, default=0.9)
    ap.add_argument("--device", type=int, default=0)
    args = ap.parse_args()

    torch.manual_seed(0)
    dev = args.device
    esmmodel, alph = esm2_t33_650M_UR50D()
    bc = alph.get_batch_converter()
    esmmodel = esmmodel.to(dev).eval()
    raymodel = raygun_8_8mil_800M().to(dev).eval()

    records = list(SeqIO.parse(args.fasta, "fasta"))
    rows = []
    with torch.no_grad():
        for rec in tqdm(records, desc="benchmark"):
            seq = str(rec.seq)
            _, _, tok = bc([(rec.id, seq)])
            emb = esmmodel(tok.to(dev), repr_layers=[33],
                           return_contacts=False)["representations"][33][:, 1:-1]

            recon = raymodel(emb, return_logits_and_seqs=True)["generated-sequences"][0]
            recon_id = sum(a == b for a, b in zip(recon, seq)) / len(seq)

            target = max(1, int(round(len(seq) * args.shrink)))
            gens = [raymodel.get_sequences_from_fixed(
                        raymodel.encoder(emb, noise=args.noise), target)[0]
                    for _ in range(args.samples)]
            idents = [aligned_identity(seq, g) for g in gens]
            plls = [get_PLL(g, esmmodel, alph, bc) / len(g) for g in gens]

            rows.append(dict(
                id=rec.id,
                description=rec.description.split("|", 1)[-1].strip(),
                length=len(seq),
                target_length=target,
                recon_identity=recon_id,
                noise_effective=noise_effective_fraction(len(seq)),
                unique_samples=len(set(gens)),
                n_samples=args.samples,
                gen_identity_mean=float(np.mean(idents)),
                gen_identity_max=float(np.max(idents)),
                any_X=any("X" in g for g in gens),
                pll_template=get_PLL(seq, esmmodel, alph, bc) / len(seq),
                pll_gen_mean=float(np.mean(plls)),
                best_sequence=gens[int(np.argmax(plls))],
            ))

    df = pd.DataFrame(rows)
    df.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"\nwrote {args.out_tsv} ({len(df)} templates)")
    print(f"reconstruction identity: median {df.recon_identity.median():.3f}, "
          f"min {df.recon_identity.min():.3f}, "
          f"below 0.90: {(df.recon_identity < 0.90).sum()}")
    print(f"generation identity at {args.shrink:.0%}: median "
          f"{df.gen_identity_mean.median():.3f}")
    print(f"templates where noise is degraded: {(df.noise_effective < 1).sum()}")
    print(f"templates producing an X residue: {df.any_X.sum()}")


if __name__ == "__main__":
    main()
