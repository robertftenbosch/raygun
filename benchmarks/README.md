# Zero-shot benchmark on RCSB templates

A check of how Raygun v0.2.5 (`raygun_8_8mil_800M`) behaves on structures pulled
straight from the PDB, rather than on SwissProt sequences. The point is to find
where the model holds up and where it quietly stops doing anything useful.

## Reproducing

```bash
raygun-fetch-rcsb benchmarks/rcsb100.fasta --limit 100 --resolution-max 2.0 \
    --method "X-RAY DIFFRACTION" --lengthinfo benchmarks/rcsb100_lengths.json --seed 42

python benchmarks/run_rcsb_benchmark.py benchmarks/rcsb100.fasta \
    benchmarks/results.tsv --samples 5 --noise 0.2 --shrink 0.9
```

100 entries yielded 112 polymer entities, of which 82 survived filtering
(8 too short, 1 too long, 21 duplicate sequences). Lengths span 52-716 residues,
median 241, across a wide range of organisms.

## Results

Reconstruction at the original length with no noise, per length class:

| length | n | median identity | min identity | identity at 90% length | unique of 5 samples |
|---|---|---|---|---|---|
| 50-99 | 7 | 1.000 | 0.981 | 0.842 | 2.1 |
| 100-199 | 21 | 1.000 | 0.976 | 0.838 | 4.6 |
| 200-399 | 39 | 0.995 | 0.984 | 0.832 | 5.0 |
| 400+ | 15 | 0.990 | 0.855 | 0.831 | 5.0 |

**Zero-shot reconstruction is excellent**: median identity 0.995 over all 82
templates, 64 of them above 0.99, and only one below the 0.90 threshold at which
the README suggests fine-tuning. That is better than the ~0.96 median the README
reports for SwissProt human and mouse, so the recommendation to skip fine-tuning
holds up on PDB entries too.

**Accuracy degrades with length.** Identity correlates negatively with sequence
length (r = -0.52). The single failure is the longest template in the set,
`1Y2M_1` at 716 residues (identity 0.855, dropping to 0.723 when generating at
90% length). If your template is long, check the reconstruction before trusting
any candidates.

**The noise parameter is inert on short templates.** The encoder derives its
sigma from the spread *within* each averaging window, and a sequence of length
`L` uses windows of `L // 50` residues, so windows holding a single residue
contribute exactly zero noise. The 7 templates below 100 residues averaged 2.1
distinct sequences out of 5 generations, against 4.9 for the rest; a 52-residue
template produced the same sequence 5 times out of 5. Sigma becomes non-zero
everywhere at 100 residues, but two-residue windows still vary little in
practice, and only above ~200 residues did every generation come out distinct.
For short templates, vary the target length instead of the noise.

**Candidates essentially never beat their template under PLL.** The generated
sequences scored 0.174 lower per residue on average, and only 1 of 82 templates
produced a candidate scoring above the natural sequence. PLL filtering is
therefore useful for ranking candidates against each other, not for judging
whether a candidate is as good as the original.

No generation in the whole run contained an `X` residue.

## Columns in results.tsv

`recon_identity` is exact positional identity at the template length;
`gen_identity_mean` is BLOSUM62 global-alignment identity of the shortened
candidates; `noise_effective` is the fraction of the fixed-length representation
that the noise parameter can perturb; `unique_samples` counts distinct sequences
among the generations.
