# Practical guide

How to get useful candidates out of Raygun, and how to tell when it is not
working. Everything here was measured on this fork; the numbers come from a
single RTX A6000 and from the 82-template benchmark in `benchmarks/`.

For what Raygun is and the paper behind it, see [README.md](README.md). For work
still outstanding, see [TODO.md](TODO.md).

## 1. What this tool is for

Raygun takes an existing protein and produces variants at a length you choose.
It is a *template* tool: it never designs from scratch, and the candidates stay
recognisably related to what you feed it (typically 0.7-0.9 identity when
shortening by 10%).

It works well when:

- your template is 200-500 residues,
- you want to shorten it by 5-20%,
- you have a downstream filter of your own (an assay, a domain search, a
  structure prediction) to rank candidates with.

It works poorly, or not at all, when:

- the template is under 100 residues (see §5: the noise parameter stops working),
- the template is over ~500 residues (reconstruction accuracy degrades),
- you want to *lengthen* a short template (it duplicates residues; see TODO.md),
- you expect PLL filtering alone to tell you which candidate is good (§7).

## 2. Setup

Needs a GPU with roughly 6 GB free, plus about 6 GB of disk for the model
weights, which download on first use into `~/.cache/torch/hub/checkpoints`:
2.5 GB for ESM-2 650M and 3.1 GB for the Raygun checkpoint.

```bash
uv venv --python 3.11 .venv          # or python -m venv .venv
VIRTUAL_ENV=.venv uv pip install torch --index-url https://download.pytorch.org/whl/cu128
VIRTUAL_ENV=.venv uv pip install -e .
.venv/bin/python -c "import torch; print(torch.cuda.is_available())"
```

Timings on an A6000, for calibration: the models take about 25 s to load per
run, generation runs at roughly 1.5 sequences/s for a 236-residue template, and
the full 82-template benchmark (5 generations plus PLL and alignment each) takes
under two minutes.

## 3. The workflow

```
templates  ->  diagnose  ->  generate  ->  filter
   §4            §5            §6          §7
```

The diagnose step is the one people skip. It costs one forward pass and tells
you whether the other three are worth running.

## 4. Getting templates

From the PDB, with matching length ranges written in the same run:

```bash
raygun-fetch-rcsb templates.fasta --limit 100 --resolution-max 2.0 \
    --method "X-RAY DIFFRACTION" --lengthinfo lengths.json --shrink 0.9
```

Or for specific structures: `--entries 4R8P,1LGH`.

If you bring your own FASTA, keep the record ids short and whitespace-free.
Biopython truncates ids at the first space, so a PDB header like
`>4R8P_1|Chains A, E|Histone H3.2|Xenopus laevis` becomes the id
`4R8P_1|Chains`, which will not match any key in your length JSON.
`raygun-fetch-rcsb` already writes clean `<entry>_<entity>` ids.

A PDB entry is not a protein: entries routinely contain DNA chains, and A, C, G
and T are all valid amino-acid codes, so nothing downstream objects to them.
They are dropped by default, both by `raygun-fetch-rcsb` and by the sampling
commands. Watch the counts in the log:

```
kept 82 of 112 entities (dropped -> too_short: 8, too_long: 1, duplicate_sequence: 21)
```

## 5. Diagnose before you generate

Both sampling commands now do this automatically and print one line per
template:

```
INFO  - 1GMM_1: reconstruction identity 1.000
WARN  - 1Y2M_1: zero-shot reconstruction identity is only 0.855, below the 0.90 ...
```

**Reconstruction identity** is how faithfully the model reproduces the template
at its own length with no noise. Across 82 PDB templates the median was 0.995
and 64 of them exceeded 0.99, so a low value is a real signal, not normal
variance. Below 0.90, fine-tune (`--finetune`) or use a shorter template: the
one failure in the benchmark was the longest template at 716 residues, and
identity correlates negatively with length (r = -0.52).

**Noise support** is warned about separately. The encoder compresses the
sequence to 50 positions by averaging over windows of `length // 50` residues,
and the noise it injects is scaled by the spread *within* each window. Below 100
residues those windows hold a single residue, whose spread is exactly zero:

| template length | noise actually does |
|---|---|
| under 50 | nothing at all, generation is deterministic |
| 50-99 | acts on a few positions near the termini only |
| 100-199 | works, but two-residue windows give little variation |
| 200+ | works fully |

Measured effect: templates under 100 residues produced 2.1 distinct sequences
out of 5 generations, against 4.9 for longer ones. A 52-residue template
produced the same sequence 5 times out of 5. **For short templates, vary the
target length instead of the noise** — that is the only knob that still does
anything.

## 6. Generating

One template per file:

```bash
raygun-sample-single --minlength 180 --maxlength 200 --noiseratio 0.3 \
    --num_raygun_samples_to_generate 50 template.fasta output-folder
```

Several at once, with a length range per record:

```bash
raygun-sample-multiple --lengthinfo lengths.json templates.fasta output-folder
```

Note `--lengthinfo`, not `--leninfo`.

Parameters worth understanding:

- `--model` selects the weights: `2.2M`, `4.4M` or `8.8M`, defaulting to `8.8M`,
  the one the README recommends. `--checkpoint mine.ckpt` uses a model you
  trained yourself instead, with `--checkpoint_encoders/--checkpoint_decoders`
  if it does not have the standard 12 of each. Both sampling commands and the
  benchmark take the same options, and each logs which weights it loaded:
  `Loading pretrained Raygun 8.8M`. Before this fork, `raygun-sample-multiple`
  silently used `4.4M` with no way to change it.
- `--noiseratio` between 0 and 0.5 for variation that stays close to the
  template, higher for more diversity. Useless below 100 residues (§5).
- `--num_raygun_samples_to_generate` is how many survive PLL filtering;
  `--sample_ratio` is the multiplier for how many get generated first. The
  default 10 means 500 generated to return 50. Lower it for short templates,
  where most generations come out identical anyway.
- `--penalizerepeats` is off by default. Worth turning on for charged,
  low-complexity templates: histone candidates still show single-residue runs of
  3-5 where the templates have 2-3. (The stronger drift toward basic
  homopolymers seen initially turned out to come from the 4.4M model rather than
  from the setting.)
- `--finetune` is rarely needed given median reconstruction of 0.995; use it
  when §5 tells you to.

## 7. Reading the output

```
output-folder/
  unfiltered_0.3_500.fasta      every generation
  unfiltered_0.3_500.pll.tsv    name, length, PLL, sequence
  filtered_0.3_50.fasta         the PLL-selected survivors
```

Candidate ids encode the settings: `h3_i_12_l_121_n_0.3` is generation 12 of
template `h3`, length 121, noise 0.3.

**Check how many distinct sequences you actually got.** The filtered file
contains the top N by PLL, which can be the same molecule many times over:

```bash
grep -v ">" filtered_0.3_50.fasta | sort -u | wc -l
```

**PLL ranks candidates against each other; it does not tell you a candidate is
good.** Generated sequences scored 0.174 lower per residue than their templates
on average, and only 1 of 82 templates produced a candidate that beat the
natural sequence. It is also a poor judge across different lengths, and for
membrane proteins it reflects hydrophobicity as much as plausibility. Add a
second filter suited to your protein: a PFAM/HMMER search, a hydrophobicity
check for TM segments, or a structure prediction.

## 8. Using it from Python

```python
import torch
from esm.pretrained import esm2_t33_650M_UR50D
from raygun.pretrained import load_pretrained

esmmodel, alph = esm2_t33_650M_UR50D()
bc = alph.get_batch_converter()
esmmodel = esmmodel.to(0).eval()

# load_pretrained picks a released model by name; load_model also accepts
# checkpoint="mine.ckpt" for one you trained yourself.
raymodel = load_pretrained("8.8M", return_lightning_module=False).to(0).eval()

with torch.no_grad():
    _, _, tok = bc([("egfp", sequence)])
    emb = esmmodel(tok.to(0), repr_layers=[33],
                   return_contacts=False)["representations"][33][:, 1:-1]

    # reconstruction, for the check in §5
    recon = raymodel(emb, return_logits_and_seqs=True)["generated-sequences"][0]

    # shortened variant
    out = raymodel(emb, target_lengths=torch.tensor([210], dtype=int),
                   noise=0.1, return_logits_and_seqs=True)["generated-sequences"][0]
```

The API has no length floor, so it will process sequences under 50 residues that
the command line rejects. That is deliberate — but the noise limits in §5 still
apply, and the fixed-length representation is shorter than 50 positions there.

## 9. When something looks wrong

| symptom | cause | what to do |
|---|---|---|
| fewer results than templates supplied | records rejected on length or as nucleotide | read the warning; it names every record and the reason. `--filter-minlength` to relax |
| all candidates identical | template too short for noise to work (§5) | vary the target length; lower `--sample_ratio` |
| candidates drift far from the template | reconstruction already poor | check the identity line; `--finetune` or use the folded domain only |
| `argparse` error about `--leninfo` | the flag is `--lengthinfo` | the old README was wrong |
| `KeyError` on a record name | length JSON keys must equal FASTA record ids | this now fails before generation, listing the missing names |
| poly-K/R runs in candidates | charged low-complexity template | `--penalizerepeats` |
| nonsense sequence, high PLL | a DNA chain was processed as protein | it is dropped by default now; do not pass `--allow-nucleotides` |
