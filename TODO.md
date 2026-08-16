# TODO

Open work on this fork. Every item below comes from something observed while
running the tool, not from a general wishlist; the evidence is noted so nobody
has to re-derive it. Items already handled are in `benchmarks/README.md` and the
git history.

## Model behaviour that is still unguarded

- [ ] **Magnification duplicates residues on short templates.** Generating
  *above* the template length turned `AERSLSGL...` into `AAEERRSSLL...` and
  `...IALGAAK` into `...IIAALLGGAAAK` for both 1LGH chains (45 and 56 aa). The
  model is padding by repeating characters rather than inserting anything.
  Unknown whether this is specific to short templates or also happens at 200+
  residues. Establish the range, then either warn or reject, the way the noise
  and reconstruction checks now do.
- [ ] **Recalibrate the PLL length correction.** Both sampling commands divide
  the raw PLL by `abs(-0.406 * len + 1.363)`. The provenance of those constants
  is not documented anywhere in the repo, and the divisor grows linearly with
  length, so scores are not comparable across the length range: a 45-residue
  candidate and a 500-residue candidate are being ranked on different scales.
  `benchmarks/results.tsv` has PLL per residue for 82 templates spanning 52-716
  residues and can be used to fit or falsify the relation. Do not change the
  constants on intuition.
- [ ] **Decide whether `--penalizerepeats` should default to on.** With it off,
  histone candidates drifted toward basic homopolymers: longest single-residue
  run went from 2-3 in the templates to 3-5 in candidates, K/R fraction rose
  1-4 percentage points, and `AVTKTQKKDGKKRRKTRKE` collapsed to `ATKQKKKRRRRRE`.
  Quantify the cost on well-behaved templates before flipping the default.
- [ ] **Deletions land in the structured core, not the disordered tails.** For
  histones H3, H4 and H2A the model kept the flexible N-terminal tails intact
  and cut into the histone-fold helices instead, which is backwards for
  miniaturisation. Worth investigating whether a region mask or a per-position
  weight is feasible; failing that, document it so users supply the folded
  domain as the template rather than the full chain.

## Repository and CI

- [ ] **Make `pytest tests/` runnable.** It currently hangs: `test_model_access.py`
  downloads four checkpoints (roughly 12 GB) and `test_ltraygun_training.py`
  trains. Mark them (`@pytest.mark.slow` / `network`) and deselect by default, so
  the fast offline suite (`test_validation.py`, `test_rcsb.py`, 46 tests, ~0.1 s)
  is what runs normally.
- [ ] **Fix the CI workflow.** `.github/workflows/python-tests.yml` runs
  `pytest tests/` on three Python versions, which means those same downloads on
  a GitHub runner. Depends on the item above. Note that Actions are not enabled
  on this fork yet, so no run has been recorded and the workflow's actual state
  is unverified.
- [ ] **Investigate the packaging metadata.** The build produces
  `UNKNOWN.egg-info` with `Name: UNKNOWN, Version: 0.0.0` despite
  `pyproject.toml` declaring `raygun` 0.2.5, which suggests a stale `setup.py`
  path is being picked up.
- [ ] **Delete the dead files** `setup.py.bak` and
  `tests/test_raygun_generation.py.bak`, once it is clear nothing depends on them.

## Benchmark extensions

- [ ] **Pin the diversity threshold properly.** The current figure (full
  diversity only above roughly 200 residues) rests on 5 generations per template
  at a single noise level. A noise x length sweep would turn that into a real
  threshold, and the warning text in `raygun/validation.py` should then quote
  the measured number.
- [ ] **Benchmark magnification.** Everything so far is `--shrink 0.9`. Ties
  directly into the duplication problem above.
- [ ] **Validate structurally, not just by sequence.** Reconstruction identity
  and PLL say nothing about whether a candidate folds. Running the candidates
  through ESMFold or AlphaFold and comparing to the template structure would
  test the actual claim the tool makes. The PDB templates already have
  experimental structures to compare against, which is why they were chosen.

## raygun-fetch-rcsb

- [ ] **Cache downloaded FASTA.** Re-running the same query re-downloads
  everything; a small on-disk cache keyed by entry id would make iteration
  cheaper and be politer to the RCSB servers.
- [ ] **Validate `--shrink` above 1.0.** `lengthinfo_for` accepts it and will
  produce magnification ranges, but that path is untested and depends on the
  duplication problem being resolved first.
