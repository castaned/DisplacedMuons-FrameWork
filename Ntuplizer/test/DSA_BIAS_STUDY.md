# DSA ntuplizer bias study

This study is intentionally isolated on `codex/dsa-ntuplizer-bias-study`. It adds
diagnostic output without changing `passTagID`, `passProbeID`, or the existing
probe-choice logic.

## Geometry convention

`dmu_dsa_side` is derived from the midpoint of the standalone track's inner and
outer positions:

- `+1`: positive global `y` (upper detector side)
- `-1`: negative global `y` (lower detector side)
- `0`: undefined midpoint side

This replaces the earlier study-only `phi`-sign label, which described track
direction rather than detector side. If `TrackExtra` is unavailable in the input
format, the inner/outer positions remain zero and the side is set to `0`.

## Per-track branches

The existing DSA arrays are supplemented by:

- `dmu_dsa_collectionIndex`
- `dmu_dsa_p`
- `dmu_dsa_qoverp`
- `dmu_dsa_qoverpError`
- `dmu_dsa_qoverpt`, calculated as `qoverp * cosh(eta)`
- `dmu_dsa_ref{x,y,z}`
- `dmu_dsa_inner{x,y,z}`
- `dmu_dsa_outer{x,y,z}`
- `dmu_dsa_side`

## Exactly-two-DSA branches

For events with exactly two reconstructed DSA muons, the ntuplizer stores:

- collection indices, sides, `pt`, and signed `qoverpt` for both tracks
- the collection-order `pt` asymmetry and its absolute value
- `evt_dsa_residual_12` and the swapped `evt_dsa_residual_21`
- upper/lower array indices, `pt`, and signed `qoverpt` when the tracks occupy opposite sides
- `evt_dsa_residual_lower_upper`, with the upper track as the reference

The collection-order branches test ordering effects. The upper/lower branches
test detector-side effects without applying a tag/probe interpretation.

## Pair selection flags

No event is discarded by the ntuplizer. Instead, three cumulative flags preserve
the raw sample and provide reproducible matched selections:

- `evt_dsa_passRawPair`: exactly two DSA tracks, opposite geometric sides, and
  angular separation greater than 2.1 radians.
- `evt_dsa_passQualityPair`: raw pair plus, for both tracks, `pt > 20 GeV`,
  `abs(eta) < 0.7`, at least 31 valid DT hits, and normalized `chi2 < 5`.
- `evt_dsa_passResolutionPair`: quality pair plus `ptError/pt < 0.5` for both
  tracks.

The plotting script uses the quality flag for the upper/lower uncertainty and
quality comparisons. It deliberately does not apply the `ptError/pt` cut while
plotting that quantity. Residual and threshold-scan plots use the resolution flag.
For ntuples produced before these flags were added, the script reconstructs the
same selections from the stored pair indices and per-track branches.

## Validation workflow

First produce a small 2022 MC or DATA ntuple with the normal CMSSW configuration.
Then run:

```bash
python3 plot_dsa_bias_study.py \
  --input '/path/to/ntuples/*.root' \
  --outdir dsa_bias_study_2022
```

At startup, the script reports how many events pass each cumulative pair flag.
It then writes normalized collection-order and upper/lower comparisons,
upper/lower quality plots (`ptError/pt`, normalized `chi2`, and DT hits), a
`qoverpt` parameterization check, forward/swapped and deterministic random-order
residuals, a two-dimensional upper-vs-lower `pt` plot, and a common minimum-`pt`
scan at 20, 30, 40, and 50 GeV. The threshold scan always reports both the
residual mean and retained event count so a low-statistics effect is not mistaken
for bias reduction.

Run MC and DATA into separate output directories and compare the same plots
before changing the production selection.
