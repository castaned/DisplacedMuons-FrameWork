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

## Geometry-pair resolution diagnostic

The plotting script also calculates the residual used by the resolution analysis
from the new quality-matched geometric pair, without using the original tag/probe
assignment:

```text
R = (|q/pT|lower - |q/pT|upper) / |q/pT|upper
```

The upper detector track is the reference and the lower detector track is the
comparison. Both tracks pass `evt_dsa_passResolutionPair`: `pt > 20 GeV`,
`abs(eta) < 0.7`, at least 31 valid DT hits, normalized `chi2 < 5`, and
`ptError/pt < 0.5`. The original `passTagID`, `hasProbe`, and `probeID` values are
not used.

The residual is shown inclusively and in upper-track `pt` bins of 20, 30, 40, 50, 65,
85, 120, 200, and 1000 GeV. Each sufficiently populated bin is fitted with one
Gaussian over a robust central range. Separate mean and sigma summary graphs are
written versus upper-track `pt`. Comparing these results with the original
tag/probe resolution analysis tests whether the observed bias enters through the
tag/probe assignment and its asymmetric preselection.

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
for bias reduction. It additionally writes
`upper_lower_resolution_residual_inclusive.png`,
`upper_lower_resolution_residual_by_upper_pt.png`, and Gaussian mean/sigma
summaries versus upper-track `pt`; all histograms, fits, and graphs are saved in
`dsa_bias_study.root`.

Run MC and DATA into separate output directories and compare the same plots
before changing the production selection.

## MC truth study

The MC MiniAOD contains `prunedGenParticles`. In the inspected 2024 campaign,
each event has one status-1 cosmic muon near the upper detector boundary and
several status-3 bookkeeping copies. The ntuplizer therefore selects:

- the status-1 muon with the largest global `y` as the detector-entry truth; and
- the same-charge status-3 muon with the largest global `y` as the initial truth.

Candidate counts and validity flags are stored so this convention can be checked
before applying it to a different campaign. DATA does not consume generator
products; its truth branches retain their reset values.

The generator branches are:

- `gen_status1_nMuon`, `gen_entry_valid`, and entry `pdgId`, charge, `pt`, `eta`,
  `phi`, and vertex coordinates;
- `gen_status3_nMuon`, `gen_initial_valid`, and initial `pdgId`, `pt`, `eta`,
  `phi`, and vertex coordinates; and
- `gen_entry_over_initial_pt`, which measures propagation loss before the muon
  reaches the upper detector boundary.

For an opposite-side DSA pair, the ntuplizer also stores
`evt_dsa_gen_response_valid`, upper and lower response to entry truth, cosine
agreement of the upper DSA direction, and cosine agreement after reversing the
lower DSA direction.

### Reconstructed-level comparison

Run the same quality-matched geometric study independently for MC and DATA. Use:

- `pt_upper_vs_lower_zoom_logz.png` to expose the populated 20--150 GeV region;
- `mean_lower_pt_vs_upper_pt_profile.png` to test displacement from the equality line;
- `lower_over_upper_pt_vs_upper_pt.png` to measure the side response ratio; and
- `symmetric_residual_vs_average_pt.png` to remove the arbitrary denominator and
  reference-side choice.

Compare the profile trends between MC and DATA. If MC reproduces the upper/lower
shift, detector traversal, energy-loss modeling, or a common reconstruction
effect is favored. If the shift is substantially different in DATA, investigate
alignment, calibration, and side-dependent reconstruction efficiency. Repeat
the comparison in `eta`, `phi`, DT-hit, and relative-`pt`-error regions to locate
where the discrepancy enters.

### Truth-product validation

Inspect one parent EDM MC file on lxplus before changing the ntuplizer:

```bash
edmDumpEventContent input_mc.root | \
  grep -Ei 'HepMCProduct|GenParticle|SimTrack|SimVertex|TrackingParticle'
```

For the GRID `USER` dataset, the repository includes a helper that selects and
inspects the first file through the global XRootD redirector:

```bash
./check_grid_mc_truth.sh
```

The default is the 2024 `llunerti` dataset. A different dataset and output file
can be supplied explicitly:

```bash
./check_grid_mc_truth.sh '/primary/processed/USER' products_2023.txt
```

Set `DAS_INSTANCE` or `XROOTD_REDIRECTOR` in the environment only when a
non-default DAS instance or redirector is required.

New MC ntuples produce `dsa_response_to_gen_entry_comparison.png`,
`dsa_response_vs_gen_entry_pt.png`, `dsa_direction_vs_gen_entry.png`, and
`gen_entry_over_initial_pt.png`. A single detector-entry generator momentum can
identify a side-dependent reconstructed response, but it cannot by itself
separate energy loss while traversing CMS from lower-side reconstruction bias.

For that separation, prefer `SimTrack`/`SimVertex`, `TrackingParticle`, or
propagated truth states near the upper and lower muon systems. The decisive
quantities are upper-reco/upper-truth and lower-reco/lower-truth response. This
requires confirming which simulation products survive in the parent input before
implementing truth matching.
