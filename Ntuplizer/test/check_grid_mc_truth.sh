#!/usr/bin/env bash

set -euo pipefail

DEFAULT_DATASET='/LooseMuCosmic_Bin-P-10to3000-T0-Minus50to0_cosmuogen/llunerti-Run3_2024_MINI-df1e99b50d14b85be33e7e4ab518ee3a/USER'
DATASET="${1:-${DEFAULT_DATASET}}"
OUTPUT="${2:-miniAOD_products.txt}"
INSTANCE="${DAS_INSTANCE:-prod/phys03}"
REDIRECTOR="${XROOTD_REDIRECTOR:-root://cms-xrd-global.cern.ch}"
TRUTH_PATTERN='HepMCProduct|GenParticle|packedGenParticle|prunedGenParticle|SimTrack|SimVertex|TrackingParticle'

for command_name in dasgoclient edmDumpEventContent; do
    if ! command -v "${command_name}" >/dev/null 2>&1; then
        echo "ERROR: ${command_name} is unavailable. Enter a CMSSW environment and run cmsenv." >&2
        exit 1
    fi
done

echo "Dataset: ${DATASET}"
echo "DAS instance: ${INSTANCE}"

LFN="$(dasgoclient -query="file dataset=${DATASET} instance=${INSTANCE}" | sed -n '1p')"
if [[ -z "${LFN}" ]]; then
    echo "ERROR: DAS returned no files for this dataset." >&2
    exit 1
fi

PFN="${REDIRECTOR}/${LFN}"
echo "Inspecting first file: ${PFN}"
echo "Writing complete product list to: ${OUTPUT}"

edmDumpEventContent "${PFN}" > "${OUTPUT}"

echo
echo "Generator and simulation products:"
if ! grep -Ei "${TRUTH_PATTERN}" "${OUTPUT}"; then
    echo "No matching generator or simulation products were found."
    echo "The complete product list remains available in ${OUTPUT}."
fi
