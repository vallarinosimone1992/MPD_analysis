#!/bin/csh

set SCRIPT_DIR = `pwd -P`
set MPD_SUITE_DIR = `cd "${SCRIPT_DIR}/.." && pwd -P`

setenv MPD_SUITE "${MPD_SUITE_DIR}"
setenv ROOT_RDF_SNAPSHOT_INFO 0

mkdir -p "${MPD_SUITE}/../RECO_DATA" "${MPD_SUITE}/../output" "${MPD_SUITE}/logs"

echo "MPD_SUITE set to ${MPD_SUITE}"
echo "RECO_DATA: ${MPD_SUITE}/../RECO_DATA"
echo "output:    ${MPD_SUITE}/../output"
echo "ROOT_RDF_SNAPSHOT_INFO=0"
