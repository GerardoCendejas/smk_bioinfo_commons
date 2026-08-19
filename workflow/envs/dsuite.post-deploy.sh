#!/usr/bin/env bash
set -euo pipefail

DSUITE_COMMIT="a547f99599d763c1760548191ea3f62cc58e8ac3"
SRC="${CONDA_PREFIX}/opt/Dsuite-${DSUITE_COMMIT:0:7}"

if [[ -x "${CONDA_PREFIX}/bin/Dsuite" ]]; then
    echo "Dsuite ya presente, nada que hacer"
    exit 0
fi

mkdir -p "${SRC}"
curl -fsSL "https://codeload.github.com/millanek/Dsuite/tar.gz/${DSUITE_COMMIT}" \
  | tar -xz --strip-components=1 -C "${SRC}"

# Flags for the compiler to find the headers and libraries in the conda environment
export CPATH="${CONDA_PREFIX}/include${CPATH:+:${CPATH}}"
export LIBRARY_PATH="${CONDA_PREFIX}/lib${LIBRARY_PATH:+:${LIBRARY_PATH}}"
export LD_RUN_PATH="${CONDA_PREFIX}/lib"

# CXX exported to allow overriding the compiler, e.g., CXX=clang++ make
make -C "${SRC}" -j "$(nproc)" CXX="${CXX:-g++}"

chmod +x "${SRC}/utils/DtriosParallel" "${SRC}/utils/dtools.py"
ln -sf "${SRC}/Build/Dsuite"         "${CONDA_PREFIX}/bin/Dsuite"
ln -sf "${SRC}/utils/DtriosParallel" "${CONDA_PREFIX}/bin/DtriosParallel"
ln -sf "${SRC}/utils/dtools.py"      "${CONDA_PREFIX}/bin/dtools.py"

# No hangs
echo "--- ldd ---"
ldd "${SRC}/Build/Dsuite" | grep -E 'libz|libstdc' || true
"${CONDA_PREFIX}/bin/Dsuite" --help >/dev/null 2>&1 || true
echo "Dsuite ${DSUITE_COMMIT:0:7} instalado en ${CONDA_PREFIX}/bin"

# This is needed for the dsuite_fbranch_plot_map rule

unset LD_LIBRARY_PATH

# Install rnaturalearthhires package from ropensci repository, not in the conda repos.
Rscript -e 'remotes::install_version("rnaturalearthhires", repos = "https://ropensci.r-universe.dev", type = "source", version = "1.0.0.9000")'
