#!/usr/bin/env bash
set -euo pipefail

GENOMICS_GENERAL_COMMIT="39243440565d138f3eedd1f99b62916254706f71"
SRC="${CONDA_PREFIX}/opt/genomics_general-${GENOMICS_GENERAL_COMMIT:0:7}"

mkdir -p "${SRC}"
curl -fsSL "https://codeload.github.com/simonhmartin/genomics_general/tar.gz/${GENOMICS_GENERAL_COMMIT}" \
  | tar -xz --strip-components=1 -C "${SRC}"

chmod +x "${SRC}/phylo/phyml_sliding_windows.py"
ln -sf "${SRC}/phylo/phyml_sliding_windows.py" "${CONDA_PREFIX}/bin/phyml_sliding_windows.py"
