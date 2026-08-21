#!/usr/bin/env bash
set -euo pipefail

GENOMICS_GENERAL_COMMIT="3d10f8e9570e500524bdcc8082dd1dcb08b6946f"
SRC="${CONDA_PREFIX}/opt/genomics_general-${GENOMICS_GENERAL_COMMIT:0:7}"

mkdir -p "${SRC}"
curl -fsSL "https://codeload.github.com/GerardoCendejas/genomics_general/tar.gz/${GENOMICS_GENERAL_COMMIT}" \
  | tar -xz --strip-components=1 -C "${SRC}"

chmod +x "${SRC}/phylo/phyml_sliding_windows.py"
ln -sf "${SRC}/phylo/phyml_sliding_windows.py" "${CONDA_PREFIX}/bin/phyml_sliding_windows.py"
