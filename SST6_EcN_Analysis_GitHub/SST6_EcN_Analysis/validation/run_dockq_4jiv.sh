#!/usr/bin/env bash
set -euo pipefail

DockQ 4JIV_predicted_rank001.pdb 4JIV_experimental.pdb \
  --mapping ABCD:BACD --n_cpu 4 --json DockQ_4JIV_results.json \
  | tee DockQ_4JIV_results.txt
