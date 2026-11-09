#!/usr/bin/env bash

#SBATCH --job-name=AF_DUF4123_monomer
#SBATCH --partition=nukwa-l40s
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/work/USER/alphafold/logs/DUF4123_monomer_%j.out
#SBATCH --error=/work/USER/alphafold/logs/DUF4123_monomer_%j.err

set -euo pipefail

JOB="DUF4123_monomer"
EXPECTED=279

INPUT="/work/USER/alphafold/input/${JOB}.fasta"
OUTPUT="/work/USER/alphafold/results/${JOB}_${SLURM_JOB_ID}"

ENV="/work/USER/colabfold/env"
PYTHON="${ENV}/bin/python"
COLABFOLD="${ENV}/bin/colabfold_batch"

export PATH="${ENV}/bin:${PATH}"

mkdir -p /work/USER/alphafold/logs
mkdir -p "${OUTPUT}"

echo "============================================================"
echo "ColabFold monomer prediction"
echo "Protein: ${JOB}"
echo "Expected length: ${EXPECTED} aa"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo "============================================================"

"${PYTHON}" - "${INPUT}" "${EXPECTED}" <<'PY'
from pathlib import Path
import sys

f = Path(sys.argv[1])
expected = int(sys.argv[2])

seq = "".join(
    x.strip()
    for x in f.read_text().splitlines()
    if x.strip() and not x.startswith(">")
)

valid = set("ACDEFGHIKLMNPQRSTVWY")

print("Observed length:", len(seq))
print("Expected length:", expected)

if ":" in seq:
    raise SystemExit("ERROR: monomer FASTA contains ':'")

if len(seq) != expected:
    raise SystemExit(
        f"ERROR: expected {expected} aa, obtained {len(seq)}"
    )

if not set(seq) <= valid:
    raise SystemExit("ERROR: invalid amino-acid characters")

print("FASTA OK")
PY

echo
echo "JAX GPU verification:"

"${PYTHON}" - <<'PY'
import jax

print("JAX version:", jax.__version__)
print("Devices:", jax.devices())

gpu = [d for d in jax.devices() if d.platform == "gpu"]

if not gpu:
    raise SystemExit("ERROR: JAX does not detect GPU")

print("JAX GPU OK")
PY

echo
echo "============================================================"
echo "Starting AlphaFold2 pTM monomer"
echo "============================================================"

"${COLABFOLD}" \
    "${INPUT}" \
    "${OUTPUT}" \
    --data /work/USER/colabfold_cache \
    --model-type alphafold2_ptm \
    --num-models 5 \
    --num-recycle 6 \
    --num-seeds 1 \
    --random-seed 0 \
    --rank plddt

echo
echo "============================================================"
echo "MONOMER FINISHED SUCCESSFULLY"
echo "============================================================"
echo "Output directory: ${OUTPUT}"
du -sh "${OUTPUT}"
echo "Finished: $(date)"
