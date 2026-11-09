#!/usr/bin/env bash

#SBATCH --job-name=AF_TssI2_PAAR2_HP1_HP2
#SBATCH --partition=nukwa-l40s
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/work/USER/alphafold/logs/TssI2_PAAR2_HP1_HP2_%j.out
#SBATCH --error=/work/USER/alphafold/logs/TssI2_PAAR2_HP1_HP2_%j.err

set -euo pipefail

JOB="TssI2_PAAR2_HP1_HP2"
INPUT="/work/USER/alphafold/input/${JOB}.fasta"
OUTPUT="/work/USER/alphafold/results/${JOB}_${SLURM_JOB_ID}"

ENV="/work/USER/colabfold/env"
PYTHON="${ENV}/bin/python"
COLABFOLD="${ENV}/bin/colabfold_batch"

export PATH="${ENV}/bin:${PATH}"

echo "============================================================"
echo "AlphaFold2-Multimer / ColabFold"
echo "TETRAMER 1: TssI II + PAAR II + HP1 + HP2"
echo "A = TssI II  801 aa"
echo "B = PAAR II    94 aa"
echo "C = HP1       516 aa"
echo "D = HP2       321 aa"
echo "TOTAL        1732 aa"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo "============================================================"

"${PYTHON}" --version
"${PYTHON}" -c "import sys; print('Python executable:', sys.executable)"

echo
echo "GPU:"
nvidia-smi

# FASTA verification
"${PYTHON}" - "${INPUT}" <<'PY'
import sys

seq = "".join(
    x.strip()
    for x in open(sys.argv[1])
    if x.strip() and not x.startswith(">")
)

chains = seq.split(":")
expected = [801, 94, 516, 321]
observed = [len(x) for x in chains]
valid = set("ACDEFGHIKLMNPQRSTVWY")

print("Observed:", observed)
print("Expected:", expected)
print("Total:", sum(observed))

if len(chains) != 4:
    raise SystemExit("ERROR: expected 4 chains")

if observed != expected:
    raise SystemExit(f"ERROR: expected {expected}, obtained {observed}")

if not all(set(c) <= valid for c in chains):
    raise SystemExit("ERROR: invalid amino-acid character")

print("FASTA OK")
PY

# GPU/JAX verification
"${PYTHON}" - <<'PY'
import jax
devices = jax.devices()
print("JAX devices:", devices)
if not any(d.platform == "gpu" for d in devices):
    raise SystemExit("ERROR: JAX does not detect GPU")
print("JAX GPU OK")
PY

mkdir -p "${OUTPUT}"

echo
echo "============================================================"
echo "Starting AlphaFold2-Multimer v3"
echo "============================================================"

"${COLABFOLD}" \
    "${INPUT}" \
    "${OUTPUT}" \
    --model-type alphafold2_multimer_v3 \
    --num-models 5 \
    --num-recycle 6 \
    --num-seeds 1 \
    --random-seed 0 \
    --pair-mode unpaired_paired \
    --pair-strategy greedy \
    --rank multimer \
    --calc-extra-ptm \
    --save-all

echo
echo "============================================================"
echo "Prediction finished"
echo "============================================================"

find "${OUTPUT}" -maxdepth 1 -type f -printf "%f\n" | sort
du -sh "${OUTPUT}"

ARCHIVE="/work/USER/alphafold/results/${JOB}_${SLURM_JOB_ID}.tar.gz"

tar -czf "${ARCHIVE}" \
    -C "$(dirname "${OUTPUT}")" \
    "$(basename "${OUTPUT}")"

echo "Archive: ${ARCHIVE}"
echo "Finished: $(date)"
echo "============================================================" 
