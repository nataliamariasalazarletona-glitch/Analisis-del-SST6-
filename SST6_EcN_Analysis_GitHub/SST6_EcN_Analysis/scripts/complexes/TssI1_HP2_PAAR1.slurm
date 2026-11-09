#!/usr/bin/env bash

#SBATCH --job-name=AF_TssI1_HP2_PAAR1
#SBATCH --partition=nukwa-l40s
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/work/USER/alphafold/logs/TssI1_HP2_PAAR1_%j.out
#SBATCH --error=/work/USER/alphafold/logs/TssI1_HP2_PAAR1_%j.err

set -euo pipefail

# ============================================================
# AlphaFold2-Multimer / ColabFold
# Complex: TssI I + HP2 + PAAR I
#
# Chain A = TssI I   652 aa
# Chain B = HP2      321 aa
# Chain C = PAAR I    31 aa
# Total             1004 aa
# ============================================================

JOB="TssI1_HP2_PAAR1"

INPUT="/work/USER/alphafold/input/${JOB}.fasta"
OUTPUT="/work/USER/alphafold/results/${JOB}_${SLURM_JOB_ID}"

ENV="/work/USER/colabfold/env"
PYTHON="${ENV}/bin/python"
COLABFOLD="${ENV}/bin/colabfold_batch"

export PATH="${ENV}/bin:${PATH}"

echo "============================================================"
echo "AlphaFold2-Multimer / ColabFold"
echo "Complex: TssI I + HP2 + PAAR I"
echo "Chain A: TssI I = 652 aa"
echo "Chain B: HP2 = 321 aa"
echo "Chain C: PAAR I = 31 aa"
echo "Total = 1004 aa"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo "============================================================"

echo
echo "Python:"
"${PYTHON}" --version
"${PYTHON}" -c "import sys; print('Python executable:', sys.executable)"

echo
echo "ColabFold:"
echo "${COLABFOLD}"
"${COLABFOLD}" --help >/dev/null

echo
echo "GPU:"
nvidia-smi

# ============================================================
# Verify FASTA
# ============================================================

"${PYTHON}" - "${INPUT}" <<'PY'
import sys

fasta = sys.argv[1]

seq = "".join(
    line.strip()
    for line in open(fasta)
    if line.strip() and not line.startswith(">")
)

chains = seq.split(":")

expected = [652, 321, 31]
observed = [len(x) for x in chains]

print()
print("FASTA verification")
print("------------------")
print("Observed:", observed)
print("Expected:", expected)
print("Total:", sum(observed))

if len(chains) != 3:
    raise SystemExit(
        f"ERROR: expected 3 chains, obtained {len(chains)}"
    )

if observed != expected:
    raise SystemExit(
        f"ERROR: expected {expected}, obtained {observed}"
    )

valid = set("ACDEFGHIKLMNPQRSTVWY")

if not all(set(chain) <= valid for chain in chains):
    raise SystemExit("ERROR: invalid amino-acid character detected.")

print("FASTA OK")
PY

# ============================================================
# Verify JAX GPU
# ============================================================

echo
echo "JAX GPU verification:"

"${PYTHON}" - <<'PY'
import jax

devices = jax.devices()
print(devices)

if not any(device.platform == "gpu" for device in devices):
    raise SystemExit(
        "ERROR: JAX does not detect a GPU. Prediction cancelled."
    )

print("JAX GPU OK")
PY

# ============================================================
# Output
# ============================================================

mkdir -p "${OUTPUT}"

# ============================================================
# AlphaFold2-Multimer
# ============================================================

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

# ============================================================
# Results
# ============================================================

echo
echo "============================================================"
echo "Prediction finished"
echo "============================================================"

echo
echo "Output directory:"
echo "${OUTPUT}"

echo
echo "Files:"
find "${OUTPUT}" -maxdepth 1 -type f -printf "%f\n" | sort

echo
echo "Disk usage:"
du -sh "${OUTPUT}"

# ============================================================
# Archive
# ============================================================

ARCHIVE="/work/USER/alphafold/results/${JOB}_${SLURM_JOB_ID}.tar.gz"

tar -czf "${ARCHIVE}" \
    -C "$(dirname "${OUTPUT}")" \
    "$(basename "${OUTPUT}")"

echo
echo "Archive:"
echo "${ARCHIVE}"

echo
echo "Finished: $(date)"
echo "============================================================"
