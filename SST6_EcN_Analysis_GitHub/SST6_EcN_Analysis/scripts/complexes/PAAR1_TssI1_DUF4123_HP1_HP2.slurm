#!/usr/bin/env bash

#SBATCH --job-name=AF_PAAR1_TssI1_DUF_HP1_HP2
#SBATCH --partition=nukwa-l40s
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/work/USER/alphafold/logs/PAAR1_TssI1_DUF4123_HP1_HP2_%j.out
#SBATCH --error=/work/USER/alphafold/logs/PAAR1_TssI1_DUF4123_HP1_HP2_%j.err

set -euo pipefail

JOB="PAAR1_TssI1_DUF4123_HP1_HP2"

INPUT="/work/USER/alphafold/input/${JOB}.fasta"
OUTPUT="/work/USER/alphafold/results/${JOB}_${SLURM_JOB_ID}"

ENV="/work/USER/colabfold/env"
PYTHON="${ENV}/bin/python"
COLABFOLD="${ENV}/bin/colabfold_batch"

export PATH="${ENV}/bin:${PATH}"

echo "============================================================"
echo "AlphaFold2-Multimer / ColabFold"
echo "PENTAMER"
echo "============================================================"
echo
echo "A = PAAR1        31 aa"
echo "B = TssI1       652 aa"
echo "C = DUF4123     279 aa"
echo "D = HP1         516 aa"
echo "E = HP2         321 aa"
echo "TOTAL          1799 aa"
echo
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

echo
echo "============================================================"
echo "FASTA verification"
echo "============================================================"

"${PYTHON}" - "${INPUT}" <<'PY'
import sys

fasta = sys.argv[1]

seq = "".join(
    x.strip()
    for x in open(fasta)
    if x.strip() and not x.startswith(">")
)

chains = seq.split(":")

expected = [31, 652, 279, 516, 321]
observed = [len(x) for x in chains]

valid = set("ACDEFGHIKLMNPQRSTVWY")

print("Observed:", observed)
print("Expected:", expected)
print("Chains:", len(chains))
print("Total:", sum(observed))

if len(chains) != 5:
    raise SystemExit(
        f"ERROR: expected 5 chains, obtained {len(chains)}"
    )

if observed != expected:
    raise SystemExit(
        f"ERROR: expected {expected}, obtained {observed}"
    )

if not all(set(c) <= valid for c in chains):
    raise SystemExit("ERROR: invalid amino-acid character")

print("FASTA OK")
PY

echo
echo "============================================================"
echo "JAX GPU verification"
echo "============================================================"

"${PYTHON}" - <<'PY'
import jax

print("JAX version:", jax.__version__)

devices = jax.devices()
print("Devices:", devices)

gpu = [d for d in devices if d.platform == "gpu"]

if not gpu:
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
echo "ALPHAFOLD FINISHED SUCCESSFULLY"
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

echo
echo "Finished: $(date)"
echo "============================================================"
