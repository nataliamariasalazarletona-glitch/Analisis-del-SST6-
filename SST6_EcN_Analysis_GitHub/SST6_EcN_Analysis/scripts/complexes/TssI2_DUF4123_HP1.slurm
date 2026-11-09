#!/usr/bin/env bash

#SBATCH --job-name=AF_TssI2_DUF4123_HP1
#SBATCH --partition=nukwa-l40s
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/work/USER/alphafold/logs/TssI2_DUF4123_HP1_%j.out
#SBATCH --error=/work/USER/alphafold/logs/TssI2_DUF4123_HP1_%j.err

set -euo pipefail

JOB="TssI2_DUF4123_HP1"

INPUT="/work/USER/alphafold/input/${JOB}.fasta"
OUTPUT="/work/USER/alphafold/results/${JOB}_${SLURM_JOB_ID}"

ENV="/work/USER/colabfold/env"
PYTHON="${ENV}/bin/python"
COLABFOLD="${ENV}/bin/colabfold_batch"

export PATH="${ENV}/bin:${PATH}"

echo "============================================================"
echo "AlphaFold2-Multimer / ColabFold"
echo "Complex: TssI2 + DUF4123 + HP1"
echo "============================================================"

echo "A = TssI2 = 801 aa"
echo "B = DUF4123 = 279 aa"
echo "C = HP1 = 516 aa"

echo "TOTAL = 1596 aa"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"

echo "============================================================"
echo "Python"
echo "============================================================"

"${PYTHON}" --version
"${PYTHON}" -c "import sys; print('Python executable:', sys.executable)"

echo
echo "============================================================"
echo "GPU"
echo "============================================================"

nvidia-smi

echo
echo "============================================================"
echo "FASTA verification"
echo "============================================================"

"${PYTHON}" - "${INPUT}" <<'PYVERIFY'
import sys

fasta = sys.argv[1]

seq = "".join(
    x.strip()
    for x in open(fasta)
    if x.strip() and not x.startswith(">")
)

chains = seq.split(":")

expected = [801, 279, 516]
observed = [len(x) for x in chains]

valid = set("ACDEFGHIKLMNPQRSTVWY")

print("Observed:", observed)
print("Expected:", expected)
print("Chains:", len(chains))
print("Total:", sum(observed))

if len(chains) != 3:
    raise SystemExit(
        "ERROR: numero incorrecto de cadenas"
    )

if observed != expected:
    raise SystemExit(
        f"ERROR: expected {expected}, obtained {observed}"
    )

if not all(set(c) <= valid for c in chains):
    raise SystemExit(
        "ERROR: invalid amino-acid character"
    )

print("FASTA OK")
PYVERIFY

echo
echo "============================================================"
echo "JAX GPU verification"
echo "============================================================"

"${PYTHON}" - <<'PYJAX'
import jax

print("JAX version:", jax.__version__)

devices = jax.devices()
print("Devices:", devices)

gpu = [
    d for d in devices
    if d.platform == "gpu"
]

if not gpu:
    raise SystemExit(
        "ERROR: JAX does not detect GPU"
    )

print("JAX GPU OK")
PYJAX

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

echo "Output directory:"
echo "${OUTPUT}"

echo
echo "Disk usage:"
du -sh "${OUTPUT}"

echo
echo "Finished: $(date)"
echo "============================================================"
