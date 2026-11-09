#!/usr/bin/env bash

#SBATCH --job-name=AF_HP2_DUF4123
#SBATCH --partition=nukwa-l40s
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/work/USER/alphafold/logs/HP2_DUF4123_%j.out
#SBATCH --error=/work/USER/alphafold/logs/HP2_DUF4123_%j.err

set -euo pipefail

JOB="HP2_DUF4123"
INPUT="/work/USER/alphafold/input/${JOB}.fasta"
OUTPUT="/work/USER/alphafold/results/${JOB}_${SLURM_JOB_ID}"

ENV="/work/USER/colabfold/env"
PYTHON="${ENV}/bin/python"
COLABFOLD="${ENV}/bin/colabfold_batch"
export PATH="${ENV}/bin:${PATH}"

echo "============================================================"
echo "AlphaFold2-Multimer / ColabFold"
echo "Complex: HP2 + DUF4123"
echo "Chain A: HP2 = 321 aa"
echo "Chain B: DUF4123 = 279 aa"
echo "Total = 600 aa"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo "============================================================"

"${PYTHON}" --version
echo "${COLABFOLD}"
nvidia-smi

"${PYTHON}" - "${INPUT}" <<'PY'
import sys
seq="".join(x.strip() for x in open(sys.argv[1])
            if x.strip() and not x.startswith(">"))
chains=seq.split(":")
observed=[len(x) for x in chains]
expected=[321,279]
valid=set("ACDEFGHIKLMNPQRSTVWY")

print("Observed:",observed)
print("Expected:",expected)
print("Total:",sum(observed))

if len(chains)!=2 or observed!=expected:
    raise SystemExit("ERROR: FASTA incorrect")

if not all(set(c)<=valid for c in chains):
    raise SystemExit("ERROR: invalid amino-acid character")

print("FASTA OK")
PY

"${PYTHON}" - <<'PY'
import jax
print(jax.devices())
if not any(d.platform=="gpu" for d in jax.devices()):
    raise SystemExit("ERROR: JAX does not detect GPU")
print("JAX GPU OK")
PY

mkdir -p "${OUTPUT}"

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

ARCHIVE="/work/USER/alphafold/results/${JOB}_${SLURM_JOB_ID}.tar.gz"

tar -czf "${ARCHIVE}" \
    -C "$(dirname "${OUTPUT}")" \
    "$(basename "${OUTPUT}")"

echo "Archive: ${ARCHIVE}"
echo "Finished: $(date)"
