#!/usr/bin/env bash

#SBATCH --job-name=AF_Control_4JIV
#SBATCH --partition=nukwa-l40s
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/work/USER/Docking_Tesis/Control_Positivo/logs/Control_4JIV_%j.out
#SBATCH --error=/work/USER/Docking_Tesis/Control_Positivo/logs/Control_4JIV_%j.err

set -euo pipefail

BASE="/work/USER/Docking_Tesis/Control_Positivo"
INPUT="${BASE}/sequences/Control_4JIV.fasta"
OUTPUT="${BASE}/results/Control_4JIV_${SLURM_JOB_ID}"

ENV="/work/USER/colabfold/env"
PYTHON="${ENV}/bin/python"
COLABFOLD="${ENV}/bin/colabfold_batch"
DATA="/work/USER/colabfold_cache"

export PATH="${ENV}/bin:${PATH}"

echo "================================================"
echo "CONTROL POSITIVO ALPHAFOLD-MULTIMER"
echo "Referencia experimental: PDB 4JIV"
echo "Complejo: gp5 x 3 + PAAR"
echo "Longitudes: 96 + 96 + 96 + 93 = 381 aa"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo "================================================"

echo "Verificando FASTA..."

"${PYTHON}" - "${INPUT}" <<'PY'
import sys
from pathlib import Path

f = Path(sys.argv[1])

seq = "".join(
    line.strip()
    for line in f.read_text().splitlines()
    if line.strip() and not line.startswith(">")
)

chains = seq.split(":")
observed = [len(chain) for chain in chains]
expected = [96, 96, 96, 93]

print("Longitudes observadas:", observed)
print("Longitudes esperadas:", expected)

if observed != expected:
    raise SystemExit("ERROR: longitudes incorrectas")

if not chains[0] == chains[1] == chains[2]:
    raise SystemExit("ERROR: las tres cadenas gp5 no coinciden")

valid = set("ACDEFGHIKLMNPQRSTVWY")

if not all(set(chain) <= valid for chain in chains):
    raise SystemExit("ERROR: aminoacidos invalidos")

print("FASTA OK")
PY

echo
echo "Verificando parametros multimer..."

test -f "${DATA}/params/download_complexes_multimer_v3_finished.txt"

echo "Parametros OK"

echo
echo "Verificando GPU..."

nvidia-smi

"${PYTHON}" - <<'PY'
import jax

print("JAX version:", jax.__version__)
print("Devices:", jax.devices())

if not any(d.platform == "gpu" for d in jax.devices()):
    raise SystemExit("ERROR: JAX no detecta GPU")

print("JAX GPU OK")
PY

mkdir -p "${OUTPUT}"

echo
echo "Iniciando AlphaFold-Multimer v3..."

"${COLABFOLD}" \
    "${INPUT}" \
    "${OUTPUT}" \
    --data "${DATA}" \
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
echo "================================================"
echo "CONTROL POSITIVO FINALIZADO"
echo "================================================"
echo "Resultados: ${OUTPUT}"
du -sh "${OUTPUT}"
echo "Finished: $(date)"
