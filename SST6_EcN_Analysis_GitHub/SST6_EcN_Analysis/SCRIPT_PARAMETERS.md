# Parámetros de predicción estructural

Las predicciones se realizaron con ColabFold y se gestionaron mediante SLURM. Los archivos `.slurm` contienen las opciones de ejecución correspondientes a cada proteína o combinación proteica.

## Configuración de las predicciones

| Parámetro | Proteínas individuales | Complejos proteicos |
|---|---|---|
| Modelo | `alphafold2_ptm` | `alphafold2_multimer_v3` |
| Modelos | 5 | 5 |
| Ciclos de reciclado | 6 | 6 |
| Número de semillas | 1 | 1 |
| Semilla aleatoria | 0 | 0 |
| Clasificación | `plddt` | `multimer` |
| Modo de emparejamiento | No aplica | `unpaired_paired` |
| Estrategia de emparejamiento | No aplica | `greedy` |

## Predicción de proteínas individuales

Los archivos `HP1_monomer.slurm`, `HP2_monomer.slurm` y `DUF4123_monomer.slurm` se encuentran en `scripts/monomers/`.

## Predicción de complejos

La carpeta `scripts/complexes/` contiene las predicciones binarias y multiproteicas entre TssI1, TssI2, DUF4123, HP1, HP2 y PAAR. Cada archivo lleva el nombre de las proteínas incluidas en la predicción.

## Evaluación del control positivo

El complejo 4JIV se predijo mediante `scripts/complexes/Control_4JIV.slurm`. La evaluación estructural se realizó con los programas de `validation/`:

- `prepare_4jiv_fasta.py`: preparación de las secuencias del complejo a partir del PDB.
- `calculate_rmsd_4jiv.py`: comparación de los átomos Cα y evaluación de las seis permutaciones de gp5.
- `run_dockq_4jiv.sh`: evaluación de interfaces con DockQ y correspondencia de cadenas `ABCD:BACD`.
