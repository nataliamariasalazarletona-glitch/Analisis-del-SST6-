# Análisis del sistema de secreción tipo VI de *Escherichia coli* Nissle 1917

Este repositorio reúne los códigos utilizados para el análisis bioinformático del sistema de secreción tipo VI (SST6) de *Escherichia coli* Nissle 1917. Comprende la predicción de estructuras proteicas, el modelado de complejos y la evaluación estructural de un control positivo. La organización también permite integrar los análisis transcriptómicos de RNA-seq.

## Organización del repositorio

```text
RNA_seq/                 Análisis transcriptómicos
scripts/
  monomers/              Predicción de proteínas individuales
  complexes/             Predicción de complejos proteicos
validation/              Evaluación estructural del control positivo 4JIV
README.md                Descripción del proyecto
SCRIPT_PARAMETERS.md     Configuración de las predicciones
requirements.txt         Dependencias de Python para validación
```

## Predicción de estructuras y complejos

Las estructuras individuales se predijeron mediante AlphaFold2-pTM y los complejos mediante AlphaFold-Multimer v3, utilizando ColabFold. Los trabajos se ejecutaron con scripts SLURM. Los parámetros de cada tipo de análisis se describen en [Parámetros de predicción](SCRIPT_PARAMETERS.md).

Los archivos de entrada se proporcionan en formato FASTA. Para los complejos, las secuencias de las distintas cadenas se separan mediante el carácter `:`. Las rutas de entrada, salida y entorno de ejecución deben configurarse según el sistema donde se realicen los cálculos.

Ejemplo de envío de un trabajo:

```bash
bash -n scripts/complexes/TssI1_DUF4123.slurm
sbatch scripts/complexes/TssI1_DUF4123.slurm
```

## Análisis de interfaces

Las interfaces de los complejos predichos se caracterizaron mediante PDBePISA, considerando las áreas de contacto, las energías de solvatación estimadas y las interacciones intermoleculares. Las estructuras se visualizaron mediante UCSF ChimeraX.

## Control positivo 4JIV

Se utilizó el complejo cristalográfico [4JIV](https://www.rcsb.org/structure/4JIV), formado por una construcción trimérica derivada de gp5 y la proteína PAAR, para evaluar la recuperación de una estructura experimental conocida. El complejo se predijo con AlphaFold-Multimer v3 y se comparó con la referencia cristalográfica mediante RMSD de átomos Cα y DockQ v2.1.3.

La preparación de las secuencias se realiza a partir de `4JIV.pdb`:

```bash
python validation/prepare_4jiv_fasta.py
```

El cálculo de RMSD considera las seis permutaciones de las cadenas equivalentes de gp5. La evaluación con DockQ utiliza el mapeo `ABCD:BACD` entre el modelo y la referencia:

```bash
python validation/calculate_rmsd_4jiv.py
bash validation/run_dockq_4jiv.sh
```

Estos dos últimos comandos se ejecutan en el directorio que contiene `4JIV_experimental.pdb` y `4JIV_predicted_rank001.pdb`.

## Programas utilizados

ColabFold, AlphaFold2, AlphaFold-Multimer, SLURM, Biopython, DockQ, PDBePISA y UCSF ChimeraX. Las dependencias de Python se especifican en `requirements.txt`.
