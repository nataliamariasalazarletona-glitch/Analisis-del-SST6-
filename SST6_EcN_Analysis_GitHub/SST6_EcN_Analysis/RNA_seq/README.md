# Análisis transcriptómicos

Código en R para el análisis de RNA-seq del sistema de secreción tipo VI (SST6) de *Escherichia coli* Nissle 1917. Se compararon tres réplicas biológicas de EcN WT y tres de EcN Ptet+, secuenciadas con lecturas pareadas de 2 × 150 pb.

El análisis completo se encuentra en `RNAseq_T6SS_EcN.R`.

## Procesamiento previo

Antes del análisis en R, las lecturas se procesaron con las siguientes herramientas:

| Etapa | Herramienta | Parámetros |
|---|---|---|
| Control de calidad | FastQC v0.12.1 y MultiQC v1.29 | — |
| Recorte de adaptadores y calidad | Trimmomatic v0.39 | `LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36` |
| Alineamiento | Bowtie2 v2.5.4 | Genoma de referencia CP022686 |
| Ordenamiento e indexación | SAMtools v1.22 | — |
| Cuantificación por gen | featureCounts v2.1.1 | Anotación genómica de CP022686 |

## Contenido del script

| Sección | Descripción |
|---|---|
| 1. Configuración | Rutas, muestras, genes de los loci T6SSa y T6SSc, parámetros y estilo gráfico |
| 2. Matriz de conteos y metadata | Importación de la salida de featureCounts |
| 3. Filtrado y normalización | `filterByExpr()` y normalización TMM |
| 4. PCA | Análisis de componentes principales sobre logCPM |
| 5. Dispersión y modelo | `estimateDisp()` y `glmQLFit()` |
| 6. Expresión diferencial | `glmQLFTest()` para Ptet+ frente a WT |
| 7. Anotación Gene Ontology | Términos GO del GFF y de eggNOG-mapper |
| 8. Sobrerrepresentación GO | Prueba hipergeométrica de los candidatos con mayor expresión |
| 9. Enriquecimiento por ranking | Enriquecimiento GO sobre el transcriptoma completo |
| 10. Arquitectura de los loci | Genes y fronteras intergénicas de T6SSa y del núcleo T6SSc |
| 11. Cobertura de RNA-seq | Cobertura por base normalizada (CPM) a partir de los BAM |
| 12. Continuidad intergénica | Razón de continuidad y porcentaje de bases con ≥ 5 lecturas |
| 13. Soporte paired-end | Fragmentos que atraviesan cada frontera intergénica |
| 14. Tablas | Cuadros 6 a 9 de la tesis |
| 15. Figura 10 | Respuesta transcriptómica global y enriquecimiento GO |
| 16. Figuras 11 a 16 | Cobertura, continuidad y soporte paired-end de SST6a y SST6c |
| 17. Figuras 26 a 28 | Robustez de cobertura, enriquecimiento por ranking, PCA y BCV |

## Métodos

**Expresión diferencial.** Los genes con expresión insuficiente se eliminaron con `filterByExpr()` y los conteos se normalizaron mediante TMM. La dispersión se estimó con `estimateDisp()` y la comparación Ptet+ frente a WT se evaluó con modelos lineales generalizados de cuasiverosimilitud (`glmQLFit()` y `glmQLFTest()`). Los valores p se ajustaron con Benjamini–Hochberg.

- Criterio confirmatorio: FDR < 0.05 y |log₂FC| ≥ 1.
- Candidatos exploratorios: p nominal < 0.01 y |log₂FC| ≥ 1.

**Enriquecimiento funcional.** Las asociaciones gen–GO se obtuvieron del GFF de CP022686 y de eggNOG-mapper. El análisis de sobrerrepresentación de los candidatos con mayor expresión se evaluó con una prueba hipergeométrica, considerando términos con al menos cinco genes. Como universo se usaron los genes evaluados que tenían anotación GO. El análisis basado en ranking ordenó los genes según signo(log₂FC) × −log₁₀(p) y estimó la significancia con 1 000 permutaciones (semilla 1917). En ambos casos se consideró FDR < 0.05 como criterio de significancia.

**Continuidad transcripcional.** La cobertura por base se calculó con `Rsamtools::pileup()` sin distinguir hebra y se normalizó por millón de lecturas mapeadas. Para cada región intergénica se calcularon:

- La razón entre su cobertura media y la cobertura media de los 100 pb terminales de los genes adyacentes.
- El porcentaje de bases con una profundidad mínima de cinco lecturas.
- El número de fragmentos paired-end (pares propios con TLEN > 0) que atraviesan la frontera con un anclaje de al menos 10 pb en ambos genes.

Estos análisis se interpretaron de forma descriptiva.

## Datos de entrada

Los datos no se incluyen en el repositorio. El script espera la siguiente estructura dentro de la carpeta del proyecto:

```text
RNA_Seq/
  Reference/
    *CP022686*.gff                          Anotación genómica
    sequence*.txt | .faa | .fasta           Proteínas enviadas a eggNOG-mapper
    eggNOG/*.emapper.annotations            Resultado de eggNOG-mapper
  Results/
    03_Alignment/<muestra>/<muestra>.sorted.bam   BAM ordenado
    03_Alignment/<muestra>/<muestra>.sorted.bam.bai   Índice del BAM
    03_Alignment/<muestra>/*flagstat*        Reporte de samtools flagstat
    04_Counts/*count*.txt                    Salida de featureCounts
```

Las muestras se denominan `EcN_Ct_1`, `EcN_Ct_2` y `EcN_Ct_3` (WT), y `EcN_ind_1`, `EcN_ind_2` y `EcN_ind_3` (Ptet+).

## Ejecución

1. Instalar los paquetes necesarios:

```r
install.packages(c("dplyr", "ggplot2", "ggrepel", "patchwork", "BiocManager"))
BiocManager::install(c("edgeR", "Rsamtools", "GenomicRanges", "IRanges"))
```

2. Modificar la ruta `project_dir` en la sección 1 del script para que apunte a la carpeta del proyecto.

3. Ejecutar el script completo:

```r
source("RNAseq_T6SS_EcN.R")
```

## Resultados

Los resultados se guardan en `Results/RNAseq_GitHub/`:

```text
Tables/     Resultados de expresión diferencial, GO, continuidad y paired-end; Cuadros 6 a 9
Figures/    Figuras 10 a 16 y 26 a 28 en formatos PNG y PDF
Coverage/   Cobertura por base de cada muestra y locus
```

Las figuras y los cuadros se nombran según su numeración en la tesis.
