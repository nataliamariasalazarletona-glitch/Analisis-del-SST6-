###############################################################################
# Análisis de RNA-seq del sistema de secreción tipo VI (SST6)
# Escherichia coli Nissle 1917 | WT (EcN_Ct) frente a Ptet+ (EcN_ind)
#
# Natalia María Salazar Letona
# Tesis de Licenciatura, Universidad del Valle de Guatemala
#
# Contenido
#   1.  Configuración
#   2.  Matriz de conteos y metadata
#   3.  Filtrado y normalización TMM
#   4.  PCA
#   5.  Dispersión y modelo quasi-likelihood
#   6.  Expresión diferencial
#   7.  Anotación Gene Ontology (GFF + eggNOG)
#   8.  Enriquecimiento GO por sobrerrepresentación
#   9.  Enriquecimiento GO basado en ranking
#   10. Arquitectura de T6SSa y del núcleo T6SSc
#   11. Cobertura de RNA-seq
#   12. Continuidad transcripcional intergénica
#   13. Soporte de fragmentos paired-end
#   14. Tablas de la tesis
#   15. Figura 10
#   16. Figuras 11 a 16
#   17. Figuras 26 a 28
###############################################################################


###############################################################################
# 1. CONFIGURACIÓN
###############################################################################

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(edgeR)
  library(Rsamtools)
  library(GenomicRanges)
  library(IRanges)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

project_dir <- "C:/Users/Natalia Salazar/Documents/RNAseq_Tesis/RNA_Seq"
setwd(project_dir)

reference_dir <- file.path(project_dir, "Reference")
counts_dir    <- file.path(project_dir, "Results", "04_Counts")
alignment_dir <- file.path(project_dir, "Results", "03_Alignment")

output_dir   <- file.path(project_dir, "Results", "RNAseq_GitHub")
tables_dir   <- file.path(output_dir, "Tables")
figures_dir  <- file.path(output_dir, "Figures")
coverage_dir <- file.path(output_dir, "Coverage")

for (d in c(tables_dir, figures_dir, coverage_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Coberturas calculadas previamente (se reutilizan si existen)
old_coverage_dir <- file.path(project_dir, "Results", "07_T6SS_Transcript_Continuity", "Coverage")

target_reference <- "NZ_CP022686.1"

sample_info <- data.frame(
  Sample    = c("EcN_Ct_1", "EcN_Ct_2", "EcN_Ct_3", "EcN_ind_1", "EcN_ind_2", "EcN_ind_3"),
  Condition = c("Control", "Control", "Control", "Induced", "Induced", "Induced"),
  Replicate = c(1, 2, 3, 1, 2, 3)
)
sample_info$Strain <- ifelse(sample_info$Condition == "Control", "WT", "Ptet+")
sample_info$Label  <- paste0(sample_info$Strain, "-", sample_info$Replicate)
sample_info$BAM    <- file.path(alignment_dir, sample_info$Sample,
                                paste0(sample_info$Sample, ".sorted.bam"))

sample_levels <- sample_info$Label

# Locus T6SSa (18 genes)
T6SSa_dictionary <- data.frame(
  Locus_Tag = c("CIW80_RS07720", "CIW80_RS07725", "CIW80_RS07730", "CIW80_RS07735",
                "CIW80_RS07740", "CIW80_RS07745", "CIW80_RS07750", "CIW80_RS07755",
                "CIW80_RS07760", "CIW80_RS07765", "CIW80_RS07770", "CIW80_RS07775",
                "CIW80_RS07780", "CIW80_RS07785", "CIW80_RS07790", "CIW80_RS07795",
                "CIW80_RS07800", "CIW80_RS07805"),
  Label = c("TssB", "TssC", "TssK", "TssL", "OmpA family", "Hcp family",
            "TssH", "Vgr family", "Tli3 immunity", "Tle3 effector", "PAAR",
            "Hypothetical", "ImcF-related", "TssA", "TssF", "TssG", "TssJ", "TssE")
)

# Núcleo T6SSc (7 genes)
T6SSc_dictionary <- data.frame(
  Locus_Tag = c("CIW80_RS00230", "CIW80_RS00235", "CIW80_RS00240", "CIW80_RS00245",
                "CIW80_RS00250", "CIW80_RS00255", "CIW80_RS00260"),
  Label = c("PAAR I", "TssI I (Vgr)", "DUF4123", "g8 / HP1",
            "g7 / HP2", "PAAR II", "TssI II (Vgr)")
)

# Genes posteriores al núcleo T6SSc (solo amplían la región de cobertura)
T6SSc_downstream_tags <- c("CIW80_RS00265", "CIW80_RS00270", "CIW80_RS25690", "CIW80_RS00290")

GO_T6SS     <- "GO:0033103"
GO_thiamine <- "GO:0009228"

confirmatory_FDR <- 0.05
exploratory_P    <- 0.01
logFC_threshold  <- 1

# Colores (Okabe-Ito)
col_strain <- c("WT" = "#999999", "Ptet+" = "#D55E00")

col_function <- c(
  "Estructura y ensamblaje"            = "#D55E00",
  "Secreción y efectores"              = "#0072B2",
  "Dominio DUF4123 (función putativa)" = "#009E73",
  "Inmunidad"                          = "#E69F00",
  "Función no determinada"             = "#999999"
)

theme_thesis <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 4),
      plot.subtitle = element_text(size = base_size - 1),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(12, 20, 12, 12)
    )
}

save_figure <- function(plot, name, width, height) {
  ggsave(file.path(figures_dir, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 400, bg = "white", limitsize = FALSE)
  ggsave(file.path(figures_dir, paste0(name, ".pdf")), plot, device = cairo_pdf,
         width = width, height = height, bg = "white", limitsize = FALSE)
}

save_table <- function(x, name) {
  write.csv(x, file.path(tables_dir, paste0(name, ".csv")), row.names = FALSE)
}

extract_GFF_attribute <- function(attributes, key) {
  pattern <- paste0("(?:^|;)", key, "=([^;]+)")
  matches <- regmatches(attributes, regexec(pattern, attributes, perl = TRUE))
  vapply(matches, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
}


###############################################################################
# 2. MATRIZ DE CONTEOS Y METADATA
###############################################################################

count_candidates <- list.files(counts_dir, pattern = "\\.(txt|tsv|csv)$",
                               full.names = TRUE, ignore.case = TRUE)
count_matches <- grep("count", basename(count_candidates), ignore.case = TRUE)

if (length(count_candidates) == 1) {
  counts_file <- count_candidates[1]
} else if (length(count_matches) == 1) {
  counts_file <- count_candidates[count_matches]
} else {
  stop("No fue posible seleccionar la matriz de conteos en ", counts_dir)
}

first_line <- readLines(counts_file, n = 1, warn = FALSE)
counts_raw <- read.delim(counts_file, header = TRUE, sep = "\t",
                         skip = ifelse(grepl("^#", first_line), 1, 0),
                         check.names = FALSE)

annotation_columns <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
sample_columns <- setdiff(names(counts_raw), annotation_columns)
stopifnot(length(sample_columns) == 6)

clean_names <- vapply(basename(sample_columns), function(x) {
  hit <- regmatches(x, regexpr("EcN_(Ct|ind)_[123]", x, ignore.case = TRUE))
  if (length(hit) == 0) return(x)
  hit <- sub("^EcN_ct_", "EcN_Ct_", hit, ignore.case = TRUE)
  sub("^EcN_ind_", "EcN_ind_", hit, ignore.case = TRUE)
}, character(1), USE.NAMES = FALSE)

counts <- as.matrix(counts_raw[, sample_columns])
storage.mode(counts) <- "numeric"
rownames(counts) <- counts_raw$Geneid
colnames(counts) <- clean_names
counts <- counts[, sample_info$Sample, drop = FALSE]

stopifnot(!anyDuplicated(rownames(counts)), !anyNA(counts), all(counts >= 0))

metadata <- data.frame(
  sample    = sample_info$Sample,
  condition = factor(sample_info$Condition, levels = c("Control", "Induced")),
  replicate = sample_info$Replicate
)
rownames(metadata) <- metadata$sample


###############################################################################
# 3. FILTRADO Y NORMALIZACIÓN TMM
###############################################################################

dge <- DGEList(counts = counts, samples = metadata)

keep <- filterByExpr(dge, group = metadata$condition)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
dge  <- normLibSizes(dge, method = "TMM")

logCPM <- cpm(dge, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)

cat("Genes después del filtrado:", nrow(dge), "\n")


###############################################################################
# 4. PCA
###############################################################################

pca_object   <- prcomp(t(logCPM), center = TRUE, scale. = FALSE)
pca_variance <- (pca_object$sdev^2 / sum(pca_object$sdev^2)) * 100

pca_df <- data.frame(
  Sample = sample_info$Sample,
  Label  = sample_info$Label,
  Strain = factor(sample_info$Strain, levels = c("WT", "Ptet+")),
  PC1    = pca_object$x[, 1],
  PC2    = pca_object$x[, 2]
)


###############################################################################
# 5. DISPERSIÓN Y MODELO QUASI-LIKELIHOOD
###############################################################################

condition <- factor(metadata$condition, levels = c("Control", "Induced"))
design <- model.matrix(~ 0 + condition)
colnames(design) <- c("Control", "Induced")
rownames(design) <- metadata$sample

dge <- estimateDisp(dge, design = design, robust = TRUE)
fit <- glmQLFit(dge, design = design, robust = TRUE)

cat("BCV común:", round(sqrt(dge$common.dispersion), 4), "\n")


###############################################################################
# 6. EXPRESIÓN DIFERENCIAL (Induced vs Control)
###############################################################################

contrast_matrix <- makeContrasts(Induced_vs_Control = Induced - Control, levels = design)
qlf <- glmQLFTest(fit, contrast = contrast_matrix[, "Induced_vs_Control"])

all_results <- topTags(qlf, n = Inf, sort.by = "PValue")$table
all_results$Gene <- rownames(all_results)
all_results <- all_results[, c("Gene", "logFC", "logCPM", "F", "PValue", "FDR")]
all_results$Locus_Tag <- sub("^gene-", "", all_results$Gene)

all_results$Confirmatory_Status <- "Not significant"
all_results$Confirmatory_Status[all_results$FDR < confirmatory_FDR &
                                  all_results$logFC >= logFC_threshold] <- "Upregulated"
all_results$Confirmatory_Status[all_results$FDR < confirmatory_FDR &
                                  all_results$logFC <= -logFC_threshold] <- "Downregulated"

all_results$Exploratory_Status <- "Not selected"
all_results$Exploratory_Status[all_results$PValue < exploratory_P &
                                 all_results$logFC >= logFC_threshold] <- "Upregulated candidate"
all_results$Exploratory_Status[all_results$PValue < exploratory_P &
                                 all_results$logFC <= -logFC_threshold] <- "Downregulated candidate"

candidate_up_genes <- all_results[all_results$Exploratory_Status == "Upregulated candidate", ]

T6SSa_DE <- merge(T6SSa_dictionary, all_results, by = "Locus_Tag", sort = FALSE)
T6SSa_DE <- T6SSa_DE[match(T6SSa_dictionary$Locus_Tag, T6SSa_DE$Locus_Tag), ]
stopifnot(nrow(T6SSa_DE) == 18, !anyNA(T6SSa_DE$logFC))

cat("Genes evaluados:", nrow(all_results),
    "| Candidatos UP:", nrow(candidate_up_genes),
    "| FDR < 0.05:", sum(all_results$FDR < confirmatory_FDR), "\n")

save_table(all_results, "01_All_Differential_Expression_Results")


###############################################################################
# 7. ANOTACIÓN GENE ONTOLOGY (GFF + eggNOG)
###############################################################################

# 7.1 GFF

gff_candidates <- list.files(reference_dir, pattern = "\\.(gff|gff3)$",
                             full.names = TRUE, ignore.case = TRUE)
gff_match <- grep("CP022686", basename(gff_candidates), ignore.case = TRUE)

if (length(gff_match) == 1) {
  gff_file <- gff_candidates[gff_match]
} else if (length(gff_candidates) == 1) {
  gff_file <- gff_candidates[1]
} else {
  stop("No fue posible identificar el GFF de CP022686 en ", reference_dir)
}

gff <- read.delim(gff_file, header = FALSE, sep = "\t", quote = "",
                  comment.char = "#", fill = TRUE)[, 1:9]
colnames(gff) <- c("SeqID", "Source", "Type", "Start", "End",
                   "Score", "Strand", "Phase", "Attributes")

# 7.2 GO del GFF

gff_GO <- gff[grepl("Ontology_term=GO:", gff$Attributes, fixed = TRUE), ]
gff_GO$Gene <- extract_GFF_attribute(gff_GO$Attributes, "locus_tag")
missing <- is.na(gff_GO$Gene)
gff_GO$Gene[missing] <- extract_GFF_attribute(gff_GO$Attributes[missing], "ID")
missing <- is.na(gff_GO$Gene)
gff_GO$Gene[missing] <- extract_GFF_attribute(gff_GO$Attributes[missing], "Parent")
gff_GO$Gene <- trimws(gff_GO$Gene)
gff_GO$Gene <- ifelse(grepl("^CIW80_RS", gff_GO$Gene), paste0("gene-", gff_GO$Gene), gff_GO$Gene)
gff_GO$Gene <- sub("^cds-", "gene-", gff_GO$Gene)

GO_split <- strsplit(extract_GFF_attribute(gff_GO$Attributes, "Ontology_term"), ",")
GFF_GO_mapping <- data.frame(Gene = rep(gff_GO$Gene, lengths(GO_split)),
                             GO = trimws(unlist(GO_split)))
GFF_GO_mapping <- unique(GFF_GO_mapping[!is.na(GFF_GO_mapping$Gene) &
                                          grepl("^GO:[0-9]{7}$", GFF_GO_mapping$GO), ])

# 7.3 eggNOG

eggnog_file <- list.files(file.path(reference_dir, "eggNOG"),
                          pattern = "emapper\\.annotations$", full.names = TRUE)
stopifnot(length(eggnog_file) == 1)

eggnog_lines <- readLines(eggnog_file, warn = FALSE)
header_index <- grep("^#query", eggnog_lines)
stopifnot(length(header_index) == 1)

eggnog <- read.delim(eggnog_file, header = TRUE, sep = "\t", skip = header_index - 1,
                     comment.char = "", quote = "", fill = TRUE, check.names = FALSE)
eggnog <- eggnog[!grepl("^#", eggnog[[1]]), ]
names(eggnog)[1] <- sub("^#", "", names(eggnog)[1])

# 7.4 Vincular eggNOG con el FASTA de proteínas enviado

fasta_candidates <- list.files(reference_dir, pattern = "\\.(faa|fasta|fa|txt)$",
                               full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
fasta_candidates <- fasta_candidates[!grepl("emapper", basename(fasta_candidates), ignore.case = TRUE)]
preferred_fasta <- grep("sequence", basename(fasta_candidates), ignore.case = TRUE)

if (length(preferred_fasta) == 1) {
  protein_fasta_file <- fasta_candidates[preferred_fasta]
} else if (length(fasta_candidates) == 1) {
  protein_fasta_file <- fasta_candidates[1]
} else {
  stop("No fue posible seleccionar el FASTA de proteínas en ", reference_dir)
}

fasta_headers <- sub("^>", "", grep("^>", readLines(protein_fasta_file, warn = FALSE), value = TRUE))
fasta_map <- data.frame(query = sub("\\s.*$", "", fasta_headers), GenBank_Locus = NA_character_)
has_locus <- grepl("\\[locus_tag=CIW80_[0-9]{5}\\]", fasta_headers)
fasta_map$GenBank_Locus[has_locus] <- sub(".*\\[locus_tag=(CIW80_[0-9]{5})\\].*", "\\1",
                                          fasta_headers[has_locus])

eggnog <- merge(eggnog, fasta_map, by = "query", all.x = TRUE, sort = FALSE)
stopifnot(mean(!is.na(eggnog$GenBank_Locus)) >= 0.95)

# 7.5 Conversión de locus GenBank a gen RefSeq presente en el GFF

gff_CDS <- gff[gff$Type == "CDS", ]
CDS_parent <- extract_GFF_attribute(gff_CDS$Attributes, "Parent")
CDS_locus  <- extract_GFF_attribute(gff_CDS$Attributes, "locus_tag")
CDS_genes  <- ifelse(!is.na(CDS_parent) & grepl("^gene-CIW80_RS", CDS_parent), CDS_parent,
                     ifelse(!is.na(CDS_locus) & grepl("^CIW80_RS", CDS_locus),
                            paste0("gene-", CDS_locus), NA))

eggnog$Gene <- paste0("gene-CIW80_RS", sub("^CIW80_", "", eggnog$GenBank_Locus))
eggnog_RNAseq <- eggnog[!is.na(eggnog$GenBank_Locus) & eggnog$Gene %in% CDS_genes &
                          eggnog$Gene %in% all_results$Gene, ]

valid_GOs <- !is.na(eggnog_RNAseq$GOs) & !trimws(eggnog_RNAseq$GOs) %in% c("", "-", "NA")
egg_split <- strsplit(eggnog_RNAseq$GOs[valid_GOs], ",")
eggnog_GO_mapping <- data.frame(Gene = rep(eggnog_RNAseq$Gene[valid_GOs], lengths(egg_split)),
                                GO = trimws(unlist(egg_split)))
eggnog_GO_mapping <- unique(eggnog_GO_mapping[grepl("^GO:[0-9]{7}$", eggnog_GO_mapping$GO), ])

# 7.6 Mapa combinado

combined_GO_mapping <- unique(rbind(GFF_GO_mapping[GFF_GO_mapping$Gene %in% all_results$Gene, ],
                                    eggnog_GO_mapping))
combined_GO_genes <- intersect(all_results$Gene, unique(combined_GO_mapping$Gene))

cat("Genes con anotación GO:", length(combined_GO_genes), "de", nrow(all_results), "\n")

save_table(combined_GO_mapping, "02_GO_Gene_Mapping_GFF_plus_eggNOG")


###############################################################################
# 8. ENRIQUECIMIENTO GO POR SOBRERREPRESENTACIÓN (prueba hipergeométrica)
###############################################################################

run_GO_ORA <- function(selected_genes, universe_genes, GO_mapping, min_genes_term = 5) {
  selected_genes <- intersect(unique(selected_genes), universe_genes)
  universe_genes <- unique(universe_genes)
  N <- length(universe_genes)
  n <- length(selected_genes)

  results <- list()
  for (go_term in unique(GO_mapping$GO)) {
    term_genes <- intersect(unique(GO_mapping$Gene[GO_mapping$GO == go_term]), universe_genes)
    genes_in_term <- intersect(selected_genes, term_genes)
    M <- length(term_genes)
    k <- length(genes_in_term)
    if (M < min_genes_term || k == 0) next

    results[[go_term]] <- data.frame(
      GO = go_term,
      Selected_Count = k, Selected_Total = n,
      Universe_Count = M, Universe_Total = N,
      Fold_Enrichment = (k / n) / (M / N),
      PValue = phyper(q = k - 1, m = M, n = N - M, k = n, lower.tail = FALSE),
      Genes = paste(genes_in_term, collapse = ";")
    )
  }

  results <- do.call(rbind, results)
  results$FDR <- p.adjust(results$PValue, method = "BH")
  results <- results[order(results$PValue), ]
  rownames(results) <- NULL
  results
}

GO_ORA_universe <- unique(combined_GO_genes)
GO_ORA_mapping  <- unique(combined_GO_mapping[combined_GO_mapping$Gene %in% GO_ORA_universe, ])

GO_ORA_UP_results <- run_GO_ORA(
  selected_genes = intersect(candidate_up_genes$Gene, GO_ORA_universe),
  universe_genes = GO_ORA_universe,
  GO_mapping = GO_ORA_mapping,
  min_genes_term = 5
)

GO_UP_significant <- GO_ORA_UP_results[GO_ORA_UP_results$FDR < 0.05, ]
print(GO_UP_significant[, c("GO", "Selected_Count", "Selected_Total",
                            "Universe_Count", "Universe_Total", "Fold_Enrichment", "FDR")])

save_table(GO_ORA_UP_results, "03_GO_ORA_UP_All_Results")

# Genes que sustentan los dos procesos enriquecidos

gff_gene_rows <- gff[gff$Type == "gene", ]
gene_names <- data.frame(
  Locus_Tag = extract_GFF_attribute(gff_gene_rows$Attributes, "locus_tag"),
  Gene_Name = extract_GFF_attribute(gff_gene_rows$Attributes, "gene")
)

driver_list <- list()
for (go_id in c(GO_T6SS, GO_thiamine)) {
  genes <- strsplit(GO_UP_significant$Genes[GO_UP_significant$GO == go_id], ";")[[1]]
  d <- all_results[match(genes, all_results$Gene), c("Gene", "Locus_Tag", "logFC", "PValue", "FDR")]
  d$GO <- go_id
  driver_list[[go_id]] <- d
}
GO_drivers <- do.call(rbind, driver_list)

GO_drivers$Gene_Label <- T6SSa_dictionary$Label[match(GO_drivers$Locus_Tag, T6SSa_dictionary$Locus_Tag)]
no_label <- is.na(GO_drivers$Gene_Label)
GO_drivers$Gene_Label[no_label] <- gene_names$Gene_Name[match(GO_drivers$Locus_Tag[no_label],
                                                              gene_names$Locus_Tag)]
no_label <- is.na(GO_drivers$Gene_Label)
GO_drivers$Gene_Label[no_label] <- GO_drivers$Locus_Tag[no_label]

save_table(GO_drivers, "04_GO_Driver_Genes")


###############################################################################
# 9. ENRIQUECIMIENTO GO BASADO EN RANKING
###############################################################################

# Ranking: signo(logFC) * -log10(P)
GO_rank_table <- all_results[is.finite(all_results$logFC) & is.finite(all_results$PValue),
                             c("Gene", "logFC", "PValue", "FDR")]
GO_rank_table$Rank_Score <- sign(GO_rank_table$logFC) *
  (-log10(pmax(GO_rank_table$PValue, .Machine$double.xmin)))
GO_rank_table$Rank_Score[GO_rank_table$logFC == 0] <- 0
GO_rank_table <- GO_rank_table[GO_rank_table$Gene %in% combined_GO_genes, ]
GO_rank_table <- GO_rank_table[!duplicated(GO_rank_table$Gene), ]
GO_rank_table <- GO_rank_table[order(GO_rank_table$Rank_Score, decreasing = TRUE), ]

GO_rank_vector <- setNames(GO_rank_table$Rank_Score, GO_rank_table$Gene)

GO_rank_mapping <- unique(combined_GO_mapping[combined_GO_mapping$Gene %in% names(GO_rank_vector), ])
GO_gene_sets <- lapply(split(GO_rank_mapping$Gene, GO_rank_mapping$GO), unique)
set_sizes <- vapply(GO_gene_sets, length, integer(1))
GO_gene_sets <- GO_gene_sets[set_sizes >= 5 & set_sizes <= 500]

calculate_ES <- function(ranked_vector, gene_set) {
  gene_set <- intersect(gene_set, names(ranked_vector))
  N  <- length(ranked_vector)
  Nh <- length(gene_set)
  if (Nh == 0 || Nh == N) return(NA_real_)

  hits <- names(ranked_vector) %in% gene_set
  hit_weights <- abs(ranked_vector)
  hit_weights[!hits] <- 0
  if (sum(hit_weights) == 0) return(0)

  running_ES <- cumsum(hit_weights / sum(hit_weights) - ifelse(hits, 0, 1 / (N - Nh)))
  if (abs(max(running_ES)) >= abs(min(running_ES))) max(running_ES) else min(running_ES)
}

GO_ES_observed <- vapply(GO_gene_sets, function(g) calculate_ES(GO_rank_vector, g), numeric(1))
GO_gene_sets   <- GO_gene_sets[!is.na(GO_ES_observed)]
GO_ES_observed <- GO_ES_observed[!is.na(GO_ES_observed)]

# Permutaciones de etiquetas génicas
n_permutations <- 1000
set.seed(1917)

GO_ES_null <- matrix(NA_real_, nrow = length(GO_gene_sets), ncol = n_permutations,
                     dimnames = list(names(GO_gene_sets), NULL))

for (perm in seq_len(n_permutations)) {
  permuted_vector <- GO_rank_vector
  names(permuted_vector) <- sample(names(GO_rank_vector), replace = FALSE)
  GO_ES_null[, perm] <- vapply(GO_gene_sets, function(g) calculate_ES(permuted_vector, g), numeric(1))
  if (perm %% 100 == 0) cat("Permutación", perm, "de", n_permutations, "\n")
}

GO_rank_results <- data.frame(GO = names(GO_ES_observed),
                              Set_Size = vapply(GO_gene_sets, length, integer(1)),
                              ES = as.numeric(GO_ES_observed),
                              NES = NA_real_, PValue = NA_real_)

for (i in seq_len(nrow(GO_rank_results))) {
  observed_ES <- GO_rank_results$ES[i]
  null_ES <- GO_ES_null[GO_rank_results$GO[i], ]

  if (observed_ES >= 0) {
    same_direction <- null_ES[null_ES >= 0]
    if (length(same_direction) == 0) next
    GO_rank_results$PValue[i] <- (sum(same_direction >= observed_ES, na.rm = TRUE) + 1) /
      (length(same_direction) + 1)
  } else {
    same_direction <- null_ES[null_ES < 0]
    if (length(same_direction) == 0) next
    GO_rank_results$PValue[i] <- (sum(same_direction <= observed_ES, na.rm = TRUE) + 1) /
      (length(same_direction) + 1)
  }
  GO_rank_results$NES[i] <- observed_ES / mean(abs(same_direction), na.rm = TRUE)
}

GO_rank_results <- GO_rank_results[!is.na(GO_rank_results$NES) & !is.na(GO_rank_results$PValue), ]
GO_rank_results$FDR <- p.adjust(GO_rank_results$PValue, method = "BH")
GO_rank_results <- GO_rank_results[order(GO_rank_results$FDR, GO_rank_results$PValue), ]

GO_target_rank <- GO_rank_results[GO_rank_results$GO %in% c(GO_T6SS, GO_thiamine), ]
print(GO_target_rank)

save_table(GO_rank_results, "05_GO_Rank_Based_Enrichment_All")


###############################################################################
# 10. ARQUITECTURA DE T6SSa Y DEL NÚCLEO T6SSc
###############################################################################

gff_genes <- gff[gff$SeqID == target_reference & gff$Type == "gene", ]
gff_genes$Locus_Tag <- extract_GFF_attribute(gff_genes$Attributes, "locus_tag")
missing <- is.na(gff_genes$Locus_Tag)
gff_genes$Locus_Tag[missing] <- sub("^gene-", "", extract_GFF_attribute(gff_genes$Attributes[missing], "ID"))

classify_function <- function(label) {
  x <- tolower(label)
  case_when(
    grepl("hp1|hp2|hypothetical|^g7|^g8", x)   ~ "Función no determinada",
    grepl("duf4123", x)                        ~ "Dominio DUF4123 (función putativa)",
    grepl("immunity|tli", x)                   ~ "Inmunidad",
    grepl("hcp|vgr|tssi|paar|effector|tle", x) ~ "Secreción y efectores",
    grepl("^tss|imcf", x)                      ~ "Estructura y ensamblaje",
    TRUE                                       ~ "Función no determinada"
  )
}

build_locus <- function(dictionary, locus_name) {
  genes <- merge(dictionary, gff_genes[, c("Locus_Tag", "Start", "End", "Strand")],
                 by = "Locus_Tag", sort = FALSE)
  stopifnot(nrow(genes) == nrow(dictionary))
  genes <- genes[order(genes$Start), ]
  genes$Locus <- locus_name
  genes$Short_Locus <- sub("^CIW80_", "", genes$Locus_Tag)
  genes$Midpoint <- (genes$Start + genes$End) / 2
  genes$Category <- classify_function(genes$Label)
  rownames(genes) <- NULL
  genes
}

T6SSa_final <- build_locus(T6SSa_dictionary, "T6SSa")
T6SSc_core  <- build_locus(T6SSc_dictionary, "T6SSc")

# Fronteras entre genes consecutivos
build_boundaries <- function(genes) {
  n <- nrow(genes)
  left  <- genes[-n, ]
  right <- genes[-1, ]
  b <- data.frame(
    Locus = genes$Locus[1],
    Boundary_Index = seq_len(n - 1),
    Left_Gene = left$Locus_Tag, Left_Label = left$Label,
    Left_Start = left$Start, Left_End = left$End,
    Right_Gene = right$Locus_Tag, Right_Label = right$Label,
    Right_Start = right$Start, Right_End = right$End,
    Gap_bp = right$Start - left$End - 1L
  )
  b$Boundary_Type <- ifelse(b$Gap_bp > 0, "Intergenic", ifelse(b$Gap_bp == 0, "Adjacent", "Overlapping"))
  b$IGR_Start <- ifelse(b$Gap_bp > 0, b$Left_End + 1L, NA)
  b$IGR_End   <- ifelse(b$Gap_bp > 0, b$Right_Start - 1L, NA)
  b$Boundary_Label <- paste(b$Left_Label, "\u2192", b$Right_Label)
  b
}

boundaries_all <- rbind(build_boundaries(T6SSa_final), build_boundaries(T6SSc_core))
intergenic_boundaries <- boundaries_all[boundaries_all$Boundary_Type == "Intergenic", ]

cat("Fronteras intergénicas T6SSa:", sum(intergenic_boundaries$Locus == "T6SSa"),
    "| T6SSc:", sum(intergenic_boundaries$Locus == "T6SSc"), "\n")

# Regiones de cobertura (1 000 pb de flanco)
T6SSc_downstream <- gff_genes[gff_genes$Locus_Tag %in% T6SSc_downstream_tags, ]
coverage_flank <- 1000

T6SS_regions <- data.frame(
  Locus = c("T6SSa", "T6SSc"),
  Start = c(max(1, min(T6SSa_final$Start) - coverage_flank),
            max(1, min(T6SSc_core$Start) - coverage_flank)),
  End   = c(max(T6SSa_final$End) + coverage_flank,
            max(c(T6SSc_core$End, T6SSc_downstream$End)) + coverage_flank)
)


###############################################################################
# 11. COBERTURA DE RNA-seq
###############################################################################

stopifnot(all(file.exists(sample_info$BAM)), all(file.exists(paste0(sample_info$BAM, ".bai"))))

# 11.1 Lecturas mapeadas (samtools flagstat)

read_flagstat_value <- function(file, pattern) {
  hit <- grep(pattern, readLines(file, warn = FALSE), value = TRUE)
  as.numeric(sub("^([0-9]+).*", "\\1", hit[1]))
}

sample_info$Flagstat <- vapply(sample_info$Sample, function(s) {
  f <- list.files(file.path(alignment_dir, s), pattern = "flagstat",
                  full.names = TRUE, ignore.case = TRUE)
  if (length(f) == 0) NA_character_ else f[1]
}, character(1))
stopifnot(!anyNA(sample_info$Flagstat))

sample_info$Total_Reads  <- vapply(sample_info$Flagstat, read_flagstat_value, numeric(1), pattern = "in total")
sample_info$Mapped_Reads <- vapply(sample_info$Flagstat, read_flagstat_value, numeric(1), pattern = "primary mapped")

# 11.2 Cobertura por base (pileup sin distinguir hebra)

pileup_param <- PileupParam(distinguish_strands = FALSE, distinguish_nucleotides = FALSE,
                            include_deletions = FALSE, include_insertions = FALSE,
                            min_base_quality = 0, min_mapq = 0)

compute_coverage <- function(bam, start, end, mapped_reads) {
  region <- GRanges(target_reference, IRanges(start, end))
  pile <- pileup(bam,
                 scanBamParam = ScanBamParam(which = region,
                                             flag = scanBamFlag(isUnmappedQuery = FALSE,
                                                                isSecondaryAlignment = FALSE,
                                                                isSupplementaryAlignment = FALSE)),
                 pileupParam = pileup_param)

  cov <- data.frame(Position = seq.int(start, end), Raw_Coverage = 0L)
  if (nrow(pile) > 0) {
    by_position <- aggregate(count ~ pos, data = pile, FUN = sum)
    idx <- match(by_position$pos, cov$Position)
    cov$Raw_Coverage[idx[!is.na(idx)]] <- by_position$count[!is.na(idx)]
  }
  cov$CPM_Coverage <- cov$Raw_Coverage / (mapped_reads / 1e6)
  cov
}

coverage_list <- list()

for (i in seq_len(nrow(sample_info))) {
  for (j in seq_len(nrow(T6SS_regions))) {
    s <- sample_info$Sample[i]
    locus <- T6SS_regions$Locus[j]
    file_name <- paste0(s, "_", locus, "_Pileup_Coverage.rds")

    old_file <- file.path(old_coverage_dir, file_name)
    new_file <- file.path(coverage_dir, file_name)

    if (file.exists(new_file)) {
      cov <- readRDS(new_file)
    } else if (file.exists(old_file)) {
      cov <- readRDS(old_file)
    } else {
      cat("Calculando cobertura:", s, locus, "\n")
      cov <- compute_coverage(sample_info$BAM[i], T6SS_regions$Start[j],
                              T6SS_regions$End[j], sample_info$Mapped_Reads[i])
      saveRDS(cov, new_file)
    }

    cov <- cov[, c("Position", "Raw_Coverage", "CPM_Coverage")]
    cov$Sample <- s
    cov$Label  <- sample_info$Label[i]
    cov$Strain <- sample_info$Strain[i]
    cov$Locus  <- locus
    coverage_list[[paste(s, locus)]] <- cov
  }
}

# 11.3 Perfil medio por condición (media móvil de 101 pb)

moving_average <- function(x, window = 101L) {
  smoothed <- as.numeric(stats::filter(x, rep(1 / window, window), sides = 2))
  smoothed[is.na(smoothed)] <- x[is.na(smoothed)]
  smoothed
}

coverage_all <- do.call(rbind, lapply(coverage_list, function(cov) {
  cov$CPM_Smoothed <- moving_average(cov$CPM_Coverage, 101L)
  cov
}))

coverage_profile <- aggregate(cbind(CPM_Coverage, CPM_Smoothed) ~ Strain + Locus + Position,
                              data = coverage_all, FUN = mean)


###############################################################################
# 12. CONTINUIDAD TRANSCRIPCIONAL INTERGÉNICA
###############################################################################
# Razón = cobertura media de la región intergénica /
#         media de la cobertura de los 100 pb terminales de cada gen vecino

gene_window <- 100L
continuity_list <- list()

for (cov in coverage_list) {
  b_locus <- intergenic_boundaries[intergenic_boundaries$Locus == cov$Locus[1], ]

  for (j in seq_len(nrow(b_locus))) {
    b <- b_locus[j, ]

    left_start  <- max(b$Left_Start, b$Left_End - gene_window + 1L)
    right_end   <- min(b$Right_End, b$Right_Start + gene_window - 1L)

    left_mean  <- mean(cov$Raw_Coverage[cov$Position >= left_start & cov$Position <= b$Left_End])
    right_mean <- mean(cov$Raw_Coverage[cov$Position >= b$Right_Start & cov$Position <= right_end])
    neighbor_mean <- mean(c(left_mean, right_mean))

    igr <- cov$Raw_Coverage[cov$Position >= b$IGR_Start & cov$Position <= b$IGR_End]
    stopifnot(length(igr) == b$Gap_bp)

    continuity_list[[length(continuity_list) + 1]] <- data.frame(
      Sample = cov$Sample[1], Label = cov$Label[1], Strain = cov$Strain[1],
      Locus = b$Locus, Boundary_Index = b$Boundary_Index,
      Boundary_Label = b$Boundary_Label, Gap_bp = b$Gap_bp,
      Left_Mean_Raw = left_mean, Right_Mean_Raw = right_mean,
      IGR_Mean_Raw = mean(igr),
      IGR_Percent_GT0 = mean(igr > 0) * 100,
      IGR_Percent_GE5 = mean(igr >= 5) * 100,
      Continuity_Ratio = if (neighbor_mean > 0) mean(igr) / neighbor_mean else NA_real_
    )
  }
}

continuity <- do.call(rbind, continuity_list)
save_table(continuity, "06_Intergenic_Continuity_All_Samples")


###############################################################################
# 13. SOPORTE DE FRAGMENTOS PAIRED-END
###############################################################################
# Fragmento = par propio con TLEN > 0. Sustenta la frontera si se ancla
# al menos 10 pb en ambos genes vecinos.

anchor_bp <- 10L
query_padding <- 500L
paired_list <- list()

for (i in seq_len(nrow(sample_info))) {
  for (j in seq_len(nrow(intergenic_boundaries))) {
    b <- intergenic_boundaries[j, ]

    required_left  <- b$Left_End - anchor_bp + 1L
    required_right <- b$Right_Start + anchor_bp - 1L

    query <- GRanges(target_reference,
                     IRanges(max(1L, b$Left_End - query_padding), b$Right_Start + query_padding))

    bam_local <- scanBam(sample_info$BAM[i], param = ScanBamParam(
      which = query, what = c("flag", "pos", "isize"),
      flag = scanBamFlag(isUnmappedQuery = FALSE, hasUnmappedMate = FALSE,
                         isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE)
    ))[[1]]

    valid <- bitwAnd(bam_local$flag, 0x1) != 0 & bitwAnd(bam_local$flag, 0x2) != 0 &
      !is.na(bam_local$isize) & bam_local$isize > 0

    fragment_start <- bam_local$pos[valid]
    fragment_end   <- bam_local$pos[valid] + bam_local$isize[valid] - 1L

    touching   <- sum(fragment_start <= required_right & fragment_end >= required_left)
    supporting <- sum(fragment_start <= required_left & fragment_end >= required_right)

    paired_list[[length(paired_list) + 1]] <- data.frame(
      Sample = sample_info$Sample[i], Label = sample_info$Label[i], Strain = sample_info$Strain[i],
      Locus = b$Locus, Boundary_Index = b$Boundary_Index,
      Boundary_Label = b$Boundary_Label, Gap_bp = b$Gap_bp,
      Boundary_Touching_Fragments = touching,
      Boundary_Supporting_Fragments = supporting,
      Support_Fraction = if (touching > 0) supporting / touching else NA_real_
    )
  }
}

paired_support <- do.call(rbind, paired_list)
save_table(paired_support, "07_PairedEnd_Boundary_Support_All_Samples")


###############################################################################
# 14. TABLAS DE LA TESIS
###############################################################################

# Cuadro 6. Alineamiento (lecturas totales y alineadas)
save_table(data.frame(Biblioteca = sample_info$Label,
                      Lecturas_totales = sample_info$Total_Reads,
                      Lecturas_alineadas = sample_info$Mapped_Reads,
                      Alineamiento_porcentaje = round(100 * sample_info$Mapped_Reads /
                                                        sample_info$Total_Reads, 2)),
           "Cuadro_06_Alineamiento")

# Cuadro 7. Expresión diferencial de T6SSa
save_table(T6SSa_DE[, c("Locus_Tag", "Label", "logFC", "logCPM", "PValue", "FDR")],
           "Cuadro_07_T6SSa_Expresion_Diferencial")

# Cuadro 8. Enriquecimiento GO
save_table(GO_UP_significant[, c("GO", "Selected_Count", "Selected_Total", "Universe_Count",
                                 "Universe_Total", "Fold_Enrichment", "PValue", "FDR")],
           "Cuadro_08_GO_Enriquecimiento")

# Cuadro 9. Fragmentos de soporte por réplica
cuadro9 <- reshape(paired_support[, c("Locus", "Boundary_Index", "Boundary_Label", "Gap_bp",
                                      "Label", "Boundary_Supporting_Fragments")],
                   idvar = c("Locus", "Boundary_Index", "Boundary_Label", "Gap_bp"),
                   timevar = "Label", direction = "wide")
names(cuadro9) <- sub("Boundary_Supporting_Fragments\\.", "", names(cuadro9))
save_table(cuadro9, "Cuadro_09_PairedEnd")


###############################################################################
# 15. FIGURA 10. RESPUESTA TRANSCRIPTÓMICA Y ENRIQUECIMIENTO GO
###############################################################################

xlab_fc <- expression(log[2]~"cambio de expresión (Ptet+ / WT)")

# 15.1 Panel A. Gráfico de volcán

volcano <- all_results
volcano$MinusLog10P <- -log10(pmax(volcano$PValue, .Machine$double.xmin))
volcano$Class <- factor(
  recode(volcano$Exploratory_Status,
         "Not selected" = "No seleccionado",
         "Downregulated candidate" = "Candidato con menor expresión",
         "Upregulated candidate" = "Candidato con mayor expresión"),
  levels = c("No seleccionado", "Candidato con menor expresión", "Candidato con mayor expresión")
)
volcano$Label <- T6SSa_dictionary$Label[match(volcano$Locus_Tag, T6SSa_dictionary$Locus_Tag)]
volcano_T6SSa <- volcano[!is.na(volcano$Label), ]

fig10A <- ggplot(volcano, aes(logFC, MinusLog10P)) +
  geom_point(aes(color = Class), size = 1.5, alpha = 0.65) +
  geom_point(data = volcano_T6SSa, shape = 21, fill = NA, color = "black", stroke = 1.6, size = 4.8) +
  geom_text_repel(data = volcano_T6SSa, aes(label = Label), size = 5,
                  box.padding = 0.6, point.padding = 0.5, force = 4, force_pull = 0.3,
                  max.overlaps = Inf, min.segment.length = 0, segment.color = "grey40", seed = 123) +
  scale_color_manual(values = c("No seleccionado" = "#999999",
                                "Candidato con menor expresión" = "#0072B2",
                                "Candidato con mayor expresión" = "#D55E00")) +
  labs(title = "Respuesta transcriptómica global de Ptet+ frente a WT",
       subtitle = "Candidatos exploratorios; ningún gen alcanzó FDR < 0.05",
       x = xlab_fc, y = expression(-log[10]~"valor p nominal"),
       color = "Clasificación exploratoria") +
  theme_thesis(15)

# 15.2 Panel B. Genes asociados a los procesos GO

drivers_plot <- GO_drivers
drivers_plot$Process <- ifelse(drivers_plot$GO == GO_T6SS, "Secreción por T6SS", "Biosíntesis de tiamina")
drivers_plot <- drivers_plot[order(drivers_plot$GO != GO_T6SS, drivers_plot$logFC), ]
drivers_plot$Gene_Label <- factor(drivers_plot$Gene_Label, levels = rev(unique(drivers_plot$Gene_Label)))
n_thiamine <- sum(drivers_plot$GO == GO_thiamine)

fig10B <- ggplot(drivers_plot, aes(logFC, Gene_Label, color = Process)) +
  geom_segment(aes(x = 0, xend = logFC, yend = Gene_Label), color = "grey70", linewidth = 0.6) +
  geom_point(size = 5.5) +
  geom_text(aes(label = sprintf("%.2f", logFC)), hjust = -0.35, color = "black", size = 5.2) +
  geom_hline(yintercept = n_thiamine + 0.5, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c("Secreción por T6SS" = "#D55E00",
                                "Biosíntesis de tiamina" = "#0072B2")) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.2))) +
  labs(title = "Genes asociados a GO",
       subtitle = "Contraste Ptet+ / WT; FDR génico \u2265 0.05",
       x = xlab_fc, y = "Gen", color = "Proceso") +
  theme_thesis(15)

# 15.3 Panel C. Enriquecimiento por sobrerrepresentación

ora_plot <- GO_UP_significant[GO_UP_significant$GO %in% c(GO_T6SS, GO_thiamine), ]
ora_plot$Process <- ifelse(ora_plot$GO == GO_T6SS,
                           "Secreción de proteínas\nmediante T6SS", "Biosíntesis de tiamina")
ora_plot$Process <- reorder(ora_plot$Process, ora_plot$Fold_Enrichment)
ora_plot$Text <- sprintf("%.2f veces | %d/%d genes | FDR = %.2e",
                         ora_plot$Fold_Enrichment, ora_plot$Selected_Count,
                         ora_plot$Selected_Total, ora_plot$FDR)

fig10C <- ggplot(ora_plot, aes(Fold_Enrichment, Process)) +
  geom_segment(aes(x = 0, xend = Fold_Enrichment, yend = Process), color = "grey65", linewidth = 0.8) +
  geom_point(aes(size = Selected_Count), shape = 21, fill = "#D55E00", color = "black", stroke = 0.5) +
  geom_text(aes(label = Text), hjust = -0.1, size = 5) +
  scale_size_continuous(name = "Genes candidatos", range = c(4, 8),
                        breaks = sort(unique(ora_plot$Selected_Count))) +
  scale_x_continuous(limits = c(0, 105), breaks = seq(0, 100, 20)) +
  labs(title = "Enriquecimiento funcional de candidatos con mayor expresión",
       subtitle = "Análisis GO por sobrerrepresentación; procesos con FDR < 0.05",
       x = "Enriquecimiento (veces)", y = NULL) +
  theme_thesis(15)

# 15.4 Panel D. Genes del locus T6SSa en orden genómico

T6SSa_plot <- T6SSa_DE
T6SSa_plot$Tag <- paste0(T6SSa_plot$Label, " (", sub("CIW80_", "", T6SSa_plot$Locus_Tag), ")")
T6SSa_plot$Tag <- factor(T6SSa_plot$Tag, levels = rev(T6SSa_plot$Tag))

fig10D <- ggplot(T6SSa_plot, aes(logFC, Tag)) +
  geom_segment(aes(x = 0, xend = logFC, yend = Tag), color = "grey65", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(aes(size = logCPM), shape = 21, fill = "#D55E00", color = "black", stroke = 0.5) +
  geom_text(aes(label = sprintf("%.2f", logFC)), hjust = -0.5, size = 5) +
  scale_size_continuous(name = "Abundancia (logCPM)", range = c(4, 9)) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.17))) +
  labs(title = "Cambios de expresión de los genes del locus T6SSa",
       subtitle = "Orden genómico; ningún gen alcanzó FDR < 0.05",
       x = xlab_fc, y = "Gen o función anotada") +
  guides(size = guide_legend(nrow = 1)) +
  coord_cartesian(clip = "off") +
  theme_thesis(15)

save_figure(fig10A, "Figura_10A_Volcan", 17, 12)
save_figure(fig10B, "Figura_10B_Genes_GO", 15, 12)
save_figure(fig10C, "Figura_10C_Enriquecimiento_GO", 14, 8)
save_figure(fig10D, "Figura_10D_Expresion_T6SSa", 16, 14)


###############################################################################
# 16. FIGURAS 11 A 16. ORGANIZACIÓN Y CONTINUIDAD DE SST6a Y SST6c
###############################################################################

locus_name <- c(T6SSa = "SST6a", T6SSc = "SST6c")
subtitle_strains <- "Escherichia coli Nissle 1917 | WT frente a Ptet+"
boundary_axis <- function(x) gsub(" \u2192 ", " - ", x, fixed = TRUE)

prepare_boundaries <- function(data, locus) {
  d <- data[data$Locus == locus, ]
  d$Boundary  <- factor(d$Boundary_Label, levels = unique(d$Boundary_Label[order(d$Boundary_Index)]))
  d$Replicate <- factor(d$Label, levels = sample_levels)
  d$Strain    <- factor(d$Strain, levels = c("WT", "Ptet+"))
  d
}

heatmap_theme <- function(axis_size) {
  theme_thesis(14) +
    theme(axis.title.x = element_text(size = 15, margin = margin(t = 12)),
          axis.title.y = element_text(size = 15, margin = margin(r = 5)),
          axis.text.x = element_text(size = axis_size, angle = 35, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 14),
          panel.grid = element_blank())
}

# 16.1 Cobertura y arquitectura genómica (Figuras 11 y 14)

gene_arrows <- function(genes) {
  max_tip <- 0.18 * median(genes$End - genes$Start)
  do.call(rbind, lapply(seq_len(nrow(genes)), function(i) {
    g <- genes[i, ]
    tip <- min(max_tip, (g$End - g$Start) * 0.35)
    x <- if (g$Strand == "+") {
      c(g$Start, g$End - tip, g$End, g$End - tip, g$Start)
    } else {
      c(g$Start + tip, g$End, g$End, g$Start + tip, g$Start)
    }
    data.frame(x = x, y = c(-0.09, -0.09, 0, 0.09, 0.09), id = i, Category = g$Category)
  }))
}

coverage_figure <- function(genes, locus, label_size) {
  span <- range(c(genes$Start, genes$End))
  limits <- span + c(-1, 1) * diff(span) * 0.05

  cov <- coverage_profile[coverage_profile$Locus == locus &
                            coverage_profile$Position >= limits[1] &
                            coverage_profile$Position <= limits[2], ]
  cov$Strain <- factor(cov$Strain, levels = c("WT", "Ptet+"))

  p_cov <- ggplot(cov, aes(Position, CPM_Smoothed, color = Strain, linetype = Strain)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = col_strain) +
    scale_linetype_manual(values = c("WT" = "solid", "Ptet+" = "dashed")) +
    coord_cartesian(xlim = limits, expand = FALSE) +
    labs(x = NULL, y = "Cobertura normalizada (CPM)", color = "Condición", linetype = "Condición") +
    theme_thesis(14) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top")

  p_genes <- ggplot() +
    geom_polygon(data = gene_arrows(genes), aes(x, y, group = id, fill = Category),
                 color = "white", linewidth = 0.25) +
    geom_text(data = genes, aes(Midpoint, 0.22, label = Label),
              angle = 45, hjust = 0, size = label_size, check_overlap = TRUE) +
    geom_text(data = genes, aes(Midpoint, -0.22, label = Short_Locus),
              angle = 45, hjust = 1, size = label_size * 0.85) +
    scale_fill_manual(values = col_function) +
    scale_y_continuous(limits = c(-0.85, 0.95)) +
    coord_cartesian(xlim = limits, expand = FALSE, clip = "off") +
    labs(x = "Posición genómica (pb)", y = NULL, fill = "Categoría funcional") +
    theme_thesis(14) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank())

  (p_cov / p_genes) +
    plot_layout(heights = c(2.2, 1.8)) +
    plot_annotation(title = paste("Cobertura de RNA-seq y organización genómica de", locus_name[[locus]]),
                    subtitle = subtitle_strains,
                    theme = theme(plot.title = element_text(face = "bold", size = 17)))
}

# 16.2 Razón de continuidad (Figuras 12 y 15)

continuity_colors <- c("<0.75" = "#440154", "0.75\u20130.94" = "#3B528B",
                       "0.95\u20131.05" = "#21918C", "1.06\u20131.29" = "#5EC962",
                       "\u22651.30" = "#FDE725")

continuity_figure <- function(locus, text_size, axis_size) {
  d <- prepare_boundaries(continuity, locus)
  d$Interval <- cut(d$Continuity_Ratio, breaks = c(-Inf, 0.75, 0.95, 1.05, 1.30, Inf),
                    labels = names(continuity_colors), right = FALSE)

  ggplot(d, aes(Boundary, Replicate, fill = Interval)) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_text(aes(label = sprintf("%.2f", Continuity_Ratio),
                  color = Interval %in% names(continuity_colors)[1:2]), size = text_size) +
    geom_hline(yintercept = 3.5, color = "grey30", linewidth = 0.8) +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
    scale_fill_manual(values = continuity_colors, drop = FALSE, name = "Razón de continuidad",
                      guide = guide_legend(title.position = "top", nrow = 1)) +
    scale_y_discrete(limits = rev(sample_levels)) +
    scale_x_discrete(labels = boundary_axis) +
    labs(title = paste("Continuidad transcripcional intergénica de", locus_name[[locus]]),
         subtitle = subtitle_strains, x = "Frontera intergénica", y = "Réplica biológica") +
    heatmap_theme(axis_size)
}

# 16.3 Soporte paired-end (Figuras 13 y 16)

support_levels <- c("1\u20139", "10\u201399", "100\u2013999", "\u22651,000")

paired_figure <- function(locus, point_size, axis_size) {
  d <- prepare_boundaries(paired_support, locus)
  d$Support <- cut(d$Boundary_Supporting_Fragments, breaks = c(0, 9, 99, 999, Inf),
                   labels = support_levels)

  ggplot(d, aes(Boundary, Replicate)) +
    geom_point(data = d[d$Boundary_Supporting_Fragments > 0, ],
               aes(shape = Support, color = Strain), size = point_size, stroke = 1) +
    geom_text(data = d[d$Boundary_Supporting_Fragments == 0, ], label = "\u00d7",
              size = 4.5, color = "#333333", fontface = "bold") +
    geom_hline(yintercept = 3.5, color = "#555555", linewidth = 0.7) +
    scale_color_manual(values = c("WT" = "#4477AA", "Ptet+" = "#EE6677"), name = "Condición") +
    scale_shape_manual(values = setNames(c(16, 17, 15, 18), support_levels),
                       drop = FALSE, name = "Fragmentos de soporte") +
    scale_y_discrete(limits = rev(sample_levels)) +
    scale_x_discrete(labels = boundary_axis) +
    labs(title = paste("Soporte de fragmentos paired-end de", locus_name[[locus]]),
         subtitle = subtitle_strains, x = "Frontera intergénica", y = "Réplica biológica",
         caption = "\u00d7 = 0 fragmentos de soporte.") +
    guides(shape = guide_legend(order = 1, override.aes = list(color = "#333333", size = point_size)),
           color = guide_legend(order = 2, override.aes = list(shape = 16, size = 4.5))) +
    heatmap_theme(axis_size) +
    theme(panel.border = element_blank(),
          panel.grid.major = element_line(color = "grey93"),
          legend.box = "horizontal",
          plot.caption = element_text(hjust = 0.5, size = 9))
}

save_figure(coverage_figure(T6SSa_final, "T6SSa", 3.8), "Figura_11_Cobertura_SST6a", 20, 10)
save_figure(continuity_figure("T6SSa", 4.8, 12), "Figura_12_Continuidad_SST6a", 18, 8)
save_figure(paired_figure("T6SSa", 5, 12), "Figura_13_PairedEnd_SST6a", 19, 9)

save_figure(coverage_figure(T6SSc_core, "T6SSc", 4.5), "Figura_14_Cobertura_SST6c", 16, 9)
save_figure(continuity_figure("T6SSc", 5.5, 14), "Figura_15_Continuidad_SST6c", 15, 8)
save_figure(paired_figure("T6SSc", 5.5, 14), "Figura_16_PairedEnd_SST6c", 16, 9)


###############################################################################
# 17. FIGURAS 26 A 28 (ANEXOS)
###############################################################################

# 17.1 Figura 26. Robustez de la cobertura intergénica

robustness_figure <- function(locus, tag) {
  d <- prepare_boundaries(continuity, locus)

  ggplot(d, aes(Boundary, Replicate, fill = IGR_Percent_GE5)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.0f%%", IGR_Percent_GE5), color = IGR_Percent_GE5 >= 70), size = 3.1) +
    geom_hline(yintercept = 3.5, color = "grey35", linewidth = 0.7) +
    scale_fill_gradientn(colors = c("#F2F2F2", "#56B4E9", "#0072B2"), values = c(0, 0.5, 1),
                         limits = c(0, 100), breaks = seq(0, 100, 25),
                         name = "Bases con \u22655\nlecturas (%)") +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
    scale_y_discrete(limits = rev(sample_levels)) +
    scale_x_discrete(labels = boundary_axis) +
    labs(title = paste("Robustez de cobertura intergénica de", locus),
         subtitle = "Porcentaje de bases de cada región intergénica con al menos cinco lecturas",
         x = "Frontera intergénica", y = "Réplica biológica", tag = tag) +
    theme_thesis(11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
          panel.grid = element_blank(), legend.position = "right")
}

fig26 <- robustness_figure("T6SSa", "A") / robustness_figure("T6SSc", "B")
save_figure(fig26, "Figura_26_Robustez_Intergenica", 15, 11)

# 17.2 Figura 27. Enriquecimiento basado en ranking

gsea_plot <- GO_target_rank
gsea_plot$Process <- ifelse(gsea_plot$GO == GO_T6SS,
                            "Secreción de proteínas mediante T6SS", "Biosíntesis de tiamina")
gsea_plot$Process <- factor(gsea_plot$Process, levels = c("Biosíntesis de tiamina",
                                                          "Secreción de proteínas mediante T6SS"))
gsea_plot$Text <- sprintf("NES = %.2f | FDR = %.3f", gsea_plot$NES, gsea_plot$FDR)

fig27 <- ggplot(gsea_plot, aes(NES, Process)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(size = Set_Size), shape = 21, fill = "#D55E00", color = "black", stroke = 0.5) +
  geom_text(aes(label = Text), hjust = -0.15, size = 3.5) +
  scale_size_continuous(name = "Genes del conjunto", range = c(4, 8),
                        breaks = sort(unique(gsea_plot$Set_Size))) +
  scale_x_continuous(limits = c(0, 3.5), breaks = seq(0, 3.5, 0.5)) +
  labs(title = "Enriquecimiento basado en ranking de procesos GO",
       subtitle = "Procesos de interés; ninguno alcanzó FDR < 0.05",
       x = "Puntuación de enriquecimiento normalizada (NES)", y = NULL) +
  theme_thesis(11)

save_figure(fig27, "Figura_27_GSEA_T6SS_Tiamina", 11, 5)

# 17.3 Figura 28. PCA y coeficiente de variación biológica

fig28A <- ggplot(pca_df, aes(PC1, PC2, color = Strain, shape = Strain)) +
  geom_hline(yintercept = 0, color = "grey90") +
  geom_vline(xintercept = 0, color = "grey90") +
  geom_point(size = 3.6) +
  geom_text(aes(label = Label), hjust = -0.12, vjust = -0.7, size = 3, show.legend = FALSE) +
  scale_color_manual(values = col_strain) +
  scale_shape_manual(values = c("WT" = 16, "Ptet+" = 17)) +
  scale_x_continuous(expand = expansion(mult = 0.18)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  labs(title = "A. Análisis de componentes principales",
       subtitle = "Distribución de las seis bibliotecas de RNA-seq",
       x = sprintf("PC1 (%.2f %%)", pca_variance[1]),
       y = sprintf("PC2 (%.2f %%)", pca_variance[2]),
       color = "Condición", shape = "Condición") +
  theme_thesis(11)

bcv_data <- data.frame(logCPM = aveLogCPM(dge),
                       BCV = sqrt(dge$tagwise.dispersion),
                       Trend = sqrt(dge$trended.dispersion))
bcv_data <- bcv_data[order(bcv_data$logCPM), ]

fig28B <- ggplot(bcv_data, aes(logCPM, BCV)) +
  geom_point(color = "grey35", size = 0.4, alpha = 0.45) +
  geom_hline(yintercept = sqrt(dge$common.dispersion), color = "#D55E00", linewidth = 0.8) +
  geom_line(aes(y = Trend), color = "#0072B2", linewidth = 0.9) +
  labs(title = "B. Coeficiente de variación biológica",
       subtitle = "Variación estimada según la abundancia génica",
       x = "Abundancia media (logCPM)", y = "Coeficiente de variación biológica") +
  theme_thesis(11)

fig28 <- (fig28A | fig28B) +
  plot_layout(widths = c(1, 1.15), guides = "collect") &
  theme(legend.position = "bottom")

save_figure(fig28, "Figura_28_PCA_BCV", 16, 7.5)

cat("\nAnálisis completo. Resultados en:\n", output_dir, "\n")
