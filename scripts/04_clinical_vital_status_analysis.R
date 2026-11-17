################################################################################
# Project: TCGA-ESCA Clinical Metadata and Vital Status Analysis
# Purpose: Exploratory differential expression by patient vital status
# Data source: Genomic Data Commons (GDC) / TCGA-ESCA
#
# IMPORTANT:
# - The RNA-seq matrix currently contains only 8 samples (4 normal, 4 tumor).
# - Clinical records are patient-level; sample barcodes are matched using the
#   first 12 characters (TCGA patient/case identifier).
# - This is an exploratory demonstration, not a validated survival analysis.
# - Vital status (Alive/Dead) is NOT the same as tissue sample type (Normal/Tumor).
################################################################################

# 1. Load packages -------------------------------------------------------------
library(TCGAbiolinks)
library(DESeq2)

# 2. Set project directory -----------------------------------------------------
# Update this path for your computer.
setwd("D:/Genetics & Bioinformatics Researchs/Clinical case TCGA-ESCA")

# 3. Download and inspect clinical data ---------------------------------------
clinicaldata_esca <- GDCquery_clinic(
  project = "TCGA-ESCA",
  type = "clinical"
)

dim(clinicaldata_esca)
str(clinicaldata_esca)
table(clinicaldata_esca$vital_status, useNA = "ifany")

write.table(
  clinicaldata_esca,
  file = "ESCA_Clinical_Data.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

# 4. Prepare patient-level vital-status metadata -------------------------------
clinical_vital <- clinicaldata_esca[
  !is.na(clinicaldata_esca$vital_status) &
    clinicaldata_esca$vital_status %in% c("Alive", "Dead"),
  c("submitter_id", "vital_status")
]

colnames(clinical_vital) <- c("PatientID", "VitalStatus")
clinical_vital$PatientID <- as.character(clinical_vital$PatientID)
clinical_vital$VitalStatus <- as.character(clinical_vital$VitalStatus)

# Ensure one clinical record per patient before matching.
clinical_vital <- clinical_vital[!duplicated(clinical_vital$PatientID), ]

if (anyDuplicated(clinical_vital$PatientID)) {
  stop("Duplicate patient IDs remain in the clinical metadata.")
}

table(clinical_vital$VitalStatus)

# 5. Import raw RNA-seq counts --------------------------------------------------
# Keep original TCGA barcodes: check.names = FALSE prevents R from replacing
# hyphens with dots. The first column contains Ensembl gene IDs.
count_matrix <- read.delim(
  "8 Samples Raw data.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

count_matrix <- as.matrix(count_matrix)
storage.mode(count_matrix) <- "numeric"

# Validate raw count data.
if (anyNA(count_matrix) ||
    any(!is.finite(count_matrix)) ||
    any(count_matrix < 0) ||
    any(count_matrix != floor(count_matrix))) {
  stop("The expression matrix must contain non-negative integer raw counts.")
}

storage.mode(count_matrix) <- "integer"

dim(count_matrix)
colnames(count_matrix)[1: ncol(count_matrix)]

# 6. Match RNA-seq samples to patient-level clinical records -------------------
# TCGA patient identifier = first 12 characters of the full sample barcode.
sample_barcodes <- colnames(count_matrix)
expression_patient_ids <- substr(sample_barcodes, 1, 12)

# Check whether the sample barcode format is compatible with clinical IDs.
if (any(grepl("\\.", sample_barcodes))) {
  stop(
    "Sample barcodes contain dots. Re-import the raw matrix with ",
    "check.names = FALSE to preserve the original TCGA barcodes."
  )
}

clinical_match <- match(
  expression_patient_ids,
  clinical_vital$PatientID
)

matched_samples <- !is.na(clinical_match)

if (sum(matched_samples) == 0) {
  stop("No RNA-seq samples matched the clinical patient IDs.")
}

count_matrix_selected <- count_matrix[
  ,
  matched_samples,
  drop = FALSE
]

clinical_selected <- clinical_vital[
  clinical_match[matched_samples],
  ,
  drop = FALSE
]

# Use full RNA-seq sample barcodes as colData row names.
rownames(clinical_selected) <- colnames(count_matrix_selected)

clinical_selected$VitalStatus <- factor(
  clinical_selected$VitalStatus,
  levels = c("Alive", "Dead")
)

# Confirm exact sample order.
stopifnot(
  identical(
    colnames(count_matrix_selected),
    rownames(clinical_selected)
  )
)

cat("Matched RNA-seq samples:", ncol(count_matrix_selected), "\n")
print(table(clinical_selected$VitalStatus))

write.table(
  clinical_selected,
  file = "ESCA_Matched_Clinical_Metadata.txt",
  quote = FALSE,
  sep = "\t",
  col.names = NA
)

# 7. Differential expression by vital status ----------------------------------
# This comparison is exploratory and may be underpowered with this small subset.
status_counts <- table(clinical_selected$VitalStatus)

if (length(status_counts) < 2 || any(status_counts < 2)) {
  stop(
    "At least two samples per vital-status group are required for this ",
    "exploratory DESeq2 workflow. Review the matched sample counts."
  )
}

# Create DESeq2 dataset.
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_selected,
  colData = clinical_selected,
  design = ~ VitalStatus
)

# Remove very low-count genes.
keep_genes <- rowSums(counts(dds)) >= 10
dds <- dds[keep_genes, ]

# Run DESeq2.
dds <- DESeq(dds)

# Explicitly compare Dead against Alive.
deseq_results <- results(
  dds,
  contrast = c("VitalStatus", "Dead", "Alive"),
  pAdjustMethod = "BH"
)

# 8. Process and export results -------------------------------------------------
all_results <- as.data.frame(deseq_results)
all_results$GeneID <- rownames(all_results)

all_results <- all_results[
  order(all_results$padj, na.last = TRUE),
]

write.table(
  all_results,
  file = "ESCA_Alive_vs_Dead_All_Results.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

significant_results <- subset(
  all_results,
  !is.na(padj) & padj < 0.05
)

write.table(
  significant_results,
  file = "ESCA_Alive_vs_Dead_Significant_DEGs.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

# 9. Add gene-symbol annotation ------------------------------------------------
# Expected annotation file: Annot_df.txt
# Expected columns after import: Symbol and BioType; row names = Ensembl IDs.
if (file.exists("Annot_df.txt")) {
  annotation <- read.delim(
    "Annot_df.txt",
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE
  )

  annotation_index <- match(
    significant_results$GeneID,
    rownames(annotation)
  )

  significant_results$Symbol <- annotation$Symbol[annotation_index]
  significant_results$BioType <- annotation$BioType[annotation_index]
} else {
  warning(
    "Annot_df.txt was not found. Significant results are exported ",
    "without gene-symbol annotation."
  )
}

write.table(
  significant_results,
  file = "ESCA_Alive_vs_Dead_Annotated_DEGs.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

# 10. Final checks --------------------------------------------------------------
cat("Total genes tested:", nrow(all_results), "\n")
cat("Significant genes (FDR < 0.05):", nrow(significant_results), "\n")
cat("Analysis complete. Interpret all findings as exploratory.\n")

################################################################################
# END
################################################################################
