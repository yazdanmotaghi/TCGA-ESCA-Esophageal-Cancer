################################################################################
# Project: Gene Expression Analysis Using TCGA-ESCA Data
# Research Area: Cancer Genomics, Transcriptomics & Bioinformatics
# Researcher: Yazdan Mottaghi
# Data Source: Genomic Data Commons (GDC)
# Project ID: TCGA-ESCA
################################################################################


# 1. Load Required Packages
library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(readxl)

# 2. Set Working Directory
setwd("D:/Genetics & Bioinformatics Researchs/TCGA-ESCA")


################################################################################
# 3. Explore Available GDC Projects
################################################################################

All_GDC_Projects <- getGDCprojects()

"TCGA-ESCA" %in% All_GDC_Projects$project_id

View(All_GDC_Projects)
All_GDC_Projects[, 1]


################################################################################
# 4. Explore TCGA-BRCA Project Summary
################################################################################

All_TCGA_ESCA_Data <- getProjectSummary("TCGA-ESCA")
View(All_TCGA_ESCA_Data)
All_TCGA_ESCA_Data$data_categories


################################################################################
# 5. Query TCGA-BRCA Gene Expression Data
################################################################################

query_TCGA_ESCA <- GDCquery(
  project = "TCGA-ESCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Solid Tissue Normal","Primary Tumor")
)


################################################################################
# 6. Retrieve Query Results
################################################################################

results <- getResults(query_TCGA_ESCA)

View(results)


################################################################################
# 7. Identify Normal Tissue Samples
################################################################################

Normal_sample_list <- subset(
  results,
  sample_type == "Solid Tissue Normal" &
    access == "open" &
    experimental_strategy == "RNA-Seq"
)

View(Normal_sample_list)

write.table(
  Normal_sample_list,
  file = "all normal samples of ESCA.txt",
  quote = FALSE,
  sep = "\t"
)

nrow(Normal_sample_list)


################################################################################
# 8. Identify Primary Tumor Samples
################################################################################

Tumor_sample_list <- subset(
  results,
  sample_type == "Primary Tumor" &
    access == "open" &
    experimental_strategy == "RNA-Seq"
)

View(Tumor_sample_list)

write.table(
  Tumor_sample_list,
  file = "all tumor samples of ESCA.txt",
  quote = FALSE,
  sep = "\t"
)

nrow(Tumor_sample_list)

length(Tumor_sample_list)
length(Normal_sample_list)


################################################################################
# 9. Select Four Tumor and Four Normal Samples
################################################################################

sample_T_4 <- Tumor_sample_list[1:4, ]
sample_N_4 <- Normal_sample_list[1:4, ]

# Keep the same sample order used in the original selection.
Selected_samples <- rbind(sample_N_4, sample_T_4)

View(Selected_samples)


################################################################################
# 10. Export Selected Samples
################################################################################

write.table(
  Selected_samples,
  file = "selected samples ESCA.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)


################################################################################
# 11. Create a GDC Query for Selected Samples
################################################################################

query_TCGA_2 <- GDCquery(
  project = "TCGA-ESCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = Selected_samples$sample.submitter_id
)


################################################################################
# 12. Download RNA-seq Data
################################################################################

GDCdownload(
  query_TCGA_2,
  method = "api",
  directory = "D:/Genetics & Bioinformatics Researchs/TCGA-ESCA",
  files.per.chunk = 8
)


################################################################################
# 13. Prepare RNA-seq Data
################################################################################

data_prepare <- GDCprepare(
  query_TCGA_2,
  directory ="D:/Genetics & Bioinformatics Researchs/TCGA-ESCA",
  save = TRUE,
  save.filename = "8 samples data.txt"
)


################################################################################
# 14. Inspect the Prepared Dataset
################################################################################

class(data_prepare)
dim(data_prepare)


################################################################################
# 15. Extract Raw Gene Expression Counts
################################################################################
library(SummarizedExperiment)
data <- assay(data_prepare)

View(data)
dim(data)


################################################################################
# 16. Export Raw Count Matrix
################################################################################

write.table(
  data,
  file = "8 Samples Raw data.txt",
  quote = FALSE,
  sep = "\t",
  col.names = NA
)


################################################################################
# 17. Create Gene Annotation File
################################################################################
# Load Package
library(readxl)

# Annotation
gene_annotation <- read_excel("Annotation.xlsx",
  col_names = FALSE)

library(TCGAbiolinks)
dim(gene_annotation)
head(gene_annotation)
View(gene_annotation)
nrow(gene_annotation)
nrow(data)

Annot_df <- gene_annotation[,1:3]
View(Annot_df)
rownames(Annot_df)<- rownames(data)
Annot_df<-Annot_df[,-1]
colnames(Annot_df) <- c("Symbol","BioType")
write.table(Annot_df,"Annot_df.txt",quote = F,sep = "\t")

#Examples
Annot_df["ENSG00000000971.16",] 
Annot_df["ENSG00000006128.12",2] 


################################################################################
# 18. Explore and Filter Raw Counts
################################################################################

class(data)
library(SummarizedExperiment)
summary(as.vector(data))

colnames(data)

boxplot(data)
hist(data)
hist(data, breaks = 1000000, xlim = c(0, 50))

data_filtered <- data[which(rowSums(data) > 10), ]

nrow(data)
nrow(data_filtered)

write.table(
  data_filtered,
  file = "data filtered.txt",
  quote = FALSE,
  sep = "\t",
  col.names = NA
)


################################################################################
# 19. Read Filtered Data and Annotation
################################################################################

Annot_df <- read.table(
  "Annot_df.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

data_filtered <- read.table(
  "data filtered.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

View(Annot_df)
View(data_filtered)


################################################################################
# 20. Add Gene Annotation to Filtered Data
################################################################################

gene_name <- Annot_df[rownames(data_filtered), "Symbol"]
gene_type <- Annot_df[rownames(data_filtered), "BioType"]

data_filtered2 <- data_filtered

data_filtered2 <- cbind(data_filtered2, gene_name)
data_filtered2 <- cbind(data_filtered2, gene_type)

data_filtered2 <- as.data.frame(data_filtered2)

class(data_filtered2)
table(data_filtered2$gene_type)

View(data_filtered2)


################################################################################
# 21. Separate RNA Types
################################################################################

mRNA_df <- subset(
  data_filtered2,
  gene_type == "protein_coding"
)

View(mRNA_df)
nrow(mRNA_df)

lncRNA_df <- subset(
  data_filtered2,
  gene_type == "lncRNA"
)

View(lncRNA_df)
nrow(lncRNA_df)

miRNA_df <- subset(
  data_filtered2,
  gene_type == "miRNA"
)

View(miRNA_df)
nrow(miRNA_df)


################################################################################
# 23. Export Expression Data
################################################################################


write.table(
  lncRNA_df,
  file = "lncRNA_df.txt",
  quote = FALSE,
  sep = "\t"
)

write.table(
  miRNA_df,
  file = "miRNA_df.txt",
  quote = FALSE,
  sep = "\t"
)

write.table(
  mRNA_df,
  file = "mRNA_df.txt",
  quote = FALSE,
  sep = "\t"
)


################################################################################
# 24. Plot RNA Type Distribution
################################################################################

RNA_table <- table(data_filtered2$gene_type)
RNA_table <- as.data.frame(RNA_table)

colnames(RNA_table) <- c("Type", "Number")

View(RNA_table)

library(ggplot2)

Plot_RNA_table <- ggplot(
  RNA_table,
  aes(x = Type, y = Number)
) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(
    title = "Barplot",
    x = "RNA Type",
    y = "Number"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  "RNA table info.png",
  plot = Plot_RNA_table,
  width = 15,
  height = 7
)


################################################################################
# 25. Create Sample Metadata
################################################################################
#metadata

Selected_samples

read.table("selected samples ESCA.txt", header = T, sep = "\t")
View(Selected_samples)
Samples_type <- c(rep("Normal",4), rep("Tumor",4))

Meta_data <- data.frame(Selected_samples,Samples_type)

View(Meta_data)
head(Meta_data)
colnames(Meta_data)

Meta_data <- data.frame(
  Sample_ID = Selected_samples$sample.submitter_id,
  Samples_type = Samples_type
)

View(Meta_data)
head(Meta_data)

colnames(Meta_data) <- c("SampleName","SampleType")

# Check sample order with data
Selected_samples$sample.submitter_id == Meta_data[,1]

data <- read.table("8 Samples Raw data.txt",sep = "\t")
rownames(data) <- data$V1
data <- data[,-1]
colnames(data) <- as.character(unlist(data[1, ]))
data <- data[-1, ]
View(data)
colnames(data)
Meta_data[,1]

# Sample ID
Sample_ID <- substr(colnames(data), 1, 16)

# Check Sample ID
all(Sample_ID == Meta_data[,1])
all(Sample_ID == Meta_data$SampleName)

write.table(Meta_data,"Metadata.txt", quote = F,sep = "\t",row.names = F)

#######################################################################################
# DEG

# Read mRNA data
mRNA_df <- read.table(
  "mRNA_df.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

View(mRNA_df)

# Extract raw counts
mRNA_Raw <- mRNA_df[, -9:-10]
View(mRNA_Raw)

# Convert counts to numeric
mRNA_Raw[,1]
mRNA_Raw[,1] <- as.numeric(mRNA_Raw[,1])

for (i in 1:ncol(mRNA_Raw)) {
  mRNA_Raw[,i] <- as.numeric(mRNA_Raw[,i])
}


# Read sample metadata
Meta_data <- read.table(
  "Metadata.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

View(Meta_data)


# DESeq2
library(DESeq2)

# Create DESeqDataSet for all mRNA genes
d1 <- DESeqDataSetFromMatrix(
  countData = mRNA_Raw,
  colData = Meta_data,
  design = ~SampleType
)

dim(d1)
colData(d1)
design(d1)


# rlog Transformation
data_rlogd1 <- rlogTransformation(d1)

library(SummarizedExperiment)

data_rlogd1_df <- data.frame(assay(data_rlogd1))

# Export rlog-transformed data
write.table(
  data_rlogd1_df,
  "d1 rlog.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

# Differential Gene Expression Analysis (DEG)

library(DESeq2)

?results

d2 <- estimateSizeFactors(d1)
d3 <- DESeq(d2)

deseq_result <- results(d3, pAdjustMethod = "fdr")

class(deseq_result)

final_deseq_results <- data.frame(deseq_result)

View(final_deseq_results)

write.table(
  final_deseq_results,
  "DEG mRNA.txt",
  quote = F,
  sep = "\t"
)


DEG <- read.table("DEG mRNA.txt", sep = "\t")
View(DEG)

Annot_df <- read.table("Annot_df.txt", sep = "\t")
View(Annot_df)

# Annotation Function

annot_fx <- function(x) {
  k <- Annot_df[x, 1]
  return(k)
}

annot_fx("ENSG00000002834.1") # Example

gene_names <- annot_fx(rownames(DEG))

length(gene_names) == nrow(DEG)

final_result_annot <- cbind(gene_names)

final_result_annot <- cbind(
  final_result_annot,
  DEG[1:6]
)
View(final_result_annot)

rownames(final_result_annot) <- rownames(DEG)

colnames(final_result_annot)[1] <- "Symbol"

View(final_result_annot)

# Sort by padj value

final_result_annot <- final_result_annot[
  order(final_result_annot$padj, decreasing = F),
]

View(final_result_annot)

# Export Final Results

write.table(
  final_result_annot,
  "Final result.txt",
  quote = F,
  sep = "\t"
)

##########################################################################
# Read Raw and Rlog Data
##########################################################################

# Read rlog-transformed data
data_rlog_df <- read.table(
  "d1 rlog.txt",
  header = TRUE,
  sep = "\t"
)

View(data_rlog_df)
##########################################################################
# Column Sum Visualization
##########################################################################

test_df <- data_rlog_df

# Calculate column sums
colsums_list <- colSums(test_df)

# Extract column names
name_list <- colnames(test_df)

# Combine results
colsums_list <- c()
name_list <- c()

for (k in 1:ncol(test_df)) {
  
  colsums_list <- append(
    colsums_list,
    sum(test_df[, k], na.rm = TRUE)
  )
  
  name_list <- append(
    name_list,
    colnames(test_df)[k]
  )
}

sum_df <- data.frame(name_list, colsums_list)

View(sum_df)

# Bar Plot
library(ggplot2)

plot <- ggplot(
  sum_df,
  aes(x = name_list, y = colsums_list)
) +
  geom_bar(stat = "identity", fill = "tomato") +
  labs(
    title = "Bar plot",
    x = "Samples",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

plot

ggsave(
  "rlog d1 plot.png",
  plot = plot,
  width = 12,
  height = 10
)


##########################################################################
# Data Visualization
# Volcano Plot
##########################################################################

library(EnhancedVolcano)

packageVersion("EnhancedVolcano")

final_result_annot <- read.table(
  "Final result.txt",
  sep = "\t",
  header = T
)

View(final_result_annot)

colnames(final_result_annot)

top_genes <- final_result_annot[1:10, 1]

View(top_genes)

my_plot <- EnhancedVolcano(
  final_result_annot,
  title = "Genes volcano plot",
  lab = final_result_annot$Symbol,
  x = "log2FoldChange",
  y = "pvalue",
  selectLab = top_genes,
  xlab = bquote(~log[2]~'Fold change'),
  pCutoff = 0.01,
  FCcutoff = 2.0,
  pointSize = 1.0,
  labSize = 4.0,
  labCol = "black",
  labFace = "bold",
  boxedLabels = TRUE,
  colAlpha = 4/5,
  legendPosition = "right",
  legendLabSize = 14,
  legendIconSize = 4.0,
  drawConnectors = F,
  widthConnectors = 1.0,
  colConnectors = "black"
)

my_plot

ggsave(
  "genes volcanoplot.png",
  plot = my_plot,
  width = 12,
  height = 10
)

##########################################################################
# Sort by padj value
final_result_annot <- final_result_annot[
  order(final_result_annot$padj, decreasing = F),
]

View(final_result_annot)
##########################################################################
# Heatmap
##########################################################################

library(pheatmap)

# Keep original rlog data
norm_data <- data_rlog_df
View(norm_data)

# Row names of normalized data
head(rownames(norm_data))

# Row names of annotation results
head(rownames(final_result_annot))

# Gene symbols
head(final_result_annot$Symbol)

# Check matching
sum(rownames(final_result_annot) %in% rownames(norm_data))

all(rownames(norm_data) == rownames(final_result_annot))

# Assign gene symbols
rownames(norm_data) <- final_result_annot$Symbol

norm_data2 <- as.data.frame(norm_data)

# Make gene symbols unique
rownames(norm_data2) <- make.unique(final_result_annot$Symbol)

# Select top 100 genes
top_genes_count <- norm_data2[
  make.unique(final_result_annot$Symbol[1:100]),
  ,
  drop = FALSE
]

dim(top_genes_count)
head(top_genes_count)

View(top_genes_count)

# Rename sample columns
colnames(top_genes_count) <- c(
  "Normal-1", "Normal-2", "Normal-3", "Normal-4",
  "Tumor-1", "Tumor-2", "Tumor-3", "Tumor-4"
)

# Heatmap
library(gplots)
library(pheatmap)

Heat_plot <- pheatmap(
  top_genes_count,
  color = greenred(256),
  fontsize = 8,
  cellwidth = 160,
  cellheight = 12,
  angle_col = 45,
  main = "Top Genes",
  fontsize_row = 12,
  fontsize_col = 15
)

png(
  "Heatmap.png",
  width = 1500,
  height = 1500
)

Heat_plot

dev.off()

##########################################################################
# PCA Plot
##########################################################################

library(ggplot2)

Meta_data

norm_data <- data_rlog_df

pc <- prcomp(t(norm_data))

pcr <- data.frame(
  pc$x[, 1:3],
  Meta_data$SampleType
)

Group <- Meta_data$SampleType

PCA_plot <- ggplot(
  pcr,
  aes(PC1, PC2, color = Group)
) +
  geom_point(size = 3)

PCA_plot

ggsave(
  "PCA plot.png",
  plot = PCA_plot,
  width = 7,
  height = 6
)


##########################################################################
# End of Visualization
##########################################################################