# =========================================================================================
# Description: Multi-omic single-cell analysis pipeline (Checkpoints, Metabolism, Interactome, TFs)
# Cohort: GSE131907 (Human Early-Stage Lung Adenocarcinoma)
# =========================================================================================

set.seed(42)

# ---------------------------------------------------------
# Step 1: Environment Setup & Initialization
# ---------------------------------------------------------
cat("\n[1] Initializing Environment and Loading Libraries...\n")
suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    library(Matrix)
    library(ggplot2)
    library(patchwork)
    library(broom)
})

# Create directory for outputs
dir.create("meta_data", showWarnings = FALSE)
options(timeout = max(3600, getOption("timeout")))

# ---------------------------------------------------------
# Step 2: Download & Decompress Data (GSE131907)
# ---------------------------------------------------------
cat("\n[2] Preparing Dataset...\n")
url_meta <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
url_matrix <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.rds.gz"

dest_meta <- "meta_data/GSE131907_cell_annotation.txt.gz"
dest_matrix <- "meta_data/GSE131907_raw_UMI_matrix.rds.gz"
dest_unzipped <- "meta_data/GSE131907_raw_UMI_matrix.rds"

if(!file.exists(dest_meta)) download.file(url_meta, dest_meta, mode = "wb")
if(!file.exists(dest_matrix) && !file.exists(dest_unzipped)) download.file(url_matrix, dest_matrix, mode = "wb")

if (!file.exists(dest_unzipped)) {
    cat("Unzipping the doubly-compressed RDS file (this takes a moment)...\n")
    if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils", repos = "http://cran.us.r-project.org")
    R.utils::gunzip(dest_matrix, destname = dest_unzipped, remove = FALSE)
}

# ---------------------------------------------------------
# Step 3: Memory-Safe Subsetting & Seurat Object Creation
# ---------------------------------------------------------
cat("\n[3] Memory-Safe Data Loading & Subsetting...\n")
human_meta <- read_tsv(dest_meta, show_col_types = FALSE) %>% as.data.frame()
rownames(human_meta) <- human_meta[[1]]                                        

annot_cols <- grep("type|cell|class|cluster", colnames(human_meta), ignore.case = TRUE, value = TRUE)
main_annot <- annot_cols[1]

human_counts <- readRDS(dest_unzipped)
total_cells <- ncol(human_counts)

# Engage Memory-Safe Subsetting (Cap at 35,000 cells)
immune_idx <- grep("T cell|Macrophage|Myeloid|NK|B cell|Lympho|Dendritic|Monocyte|Neutrophil|Mast|CD4|CD8", human_meta[[main_annot]], ignore.case = TRUE)
epi_idx <- grep("Epithelial|Tumor|Malignant", human_meta[[main_annot]], ignore.case = TRUE)

keep_immune <- immune_idx[sample(length(immune_idx), min(25000, length(immune_idx)))]
keep_epi <- epi_idx[sample(length(epi_idx), min(10000, length(epi_idx)))]
keep_idx <- c(keep_immune, keep_epi)
keep_barcodes <- rownames(human_meta)[keep_idx]

# Subset dense dataframe BEFORE coercion to avoid RAM crash
human_counts <- human_counts[, keep_barcodes]
human_meta <- human_meta[keep_barcodes, ]

cat("Converting to Sparse Matrix and Building Seurat Object...\n")
human_counts <- Matrix::Matrix(as.matrix(human_counts), sparse = TRUE)
seu_human <- CreateSeuratObject(counts = human_counts, project = "Human_EarlyLUAD", meta.data = human_meta)
seu_human[["percent.mt"]] <- PercentageFeatureSet(seu_human, pattern = "^MT-")

rm(human_counts)
gc()

# ---------------------------------------------------------
# Step 4: Normalization & Surface Checkpoint Scoring
# ---------------------------------------------------------
cat("\n[4] Normalizing and Scoring Checkpoint Modules...\n")
seu_human <- NormalizeData(seu_human, verbose = FALSE)
seu_human$Cell_Type <- seu_human$Cell_type

all_human_genes <- rownames(seu_human)
human_modules <- list(
    Innate_Suppression = c("HAVCR2", "TGFB1", "IL10", "ARG1", "LGALS9"), 
    Ceacam1_Module = c("CEACAM1"),
    Adaptive_Evasion = c("CD274", "CD276", "CTLA4", "PDCD1"), 
    Cytotoxic = c("GZMB", "PRF1", "IFNG")
)

match_human_genes <- function(genes) {
    matched <- c()
    for (g in genes) {
        m <- grep(paste0("^", g, "$"), all_human_genes, ignore.case = TRUE, value = TRUE)
        if (length(m) > 0) matched <- c(matched, m[1])
    }
    return(matched)
}

valid_human_modules <- lapply(human_modules, match_human_genes)
if(length(valid_human_modules$Innate_Suppression)>0) seu_human <- AddModuleScore(seu_human, features=list(valid_human_modules$Innate_Suppression), name="Mod_Innate")
if(length(valid_human_modules$Adaptive_Evasion)>0) seu_human <- AddModuleScore(seu_human, features=list(valid_human_modules$Adaptive_Evasion), name="Mod_Adaptive")
if(length(valid_human_modules$Cytotoxic)>0) seu_human <- AddModuleScore(seu_human, features=list(valid_human_modules$Cytotoxic), name="Mod_Cytotoxic")

# ---------------------------------------------------------
# Step 5: Checkpoint Landscape Plotting
# ---------------------------------------------------------
cat("\n[5] Generating Checkpoint Landscape Plots...\n")
seu_human@meta.data$Cell_Type <- factor(seu_human@meta.data$Cell_Type, levels = c("Myeloid cells", "MAST cells", "NK cells", "T lymphocytes", "B lymphocytes", "Epithelial cells"))
seu_plot <- subset(seu_human, subset = Cell_Type %in% levels(seu_human@meta.data$Cell_Type))

p1 <- VlnPlot(seu_plot, features = "Mod_Innate1", group.by = "Cell_Type", pt.size = 0) + geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) + theme_classic(base_size = 12) + labs(title = "Human Innate Suppression (HAVCR2/TIM-3+)", x = "", y = "Module Score") + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_brewer(palette = "Set2")
p2 <- VlnPlot(seu_plot, features = "Mod_Adaptive1", group.by = "Cell_Type", pt.size = 0) + geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) + theme_classic(base_size = 12) + labs(title = "Human Adaptive Checkpoint (PDCD1/PD-1+)", x = "", y = "Module Score") + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_brewer(palette = "Set1")

ggsave("meta_data/Human_Checkpoint_Compartmentalization.pdf", plot = p1 | p2, width = 12, height = 6, dpi = 300)
ggsave("meta_data/Human_Checkpoint_Compartmentalization.png", plot = p1 | p2, width = 12, height = 6, dpi = 300)

# ---------------------------------------------------------
# Step 6: Immuno-Metabolic Crosstalk & Correlation
# ---------------------------------------------------------
cat("\n[6] Scoring Metabolic Modules & Testing Correlation...\n")
metabolic_modules <- list(
    Glycolysis = c("ENO1", "HK2", "ALDOA", "GAPDH", "LDHA", "PKM"),
    OXPHOS = c("COX5B", "ATP5B", "NDUFB8", "SDHB", "UQCRC2")
)
valid_metabolic <- lapply(metabolic_modules, match_human_genes)

if(length(valid_metabolic$Glycolysis)>0) seu_human <- AddModuleScore(seu_human, features=list(valid_metabolic$Glycolysis), name="Met_Glycolysis")
if(length(valid_metabolic$OXPHOS)>0) seu_human <- AddModuleScore(seu_human, features=list(valid_metabolic$OXPHOS), name="Met_OXPHOS")

# Statistical Correlation (T cells: PD1 vs Glycolysis)
t_cells <- seu_human@meta.data %>% filter(Cell_Type == "T lymphocytes")
res_glycolysis <- cor.test(t_cells$Mod_Adaptive1, t_cells$Met_Glycolysis1, method = "pearson")
cat(sprintf("Pearson Correlation (PD1 vs Glycolysis): R = %.3f, P-value = %.2e\n", res_glycolysis$estimate, res_glycolysis$p.value))

# ---------------------------------------------------------
# Step 7: Paracrine Ligand-Receptor Interactome (Outsourcing)
# ---------------------------------------------------------
cat("\n[7] Modeling Paracrine Ligand-Receptor Networks...\n")
lr_pairs <- data.frame(
    Pathway = c("PD1_Axis", "TIM3_Axis", "CTLA4_Axis"),
    Ligand = c("CD274", "LGALS9", "CD86"),  
    Receptor = c("PDCD1", "HAVCR2", "CTLA4"), 
    stringsAsFactors = FALSE
) %>% filter(Ligand %in% all_human_genes & Receptor %in% all_human_genes)

target_cells <- subset(seu_human, subset = Cell_Type %in% c("Myeloid cells", "T lymphocytes", "Epithelial cells"))
expr_data <- AverageExpression(target_cells, features = unique(c(lr_pairs$Ligand, lr_pairs$Receptor)), group.by = "Cell_Type")$RNA %>% as.data.frame() %>% rownames_to_column("Gene")

crosstalk_results <- lr_pairs %>% rowwise() %>% mutate(
    Myeloid_to_T_Score = expr_data[expr_data$Gene == Ligand, "Myeloid cells"] * expr_data[expr_data$Gene == Receptor, "T lymphocytes"],
    Epi_to_T_Score = expr_data[expr_data$Gene == Ligand, "Epithelial cells"] * expr_data[expr_data$Gene == Receptor, "T lymphocytes"]
) %>% arrange(desc(Myeloid_to_T_Score))

write_tsv(crosstalk_results, "meta_data/LR_Crosstalk_Results.tsv")
print(crosstalk_results)

# ---------------------------------------------------------
# Step 8: Transcription Factor (TF) Profiling & Final Save
# ---------------------------------------------------------
cat("\n[8] Scoring TF Modules & Saving Final Object...\n")
tf_modules <- list(
    TF_Exhaustion = c("TOX", "PRDM1", "BATF", "IRF4"),
    TF_Suppressive_Myeloid = c("STAT3", "CEBPB", "HIF1A", "PPARG"),
    TF_Stemness = c("TCF7", "LEF1", "SELL")
)
valid_tfs <- lapply(tf_modules, match_human_genes)

if(length(valid_tfs$TF_Exhaustion)>0) seu_human <- AddModuleScore(seu_human, features=list(valid_tfs$TF_Exhaustion), name="TF_Exhaustion")
if(length(valid_tfs$TF_Suppressive_Myeloid)>0) seu_human <- AddModuleScore(seu_human, features=list(valid_tfs$TF_Suppressive_Myeloid), name="TF_Supp_Myeloid")
if(length(valid_tfs$TF_Stemness)>0) seu_human <- AddModuleScore(seu_human, features=list(valid_tfs$TF_Stemness), name="TF_Stemness")

saveRDS(seu_human, "meta_data/seu_human_final_master.rds")
cat("\nPipeline Complete! Master Human Object saved to 'meta_data/seu_human_final_master.rds'.\n")
