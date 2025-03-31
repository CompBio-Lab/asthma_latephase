library(anndata)
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(muscat)
library(scater)
library(limma)
library(ComplexHeatmap)

options(future.globals.maxSize = 16000 * 1024^2)
scperturb <- read_h5ad(here::here("data", "preprocessedData", "GSE279945_sc_counts_processed.h5ad"))
scperturb <- CreateSeuratObject(counts = t(as.matrix(scperturb$X)), 
                           meta.data = scperturb$obs)

# create singleCellExperiment object 
sce <- SingleCellExperiment(
  assays = list(counts = scperturb@assays$RNA$counts),
  colData = scperturb@meta.data
)

# remove undetected genes
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

# calculate per-cell quality control (QC) metrics
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes 
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

# compute sum-factors & normalize
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)

metadata <- colData(sce)
df <- t(sce@assays@data$logcounts) %>% 
  as.data.frame() %>% 
  mutate(cell_type_orig = metadata$cell_type_orig,
         sm_name = metadata$sm_name,
         dose_uM = metadata$dose_uM,
         donor_id = metadata$donor_id,
         plate_name = metadata$plate_name,
         library_id = metadata$library_id) %>% 
  dplyr::group_by(cell_type_orig, sm_name, dose_uM, donor_id,
                  plate_name, library_id) %>% 
  summarise_each(funs(mean))

df$sm_name <- relevel(df$sm_name, ref = "Dimethyl Sulfoxide")
latephase_exacer_biomarkers <- readRDS(here::here("results", "latephase_exacer_biomarkers.rds"))
result_drug <- df %>% 
  group_by(cell_type_orig) %>% 
  nest() %>% 
  mutate(result = purrr::map(data, ~{
    block <- .$donor_id
    design <- model.matrix(~sm_name, data = .)
    eset <- t(.[, -c(1:5)])
    v <- voom(eset, design, plot=FALSE)
    dupcor <- duplicateCorrelation(v, design, block = block)
    vobj <- voom(eset, design, plot = FALSE, block = block, correlation = dupcor$consensus)
    # fit <- lmFit(eset, design)
    # fit <- eBayes(fit, trend = TRUE)
    fit <- lmFit(vobj, design, block = block, correlation = dupcor$consensus)
    top_all <- lapply(colnames(design)[2:144], function(i){
      top <- topTable(fit, coef = i, adjust.method = "BH", n = nrow(fit)) %>% 
        mutate(gene = rownames(.)) %>% 
        filter(gene %in% latephase_exacer_biomarkers) %>% 
        mutate(adj.P.Val = p.adjust(P.Value, "BH"))
      top$sm_name <- gsub("sm_name", "", i)
      top
    }) %>% 
      do.call(rbind, .)
    top_all
  }))
result_drug2 <- result_drug %>% 
  dplyr::select(-data) %>% 
  unnest(result)

sig_drugs <- result_drug2 %>% 
  mutate(FC = ifelse(logFC > 0, "UP", "DOWN")) %>% 
  dplyr::group_by(sm_name, cell_type_orig, FC) %>% 
  filter(adj.P.Val < 0.1)

sig_drugs


sig_drugs %>% 
  summarise(n = n()) %>% 
  spread(FC, n) %>% 
  arrange(desc(DOWN))

lar_sc_toptable <- readRDS(here::here("results", "lar_sc_exacer_toptable.rds"))

dis_sig_up <- lar_sc_toptable %>% 
  mutate(FC = ifelse(logFC > 0, "UP", "DOWN")) %>% 
  dplyr::group_by(cluster_id, FC) %>% 
  filter(p_adj.loc < 0.1) %>% 
  filter(FC == "UP", contrast == "Ag - Bln")
dis_sig_down <- lar_sc_toptable %>% 
  mutate(FC = ifelse(logFC > 0, "UP", "DOWN")) %>% 
  dplyr::group_by(cluster_id, FC) %>% 
  filter(p_adj.loc < 0.1) %>% 
  filter(FC == "DOWN", contrast == "Ag - Bln")
# 
# jaccard <- function(a, b) {
#   intersection = length(intersect(a, b))
#   union = length(a) + length(b) - intersection
#   return (intersection/union)
# }
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  return (length(intersection))
}

dis_up_drug_down <- result_drug2 %>% 
  mutate(FC = ifelse(logFC > 0, "UP", "DOWN")) %>% 
  dplyr::group_by(sm_name, cell_type_orig, FC) %>% 
  filter(adj.P.Val < 0.2) %>% 
  filter(FC == "DOWN") %>% 
  group_by(sm_name) %>% 
  nest() %>% 
  mutate(result = purrr::map(data, ~{
    drug_sigs <- split(.$gene, .$cell_type_orig)
    disease_sigs <- split(dis_sig_up$gene, dis_sig_up$cluster_id)
    data.frame(bcells = jaccard(drug_sigs$`B cells`, disease_sigs$`B cells`),
               nkcells = jaccard(drug_sigs$`NK cells`, disease_sigs$`NK cells`),
               mnpcells = jaccard(drug_sigs$`Myeloid cells`, disease_sigs$MNP),
               cd4cells = jaccard(drug_sigs$`T cells CD4+`, disease_sigs$`CD4 T cells`),
               cd8cells = jaccard(drug_sigs$`T cells CD8+`, disease_sigs$`CD8 T cells`))
  })) %>% 
  dplyr::select(-data) %>% 
  unnest(result)
dis_up_drug_down2 <- dis_up_drug_down[, -1]
rownames(dis_up_drug_down2) <- dis_up_drug_down$sm_name
Heatmap(dis_up_drug_down2)

dis_down_drug_up <- result_drug2 %>% 
  mutate(FC = ifelse(logFC > 0, "UP", "DOWN")) %>% 
  dplyr::group_by(sm_name, cell_type_orig, FC) %>% 
  filter(adj.P.Val < 0.1) %>% 
  filter(FC == "UP") %>% 
  group_by(sm_name) %>% 
  nest() %>% 
  mutate(result = purrr::map(data, ~{
    drug_sigs <- split(.$gene, .$cell_type_orig)
    disease_sigs <- split(dis_sig_down$gene, dis_sig_down$cluster_id)
    data.frame(bcells = jaccard(drug_sigs$`B cells`, disease_sigs$`B cells`),
               nkcells = jaccard(drug_sigs$`NK cells`, disease_sigs$`NK cells`),
               mnpcells = jaccard(drug_sigs$`Myeloid cells`, disease_sigs$MNP),
               cd4cells = jaccard(drug_sigs$`T cells CD4+`, disease_sigs$`CD4 T cells`),
               cd8cells = jaccard(drug_sigs$`T cells CD8+`, disease_sigs$`CD8 T cells`))
  })) %>% 
  dplyr::select(-data) %>% 
  unnest(result)
dis_down_drug_up2 <- dis_down_drug_up[, -1]
rownames(dis_down_drug_up2) <- dis_down_drug_up$sm_name
Heatmap(dis_down_drug_up2)



a=result_drug2 %>% 
  filter(sm_name == "Belinostat") %>% 
  dplyr::select(cell_type_orig, logFC, gene) %>% 
  spread(gene, logFC)
a2 <- as.matrix(a[, -1])
rownames(a2) <- a$cell_type_orig
Heatmap(a2, name="logFC")

lar_sc_exacer <- readRDS(here::here("results", "lar_sc_exacer.rds"))

drug_perturb <- result_drug2 %>% 
  dplyr::group_by(sm_name) %>% 
  nest() %>% 
  mutate(result = purrr::map(data, ~{
    b = (.) %>% 
      mutate(logFC = logFC) %>% 
      dplyr::select(cell_type_orig, gene, logFC) %>% 
      dplyr::rename(cluster_id = cell_type_orig)
    # b2 <- as.matrix(b[,-1])
    cellnames <- as.character(b$cluster_id)
    cellnames[cellnames == "Myeloid cells"] <- "MNP"
    cellnames[cellnames == "T cells CD4+"] <- "CD4 T cells"
    cellnames[cellnames == "T cells CD8+"] <- "CD8 T cells"
    b$cluster_id <- cellnames
    
   cordat <- b %>% 
      inner_join(lar_sc_exacer, by=c("gene", "cluster_id")) %>% 
      group_by(cluster_id) %>% 
      summarise(cor = cor(logFC.x, logFC.y))
   cordat
  })) %>% 
  dplyr::select(-data) %>% 
  unnest(result)

drug_perturb2 <- result_drug2 %>% 
  dplyr::group_by(sm_name) %>% 
  nest() %>% 
  mutate(result = purrr::map(data, ~{
    b = (.) %>% 
      mutate(logFC=logFC) %>% 
      dplyr::select(cell_type_orig, gene, logFC) %>% 
      dplyr::rename(cluster_id = cell_type_orig)
    # b2 <- as.matrix(b[,-1])
    cellnames <- as.character(b$cluster_id)
    cellnames[cellnames == "Myeloid cells"] <- "MNP"
    cellnames[cellnames == "T cells CD4+"] <- "CD4 T cells"
    cellnames[cellnames == "T cells CD8+"] <- "CD8 T cells"
    b$cluster_id <- cellnames
    
    cordat <- b %>% 
      inner_join(lar_sc_exacer, by=c("gene", "cluster_id")) %>% 
      dplyr::rename(logFC_drug=logFC.x) %>% 
      dplyr::rename(logFC_disease=logFC.y)
    cordat
  })) %>% 
  dplyr::select(-data) %>% 
  unnest(result)

a=drug_perturb %>% 
  spread(sm_name, cor)
a2 <- as.matrix(a[,-1])
rownames(a2) <- a$cluster_id
Heatmap(a2)

smname <- na.omit(drug_perturb) %>% 
  filter(cor > 0) %>% 
  group_by(sm_name) %>% 
  summarise(n = n()) %>% 
  filter(n > 3) %>% 
  pull(sm_name)

smname <- na.omit(drug_perturb) %>% 
  filter(cor < 0) %>% 
  group_by(sm_name) %>% 
  summarise(n = n()) %>% 
  filter(n > 4) %>% 
  pull(sm_name)

ranked_drugs <- drug_perturb %>% 
  filter(sm_name %in% smname) %>% 
  group_by(sm_name) %>% 
  summarise(cor = mean(cor, na.rm = TRUE)) %>% 
  arrange(cor)

drug_perturb2 %>% 
  filter(sm_name %in% c(head(ranked_drugs$sm_name))) %>% 
  ggplot(aes(x = logFC_drug, y = logFC_disease, col = sm_name)) +
  geom_point() +
  ggh4x::facet_grid2(vars(cluster_id), vars(sm_name),
                     scales = "free", independent = "all", space = "free") + 
  stat_smooth(method="lm") +
  theme(legend.position = "none") +
  geom_hline(yintercept =0) +
  geom_vline(xintercept =0) +
  ggrepel::geom_text_repel(aes(label = gene))

lar_sc_exacer_logfc <- readRDS(here::here("results", "lar_sc_exacer_logfc.rds"))
b = result_drug2 %>% 
  filter(sm_name == "BI-D1870") %>%  
  dplyr::select(cell_type_orig, gene, logFC) %>% 
  tidyr::spread(gene, logFC)
b2 <- as.matrix(b[,-1])
cellnames <- as.character(b$cell_type_orig)
cellnames[cellnames == "Myeloid cells"] <- "MNP"
cellnames[cellnames == "T cells CD4+"] <- "CD4 T cells"
cellnames[cellnames == "T cells CD8+"] <- "CD8 T cells"
rownames(b2) <- cellnames
b2[is.na(b2)] <- 0
com_cells <- intersect(rownames(lar_sc_exacer_logfc), rownames(b2))
com_genes <- intersect(colnames(lar_sc_exacer_logfc), colnames(b2))
disease_sig <- lar_sc_exacer_logfc[com_cells, com_genes]
drug_sig <- b2[com_cells, com_genes]
both <- rbind(disease_sig, drug_sig)
cc <- rownames(both)
signature=rep(c("disease", "drug"), each = 5)

library(circlize)
col_fun = colorRamp2(c(min(both), 0, max(both)), c("green", "black", "red"))
row_ha <- rowAnnotation(signature=signature[order(cc)], col = list(
  signature = c("disease"="#D55E00", "drug"="#56B4E9")
))
Heatmap(both[order(cc),], right_annotation = row_ha, col = col_fun,
        border = TRUE, name="logFC", cluster_rows = FALSE)

Heatmap(drug_sig, border = TRUE, name="logFC", cluster_rows = FALSE)
Heatmap(disease_sig, border = TRUE, name="logFC", cluster_rows = FALSE)

a <- lar_sc_toptable %>% 
  filter(p_adj.loc < 0.1)
b = split(a$gene, a$cluster_id)

x = c(1, 2, 3)
y = c(1.1, 10, 3.2)
mean(y) - mean(x)
mean(c(10-1, 1.1-3, 3.2-2))

cc <- "B cells"
dat <- t(data.frame("Asthma" = disease_sig[cc, intersect(colnames(disease_sig), b[[cc]])], 
                    "R428" = drug_sig[cc, intersect(colnames(drug_sig), b[[cc]])]))
col_fun = colorRamp2(c(min(dat), 0, max(dat)), c("green", "white", "red"))
Heatmap(dat, 
        border = TRUE, name="logFC", cluster_rows = FALSE, col = col_fun)

Heatmap(t(dat[1,,drop=FALSE][,order(dat[1,])]), border = TRUE, name="logFC", 
        cluster_columns = FALSE)
Heatmap(t(dat[2,,drop=FALSE][,order(dat[1,])]), border = TRUE, name="logFC", 
        cluster_columns = FALSE)


# aggregation of single-cell to pseudo-bulk data 
pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cell_type_orig", "donor_id", "sm_name"))

# construct design & contrast matrix
ei <- metadata(sce)$experiment_info
ei$group_id <- relevel(ei$group_id, ref = "Bln")
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
contrast <- makeContrasts(Ag-Bln,Dil-Bln,Ag-Dil, levels = mm)

# run DS analysis
res <- pbDS(pb, design = mm, contrast = contrast, method="limma-voom", 
            min_cells = 3)

# access results table
tbl <- lapply(res$table, function(i){
  do.call(rbind, i)
}) %>% 
  do.call(rbind, .)

tbl2 <- tbl %>% 
  filter(gene %in% latephase_exacer_biomarkers) %>% 
  group_by(contrast, cluster_id) %>% 
  mutate(p_adj.loc = p.adjust(p_val, "BH")) %>% 
  arrange(p_val) %>% 
  mutate(n = 1:n())
tbl2






data <- read_h5ad("/Users/asingh/Downloads/GSE193816_all_cells_data.h5ad")
data <- CreateSeuratObject(counts = t(as.matrix(data$X)), 
                           meta.data = data$obs,
                           min.features = 500, 
                           min.cells = 30)


table(data@meta.data$group, data@meta.data$condition)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", 
                           "percent_mito"), ncol = 5)

options(future.globals.maxSize = 8000 * 1024^2)
data <- SCTransform(data) %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

plot_id <- DimPlot(data, reduction = "umap.rna", label = TRUE, 
                   label.size = 2.5, repel = TRUE, group.by = "id") + 
  ggtitle("RNA - id")
plot_cc <- DimPlot(data, reduction = "umap.rna", label = TRUE, 
                   label.size = 2.5, repel = TRUE, group.by = "cluster") + 
  ggtitle("RNA - cell-types")


# extract SCT counts 
sct_counts <- GetAssayData(data, assay = "SCT", layer = "counts")

# extract metadata 
metadata <- data@meta.data %>% 
  dplyr::filter(group == "AA") %>% 
  mutate(sample_id = paste(id, condition, sep="_"))

# create singleCellExperiment object 
sce <- SingleCellExperiment(
  assays = list(counts = sct_counts[, rownames(metadata)]),
  colData = metadata
)

# remove undetected genes
# sce <- sce[rowSums(counts(sce) > 0) > 0, ]
# dim(sce)

# calculate per-cell quality control (QC) metrics
# qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
# ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
# sce <- sce[, !ol]
# dim(sce)

# remove lowly expressed genes 
# sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
# dim(sce)

# compute sum-factors & normalize
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)

(sce <- prepSCE(sce, 
                kid = "cluster", # subpopulation assignments
                gid = "condition",  # group IDs (ctrl/stim)
                sid = "sample_id", # sample IDs (ctrl/stim.1234)
                drop = FALSE))  # drop all other colData columns

# store cluster + sample IDs as well as number of clusters and samples 
nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample
df <- t(table(sce$cluster_id, sce$sample_id))
df

# aggregation of single-cell to pseudo-bulk data 
pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))

library(limma)
# construct design & contrast matrix
ei <- metadata(sce)$experiment_info
ei$group_id <- relevel(ei$group_id, ref = "Bln")
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
contrast <- makeContrasts(Ag-Bln, levels = mm)

# run DS analysis
res <- pbDS(pb, design = mm, contrast = contrast, method="limma-voom", 
            min_cells = 3)

# acess results table
tbl <- lapply(names(res$table$`Ag - Bln`), function(i){
  res$table$`Ag - Bln`[[i]] %>% 
    mutate(cell = i) %>% 
    filter(gene %in% latephase_exacer_biomarkers) %>% 
    mutate(p_adj.loc = p.adjust(p_val, "BH")) %>% 
    arrange(p_val) %>% 
    mutate(n = 1:n())
}) %>% 
  do.call(rbind, .)

tbl %>% 
  ggplot(aes(x = n, y = p_adj.loc, color = cell)) +
  geom_line() +
  scale_x_log10() +
  theme_classic() +
  facet_wrap(~cell) + 
  ggtitle("Adjusted P-values for cell types (RNA)")


tbl_sig <- tbl %>% 
  group_by(cell) %>% 
  dplyr::filter(p_adj.loc < 0.2)

split(tbl_sig$gene, tbl_sig$cell)

tbl_sig %>% 
  ggplot(aes(x = logFC, y = -log10(p_val), color = cell)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = gene)) +
  facet_wrap(~cell)


## mnp
mnp <- read_h5ad("/Users/asingh/Downloads/GSE193816_mnp_data.h5ad")
mnp <- CreateSeuratObject(counts = t(as.matrix(mnp$X)), 
                           meta.data = mnp$obs,
                           min.features = 500, 
                           min.cells = 30)


table(mnp@meta.data$group, mnp@meta.data$condition)
VlnPlot(mnp, features = c("nFeature_RNA", "nCount_RNA", 
                           "percent_mito"), ncol = 5)

options(future.globals.maxSize = 8000 * 1024^2)
mnp <- SCTransform(mnp) %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

plot_id <- DimPlot(mnp, reduction = "umap.rna", label = TRUE, 
                   label.size = 2.5, repel = TRUE, group.by = "id") + 
  ggtitle("RNA - id")
plot_cc <- DimPlot(mnp, reduction = "umap.rna", label = TRUE, 
                   label.size = 2.5, repel = TRUE, group.by = "annotation") + 
  ggtitle("RNA - cell-types")


# extract SCT counts 
sct_counts <- GetAssayData(mnp, assay = "SCT", layer = "counts")

# extract metadata 
metadata <- mnp@meta.data %>% 
  dplyr::filter(group == "AA") %>% 
  mutate(sample_id = paste(id, condition, sep="_"))

# create singleCellExperiment object 
sce <- SingleCellExperiment(
  assays = list(counts = sct_counts[, rownames(metadata)]),
  colData = metadata
)

# compute sum-factors & normalize
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)

(sce <- prepSCE(sce, 
                kid = "annotation", # subpopulation assignments
                gid = "condition",  # group IDs (ctrl/stim)
                sid = "sample_id", # sample IDs (ctrl/stim.1234)
                drop = FALSE))  # drop all other colData columns

# store cluster + sample IDs as well as number of clusters and samples 
nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample
df <- t(table(sce$cluster_id, sce$sample_id))
df

# aggregation of single-cell to pseudo-bulk data 
pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))

library(limma)
# construct design & contrast matrix
ei <- metadata(sce)$experiment_info
ei$group_id <- relevel(ei$group_id, ref = "Bln")
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
contrast <- makeContrasts(Ag-Bln, levels = mm)

# run DS analysis
res <- pbDS(pb, design = mm, contrast = contrast, method="limma-voom", 
            min_cells = 3)

# acess results table
tbl <- lapply(names(res$table$`Ag - Bln`), function(i){
  res$table$`Ag - Bln`[[i]] %>% 
    mutate(cell = i) %>% 
    filter(gene %in% latephase_exacer_biomarkers) %>% 
    mutate(p_adj.loc = p.adjust(p_val, "BH")) %>% 
    arrange(p_val) %>% 
    mutate(n = 1:n())
}) %>% 
  do.call(rbind, .)

tbl %>% 
  ggplot(aes(x = n, y = p_adj.loc, color = cell)) +
  geom_line() +
  scale_x_log10() +
  theme_classic() +
  facet_wrap(~cell) + 
  ggtitle("Adjusted P-values for cell types (RNA)")


tbl_sig <- tbl %>% 
  group_by(cell) %>% 
  dplyr::filter(p_adj.loc < 0.05)

split(tbl_sig$gene, tbl_sig$cell)

tbl_sig %>% 
  ggplot(aes(x = logFC, y = -log10(p_val), color = cell)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = gene)) +
  facet_wrap(~cell)

## t cells
tcell <- read_h5ad("/Users/asingh/Downloads/GSE193816_t_cell_data.h5ad")
tcell <- CreateSeuratObject(counts = t(as.matrix(tcell$X)), 
                          meta.data = tcell$obs,
                          min.features = 500, 
                          min.cells = 30)

library(ggpubr)
p <- tcell@meta.data %>% 
  dplyr::select(id, group, condition, annotation) %>%
  dplyr::group_by(id, group, condition, annotation) %>% 
  filter(group == "AA") %>% 
  mutate(n = n()) %>% 
  ggboxplot(x = "condition", y = "n", color = "condition", 
            palette = "npg",
            add = "jitter",
            facet.by = "annotation", 
            short.panel.labs = FALSE,
            scales = "free")


library(ggplot2)
library(ggpubr)
library(tidyverse)
library(ggsignif)

df = tcell@meta.data %>% 
  dplyr::select(id, group, condition, annotation) %>%
  dplyr::group_by(id, group, condition, annotation) %>% 
  filter(group == "AA") %>% 
  mutate(n = n())

# annotation table with adjusted pvals and y-position of the labels
anno_df = compare_means(n~condition, group.by = "annotation", data = df) %>%
  mutate(y_pos = 40)

df %>% 
  mutate(condition = factor(as.character(condition), 
                            levels = c("Bln", "Dil", "Ag"))) %>% 
  ggplot(aes(x = condition, y = n)) +
  geom_boxplot() +
  facet_wrap(~annotation, scales = "free") +
  stat_compare_means(comparisons = list(c("Bln", "Ag")))







table(tcell@meta.data$group, tcell@meta.data$condition)
VlnPlot(tcell, features = c("nFeature_RNA", "nCount_RNA", 
                          "percent_mito"), ncol = 5)

options(future.globals.maxSize = 8000 * 1024^2)
tcell <- SCTransform(tcell) %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

plot_id <- DimPlot(tcell, reduction = "umap.rna", label = TRUE, 
                   label.size = 2.5, repel = TRUE, group.by = "id") + 
  ggtitle("RNA - id")
plot_cc <- DimPlot(tcell, reduction = "umap.rna", label = TRUE, 
                   label.size = 2.5, repel = TRUE, group.by = "annotation") + 
  ggtitle("RNA - cell-types")


# extract SCT counts 
sct_counts <- GetAssayData(tcell, assay = "SCT", layer = "counts")

# extract metadata 
metadata <- tcell@meta.data %>% 
  dplyr::filter(group == "AA") %>% 
  mutate(sample_id = paste(id, condition, sep="_"))

# create singleCellExperiment object 
sce <- SingleCellExperiment(
  assays = list(counts = sct_counts[, rownames(metadata)]),
  colData = metadata
)

# compute sum-factors & normalize
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)

(sce <- prepSCE(sce, 
                kid = "annotation", # subpopulation assignments
                gid = "condition",  # group IDs (ctrl/stim)
                sid = "sample_id", # sample IDs (ctrl/stim.1234)
                drop = FALSE))  # drop all other colData columns

# store cluster + sample IDs as well as number of clusters and samples 
nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample
df <- t(table(sce$cluster_id, sce$sample_id))
df

# aggregation of single-cell to pseudo-bulk data 
pb <- aggregateData(sce,
                    assay = "counts", fun = "mean",
                    by = c("cluster_id", "sample_id"))

library(limma)
# construct design & contrast matrix
ei <- metadata(sce)$experiment_info
ei$group_id <- relevel(ei$group_id, ref = "Bln")
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
contrast <- makeContrasts(Ag-Bln,Dil-Bln,Ag-Dil, levels = mm)

# run DS analysis
res <- pbDS(pb, design = mm, contrast = contrast, method="limma-voom", 
            min_cells = 3)

# acess results table
tbl <- lapply(res$table, function(i){
  do.call(rbind, i)
}) %>% 
  do.call(rbind, .)

tbl2 <- tbl %>% 
  filter(gene %in% latephase_exacer_biomarkers) %>% 
  group_by(contrast, cluster_id) %>% 
  mutate(p_adj.loc = p.adjust(p_val, "BH")) %>% 
  arrange(p_val) %>% 
  mutate(n = 1:n())
  


tbl <- lapply(names(res$table$`Ag - Bln`), function(i){
  res$table$`Ag - Bln`[[i]] %>% 
    mutate(cell = i) %>% 
    filter(gene %in% latephase_exacer_biomarkers) %>% 
    mutate(p_adj.loc = p.adjust(p_val, "BH")) %>% 
    arrange(p_val) %>% 
    mutate(n = 1:n())
}) %>% 
  do.call(rbind, .)

tbl2 %>% 
  ggplot(aes(x = n, y = p_adj.loc, color = cluster_id)) +
  geom_line() +
  scale_x_log10() +
  theme_classic() +
  facet_wrap(~contrast) + 
  ggtitle("Adjusted P-values for cell types (RNA)")


tbl_sig <- tbl %>% 
  group_by(cell) %>% 
  dplyr::filter(p_adj.loc < 0.05)

split(tbl_sig$gene, tbl_sig$cell)

tbl2 %>% 
  dplyr::filter(p_adj.loc < 0.05) %>% 
  ggplot(aes(x = logFC, y = -log10(p_val), color = cluster_id)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = gene)) +
  facet_grid(contrast~cluster_id)
