################################################################################
################################################################################
rm(list = ls()); gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
              "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(AnnotationHub)
library(future)
options(future.globals.maxSize = 512 * 1024^3) # setting the maximum memory for 400 GB
plan("multicore", workers = 30)
plan()
registerDoMC(30)
library(ggstatsplot)
library(patchwork)
################################################################################
################################################################################
#### NOTE
# MOST OF THESE SCRIPTS RELY ON THE FRAGMENTS FILE PATHS INTEGRATED IN THE OBJECT
# THAT'S WHY I RECOMMEND RUNNUNG THE SCRIPT USING ONE MACHINE/LOCATION
####
################################################################################
################################################################################
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/neurosurgery")
setwd(project.dir)
################################################################################
################################################################################
## make samples metadata file.
samples.meta <- data.frame(sex = c("M", "M", "F", "F", "M", "M",
                                   "M", "M", "F", "F", "F", "F"),
                           age = c(52,52,63,63,31,31, 
                                   62,62,58,58,
                                   48,48),
                           anes = c(T, T, F, F, F, F,
                                    F,F,F,F,
                                    T,T),
                           sampleID_3 = c("E3A", "E3B", "E1A", "E1B", "E2A", "E2B",
                                          "G1A", "G1B", "G2A", "G2B",
                                          "E4A", "E4B"),
                           sampleID_4 = rep(c("E3", "E1", "E2", "G1", "G2", "E4"), each = 2),
                           time_from_baseline = c(0, 30, 0, 26, 0, 25,
                                                  0,22,0,22,
                                                  0,35),
                           condition = c(rep("epilepsy", 6),
                                         rep("glioma", 4),
                                         rep("epilepsy", 2)))

################################################################################
################################################################################
# define list of samples to process
# samples is for multiome data
samples <- samples.meta$sampleID_3

################################################################################
################################################################################
################################################################################
################################################################################
a.ref.so <- read_rds(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/data/Allen/human-m1-10x/processed-seurat-object.rds"))
gc()

celltypes.meta <- read_rds(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/data/Allen/human-m1-10x/celltypes-metadata.rds"))
gc()
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
#### STEP 1: read and make Seurat Objects ####

#### Annotation ####
# get gene annotations for hg38
# this annotation will be added to every seurat object for ATAC data
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"


##################################################
## This one is for the multiome data (GEX+ATAC) ##
##################################################
foreach(ss = 1:length(samples)) %dopar% {
  sample <- samples[ss]
  ############# PREPARE
  #### get files
  counts <- Read10X_h5(paste0(project.dir, "/",
                              "data/raw/ME-processing/01/",
                              sample, 
                              "/outs/filtered_feature_bc_matrix.h5"))
  fragpath <- paste0(project.dir, "/",
                     "data/raw/ME-processing/01/",
                     sample, 
                     "/outs/atac_fragments.tsv.gz")
  ############# BUILD
  #### make seurat object
  # create a Seurat object containing the RNA data
  sample.so <- CreateSeuratObject(counts = counts$`Gene Expression`,
                                  assay = "RNA")
  
  # create ATAC assay and add it to the object
  sample.so[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,
                                              sep = c(":", "-"),
                                              fragments = fragpath,
                                              annotation = annotation)
  sample.so
  ############# QC
  #### ATAC mostly
  DefaultAssay(sample.so) <- "ATAC"
  # compute nucleosome signal score per cell
  sample.so <- NucleosomeSignal(sample.so)
  # compute TSS enrichment score per cell
  sample.so <- TSSEnrichment(sample.so)
  # add mt percent
  DefaultAssay(sample.so) <- "RNA"
  sample.so[["percent.mt"]] <- PercentageFeatureSet(sample.so, pattern = "^MT-")
  #### SAVE
  # save the seurat object
  write_rds(sample.so, paste0("data/derivatives/ME-processing/raw-seurat-objects/",
                              sample, ".rds"), compress = "gz")
  #### SUMMARY
  # make a plot to summarize the metrics extracted
  p1 <- DensityScatter(sample.so, 
                       x = 'nCount_ATAC', y = 'TSS.enrichment', 
                       log_x = TRUE, 
                       quantiles = TRUE)
  p2 <- VlnPlot(sample.so,
                features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "percent.mt"),
                ncol = 4,
                pt.size = 0)
  wrap_plots(p1,p2, ncol = 1)
  ggsave(paste0("figs/qc/ME-processing/01_ATAC-", sample, ".png"), bg = "white",
         width = 12, height = 14, units = "in", dpi = 360)
  ############# DONE
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
#### STEP 2: call peaks using MACS3 and filter out low quality cells ####

# I do the filtering for each Seurat object independently
# Then, return them all in a big list, and save

all.samples.so <- list()

for(ss in 1:length(samples)) {
  sample <- samples[ss]
  ############# PREPARE
  sample.so <- read_rds(paste0("data/derivatives/ME-processing/raw-seurat-objects/",
                               sample, ".rds"))
  
  ### add extra columns for the metadata
  md1 <- data.frame(c = colnames(sample.so),
                    sampleID_3 = sample,
                    participant = parse_number(sample),
                    stimulation = ifelse(grepl("B", sample), "stim", "nostim")) %>%
    column_to_rownames("c")
  sample.so <- AddMetaData(sample.so, metadata = md1)
  DefaultAssay(sample.so) <- "ATAC"
  ############# FILTER
  ## filter out low quality cells
  sample.clean <- subset(x = sample.so,
                         subset = nCount_ATAC < 100000 &
                           nCount_ATAC > 1000 &
                           nucleosome_signal < 2 &
                           TSS.enrichment > 1)
  sample.clean$dataset = sample
  ############# PEAKS CALLING
  # call peaks using MACS3
  peaks.macs <- CallPeaks(sample.clean, 
                          macs2.path = "/home/msmuhammad/.local/bin/macs3",
                          format = "BEDPE", additional.args = "-g hs --keep-dup all",
                          outdir = "data/raw/ME-processing/02/")
  # keeps peaks on standard chromosomes only and not in genomic blacklist 
  peaks.macs <- keepStandardChromosomes(peaks.macs, pruning.mode = "coarse")
  peaks.macs <- subsetByOverlaps(x = peaks.macs, 
                                 ranges = blacklist_hg38_unified, # this is manually curated
                                 invert = T)
  ############# QUANTIFICATION
  # quantify counts
  counts.macs <- FeatureMatrix(fragments = Fragments(sample.clean),
                               features = peaks.macs,
                               cells = colnames(sample.clean))
  ############# BUILD AGAIN
  # create a new assay for ATAC with MACS3 peaks
  sample.clean[["MACS3_peaks"]] <- CreateChromatinAssay(counts = counts.macs,
                                                        fragments = fragpath,
                                                        annotation = annotation)
  ############# RETURN
  all.samples.so[[sample]] <- sample.clean
  ############# DONE
}
write_rds(all.samples.so, "data/derivatives/ME-processing/raw-seurat-objects/list-of-all-raw-seurat-objects.rds",
          compress = "gz")

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
#### independent processing for each sample

# this processing includes processing the MACS3 peaks for ATAC and RNA 

# calculate gene activity of ATAC
registerDoMC(cores = 4)
foreach(ss = 1:(length(all.samples)-2)) %dopar% {
  sample <- all.samples[ss]
  sample.so <- all.samples.so[[ss]]
  # drop cells with terrible RNA data
  sample.so <- subset(sample.so,
                      subset = nCount_RNA > 1000 &
                        nCount_RNA < 40000 &
                        percent.mt < 20)
  #### Step 1
  # normalize ATAC data
  DefaultAssay(sample.so) <- "MACS3_peaks"
  sample.so <- FindTopFeatures(sample.so, min.cutoff = "q5")
  sample.so <- RunTFIDF(sample.so)
  sample.so <- RunSVD(sample.so)
  #### Step 2
  # get gene activity matrix without filtering based on gene length
  DefaultAssay(sample.so) <- "MACS3_peaks"
  gene.activities.2 <- GeneActivity(sample.so, 
                                    extend.upstream = 2000,
                                    extend.downstream = 1000, 
                                    process_n = 100000,
                                    biotypes = "protein_coding",
                                    max.width = NULL,
                                    verbose = T)
  # Add the gene activity matrix as a new assay
  sample.so[["Gene.Activity"]] <- CreateAssayObject(counts = gene.activities.2)
  # Set default assay to Gene Activity
  DefaultAssay(sample.so) <- "Gene.Activity"
  # normalize gene activities
  sample.so <- NormalizeData(sample.so)
  sample.so <- ScaleData(sample.so, features = rownames(sample.so))
  sample.so <- FindVariableFeatures(sample.so, assay = "Gene.Activity", 
                                    selection.method = "vst", 
                                    nfeatures = 2000)
  #### Step 3
  # normalize RNA data
  DefaultAssay(sample.so) <- "RNA"
  sample.so <- NormalizeData(sample.so)
  
  sample.so <- ScaleData(sample.so, features = rownames(sample.so))
  sample.so <- FindVariableFeatures(sample.so, assay = "RNA", 
                                    selection.method = "vst", 
                                    nfeatures = 2000)
  #### Step 4
  # identify cell types based on RNA
  sample.so <- RunPCA(sample.so, npcs = 30, verbose = FALSE)
  sample.so <- RunUMAP(sample.so, reduction = "pca", dims = 1:30)
  sample.so <- FindNeighbors(sample.so, reduction = "pca", dims = 1:30)
  sample.so <- FindClusters(sample.so, resolution = 0.5)
  
  transfer.anchors.1 <- FindTransferAnchors(reference = a.ref.so, 
                                            query = sample.so, 
                                            features = VariableFeatures(object = a.ref.so),
                                            reference.reduction = "pca")
  celltype.predictions.1 <- TransferData(anchorset = transfer.anchors.1, 
                                         refdata = a.ref.so$subclass_label,
                                         dims = 1:30) %>%
    rownames_to_column("cell") %>%
    pivot_longer(cols = -c(1:2)) %>%
    mutate(RNA_predicted_subclass_percentage = value * 100,
           name = sub("prediction\\.score\\.", "", name)) %>%
    slice_max(by = cell, order_by = value, n = 1) %>%
    dplyr::filter(name != "max") %>% 
    mutate(cc = name ==predicted.id) %>%
    column_to_rownames("cell") %>%
    dplyr::select(RNA_predicted_subclass=predicted.id, 
                  RNA_predicted_subclass_percentage)
  celltype.predictions.12 <- TransferData(anchorset = transfer.anchors.1, 
                                          refdata = a.ref.so$cluster_label_v2,
                                          dims = 1:30) %>%
    rownames_to_column("cell") %>%
    pivot_longer(cols = -c(1:2)) %>%
    mutate(RNA_predicted_cluster_percentage = value * 100,
           name = sub("prediction\\.score\\.", "", name)) %>%
    slice_max(by = cell, order_by = value, n = 1) %>%
    dplyr::filter(name != "max") %>% 
    mutate(cc = name ==predicted.id) %>%
    column_to_rownames("cell") %>%
    dplyr::select(RNA_predicted_cluster=predicted.id, 
                  RNA_predicted_cluster_percentage)
  sample.so <- AddMetaData(sample.so, metadata = celltype.predictions.1)
  sample.so <- AddMetaData(sample.so, metadata = celltype.predictions.12)
  #### Step 5
  # identify cell types based on ATAC
  DefaultAssay(sample.so) <- "Gene.Activity"
  sample.so <- RunPCA(sample.so, npcs = 30, verbose = FALSE)
  sample.so <- RunUMAP(sample.so, reduction = "pca", dims = 1:30)
  sample.so <- FindNeighbors(sample.so, reduction = "pca", dims = 1:30)
  sample.so <- FindClusters(sample.so, resolution = 0.5)
  transfer.anchors.2 <- FindTransferAnchors(reference = a.ref.so, 
                                            query = sample.so, 
                                            features = VariableFeatures(object = a.ref.so),
                                            reference.reduction = "pca")
  celltype.predictions.2 <- TransferData(anchorset = transfer.anchors.2, 
                                         refdata = a.ref.so$subclass_label,
                                         dims = 1:30) %>%
    rownames_to_column("cell") %>%
    pivot_longer(cols = -c(1:2)) %>%
    mutate(ATAC_predicted_subclass_percentage = value * 100,
           name = sub("prediction\\.score\\.", "", name)) %>%
    slice_max(by = cell, order_by = value, n = 1) %>%
    dplyr::filter(name != "max") %>% 
    mutate(cc = name ==predicted.id) %>%
    column_to_rownames("cell") %>%
    dplyr::select(ATAC_predicted_subclass=predicted.id, 
                  ATAC_predicted_subclass_percentage)
  celltype.predictions.22 <- TransferData(anchorset = transfer.anchors.2, 
                                          refdata = a.ref.so$cluster_label_v2,
                                          dims = 1:30) %>%
    rownames_to_column("cell") %>%
    pivot_longer(cols = -c(1:2)) %>%
    mutate(ATAC_predicted_cluster_percentage = value * 100,
           name = sub("prediction\\.score\\.", "", name)) %>%
    slice_max(by = cell, order_by = value, n = 1) %>%
    dplyr::filter(name != "max") %>% 
    mutate(cc = name ==predicted.id) %>%
    column_to_rownames("cell") %>%
    dplyr::select(ATAC_predicted_cluster=predicted.id, 
                  ATAC_predicted_cluster_percentage)
  sample.so <- AddMetaData(sample.so, metadata = celltype.predictions.2)
  sample.so <- AddMetaData(sample.so, metadata = celltype.predictions.22)
  #### Step 6
  # compare RNA to ATAC activity
  ga <- sample.so@assays$Gene.Activity$scale.data %>% as.data.frame()
  rna <- sample.so@assays$RNA$scale.data %>% as.data.frame()
  min_cells <- 0.75 * ncol(sample.so) # at least 25% of cells
  ga.tf <- ga %>%
    mutate_all(.funs = function(x) ifelse(x==0, T, F))
  ga.25 <- apply(ga, 2, quantile, probs = 0.25)
  ga.25.2 <- sweep(ga, 2, ga.25, FUN = ">") 
  ga.25.stats <- data.frame(gene = rownames(ga),
                            zz = rowSums(ga.tf),
                            pct_25 = rowSums(ga.25.2)) %>%
    dplyr::filter(zz < min_cells,
                  pct_25 <min_cells)
  rna.tf <- rna %>%
    mutate_all(.funs = function(x) ifelse(x==0, T, F))
  rna.25 <- apply(rna, 2, quantile, probs = 0.25)
  rna.25.2 <- sweep(rna, 2, ga.25, FUN = ">") 
  rna.25.stats <- data.frame(gene = rownames(rna),
                             zz = rowSums(rna.tf),
                             pct_25 = rowSums(rna.25.2)) %>%
    dplyr::filter(zz < min_cells,
                  pct_25 <min_cells)
  int <- rna.25.stats$gene[rna.25.stats$gene %in% ga.25.stats$gene]
  ga.4 <- ga[int,]
  rna.4 <- rna[int,]
  all(rownames(ga.4)==rownames(rna.4))
  all(colnames(ga.4)==colnames(rna.4))
  # get correlation between the gene expression vector per cell from the RNA and the ATAC
  corr.2 <- foreach(kk = 1:ncol(ga.4), .combine = rbind) %dopar% {
    cc <- cor.test(ga.4[,kk], rna.4[,kk], method = "spearman")
    data.frame(cell = colnames(ga.4)[kk],
               r = cc$estimate,
               pval = cc$p.value)
  }
  corr.2 <- corr.2 %>%
    mutate(count_intersecting_genes_GEX_ATAC = length(int))
  sample.so <- AddMetaData(sample.so, metadata = corr.2[,-1])
  #### Step 7
  # return
  write_rds(sample.so,
            paste0("data/derivatives/ME-processing/processed-seurat-objects/separate/", 
                   sample, ".rds"),
            compress = "gz")
  gc()
}
gc()


all.samples.so.3 <- lapply(samples, 
                           function(x) {
                             read_rds(paste0("data/derivatives/ME-processing/processed-seurat-objects/separate/", x, ".rds"))
                           })
names(all.samples.so.3) <- samples
## pload is a customized function. check my workbench GitHub repo for more info
psave(all.samples.so.3, file = "data/derivatives/ME-processing/processed-seurat-objects/list-of-independently-processed-multiome-seurat-objects-w-labels.rda")
pload("data/derivatives/ME-processing/processed-seurat-objects/list-of-independently-processed-multiome-seurat-objects-w-labels.rda.pxz")
gc()
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# combine results/metadata of all samples here
samples.meta <- foreach(ss=1:length(all.samples.so.3), .combine = full_join) %dopar% {
  cbind(all.samples.so.3[[ss]]@meta.data, 
        cell = colnames(all.samples.so.3[[ss]])) %>%
    relocate(cell)
}
################################################################################
################################################################################
################################################################################
################################################################################
################################# Integration ##################################
################################################################################
################################################################################
################################################################################
################################################################################
## integrate based on RNA
all.samples.so.4 <- all.samples.so.3 # dropped samples that you don't want for integration here if needed

gc()
for (i in 1:length(all.samples.so.4)) {
  DefaultAssay(all.samples.so.4[[i]]) <- "RNA"
}

# select features that are highly variable across all samples (mutliome only)
int.feat <- SelectIntegrationFeatures(object.list = all.samples.so.4, 
                                      assay = rep("RNA", length(all.samples.so.4)))
# find anchors (HUGE and takes a LOT of time)
int.anchors <- FindIntegrationAnchors(object.list = all.samples.so.4,
                                      anchor.features = int.feat, reduction = "rpca",
                                      assay = rep("RNA", length(all.samples.so.4)))

# integrate (HUGE and takes a LOT of time)
integrated.samples.so <- IntegrateData(anchorset = int.anchors, verbose = T, 
                                       dims = 1:20, preserve.order = T)

# process integrated object
DefaultAssay(integrated.samples.so) <- "integrated"
integrated.samples.so <- ScaleData(integrated.samples.so, verbose = FALSE)
integrated.samples.so <- RunPCA(integrated.samples.so, npcs = 30, verbose = FALSE)
integrated.samples.so <- RunUMAP(integrated.samples.so, reduction = "pca", dims = 1:30)
integrated.samples.so <- FindNeighbors(integrated.samples.so, reduction = "pca", dims = 1:30)
integrated.samples.so <- FindClusters(integrated.samples.so, resolution = 0.5)
#### predict celltype based on the integration 
####  this is not great since it's only looking at the top variable features and relies on them
int.transfer.anchors <- FindTransferAnchors(reference = a.ref.so,
                                            query = integrated.samples.so,
                                            features = VariableFeatures(object = a.ref.so),
                                            reference.reduction = "pca")
# subclass prediction
int.predictions.1 <- TransferData(anchorset = int.transfer.anchors,
                                  refdata = a.ref.so$subclass_label,
                                  dims = 1:30) %>%
  rownames_to_column("cell") %>%
  pivot_longer(cols = -c(1:2)) %>%
  mutate(integration_predicted_subclass_percentage = value * 100,
         name = sub("prediction\\.score\\.", "", name)) %>%
  slice_max(by = cell, order_by = value, n = 1) %>%
  dplyr::filter(name != "max") %>% 
  mutate(cc = name ==predicted.id) %>%
  column_to_rownames("cell") %>%
  dplyr::select(integration_predicted_subclass=predicted.id, 
                integration_predicted_subclass_percentage)
# cluster prediction
int.predictions.12 <- TransferData(anchorset = int.transfer.anchors, 
                                   refdata = a.ref.so$cluster_label_v2,
                                   dims = 1:30) %>%
  rownames_to_column("cell") %>%
  pivot_longer(cols = -c(1:2)) %>%
  mutate(integration_predicted_cluster_percentage = value * 100,
         name = sub("prediction\\.score\\.", "", name)) %>%
  slice_max(by = cell, order_by = value, n = 1) %>%
  dplyr::filter(name != "max") %>% 
  mutate(cc = name ==predicted.id) %>%
  column_to_rownames("cell") %>%
  dplyr::select(integration_predicted_cluster=predicted.id, 
                integration_predicted_cluster_percentage)
integrated.samples.so <- AddMetaData(integrated.samples.so, metadata = int.predictions.1)
integrated.samples.so <- AddMetaData(integrated.samples.so, metadata = int.predictions.12)
DefaultAssay(integrated.samples.so) <- "RNA"
integrated.samples.so <- JoinLayers(object = integrated.samples.so)

### save
write_rds(integrated.samples.so, 
          "data/derivatives/ME-processing/processed-seurat-objects/integrated-multiome-samples.rds",
          compress = "gz")

# save metadata
# check labeling differences between RNA independent and integration
int.metadata <- integrated.samples.so@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  inner_join(integrated.samples.so@meta.data %>% 
               as.data.frame() %>%
               rownames_to_column("cell")) %>%
  mutate(cell_2 = sub("_.*", "", cell)) %>%
  dplyr::select(-c(participant, dataset, `orig.ident`)) %>%
  left_join(samples.meta)
write_rds(int.metadata, "data/derivatives/RNA-cells-metadata.rds", compress = "gz")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# filter the cells to keep ones w good quality per cluster

# drop cells with low celltype prediction confidence
int.metadata.filtered <- int.metadata %>%
  group_by(RNA_predicted_cluster) %>%
  dplyr::filter(RNA_predicted_cluster_percentage > 75)
write_rds(int.metadata.filtered, "data/derivatives/cells-metadata-filtered.rds")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################