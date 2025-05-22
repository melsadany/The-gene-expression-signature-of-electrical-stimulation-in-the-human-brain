################################################################################
################################################################################
rm(list = ls()); gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
              "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
library(Signac)
library(Seurat)
library(patchwork)
library(future)
options(future.globals.maxSize = 500000 * 1024^2) # setting the maximum memory for 50 GB
plan("multicore", workers = 30)
plan()
registerDoMC(30)
source("src/98_manuscript-color-palettes.R")
################################################################################
################################################################################
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/neurosurgery")
setwd(project.dir)
################################################################################
################################################################################
integrated.samples.so <- read_rds("data/derivatives/ME-processing/processed-seurat-objects/integrated-multiome-samples.rds")

int.metadata <- read_rds("data/derivatives/RNA-cells-metadata.rds")
celltypes.meta <- read_rds(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/data/Allen/human-m1-10x/celltypes-metadata.rds"))
##
samples.meta <- read_csv("data/samples-metadata_ME.csv")
pc.genes <- read_rds("../../data/genomics/bioMart-protein-coding-genes.rds")

################################################################################
################################################################################
################################################################################
################################################################################
# write raw rna counts for imputation

int.rna <- GetAssayData(integrated.samples.so, 
                        assay = "RNA", layer = "counts") %>%
  as.data.frame()
rownames(int.rna) <- rownames(integrated.samples.so@assays$RNA@features)
colnames(int.rna) <- rownames(integrated.samples.so@assays$RNA@cells)

# keep protein-coding genes only
int.rna <- t(int.rna) %>%
  as.data.frame() %>% 
  dplyr::select(any_of(pc.genes$external_gene_name))

# drop genes that are not expressed at all from all cells in all samples
g.su <- as.data.frame(colSums(int.rna)) %>% 
  rownames_to_column("gene") %>%
  rename("sum" = 2) %>%
  dplyr::filter(sum > 0)

# keeps cells with at least 1500 total count per cell
c.su <- as.data.frame(rowSums(int.rna)) %>% 
  rownames_to_column("cell") %>%
  rename("sum" = 2)

dim(int.rna)
int.rna <- int.rna[c.su$cell,] %>%
  dplyr::select(g.su$gene)
dim(int.rna)


# save stimulated cells and unstimulated cells separate
stim.cells <- int.metadata$cell[which(int.metadata$stimulation == "stim")]
nostim.cells <- int.metadata$cell[which(int.metadata$stimulation == "nostim")]

### for DCA
t1 <- as.data.frame(t(int.rna[rownames(int.rna) %in% stim.cells,])) %>%
  rownames_to_column("gene")
write_csv(t1, "data/derivatives/rna-imputation/dca/dca-input-integrated-stimulated-rna.csv")
write_csv(colnames(t1) %>% as.data.frame(), 
          "data/derivatives/rna-imputation/dca/dca-input-integrated-stimulated-rna_cells.csv")
rm(t1);gc()
t1 <- as.data.frame(t(int.rna[rownames(int.rna) %in% nostim.cells,])) %>%
  rownames_to_column("gene")
write_csv(t1, "data/derivatives/rna-imputation/dca/dca-input-integrated-unstimulated-rna.csv")
write_csv(colnames(t1) %>% as.data.frame(),
          "data/derivatives/rna-imputation/dca/dca-input-integrated-unstimulated-rna_cells.csv")

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# check the imputation by DCA
# read imputed
imp.rna.dca.s <- fread("data/derivatives/rna-imputation/dca/out_stimulated/mean.tsv")
imp.rna.dca.ns <- fread("data/derivatives/rna-imputation/dca/out_unstimulated/mean.tsv")

tmp <- inner_join(imp.rna.dca.s, imp.rna.dca.ns) %>%
  rename("gene" = 1)
psave(tmp, file = "data/derivatives/rna-imputation/dca/dca-output-all.RData")
rm(tmp);gc()

imp.rna.dca <- inner_join(imp.rna.dca.s, imp.rna.dca.ns) %>%
  rename("gene" = 1) %>%
  pivot_longer(cols = -c("gene"), names_to = "cell", values_to = "imputed_DCA")
psave(imp.rna.dca, file = "data/derivatives/rna-imputation/dca/dca-output-all-long.RData")
rm(imp.rna.dca.s);rm(imp.rna.dca.ns);gc()
################################################################################
################################################################################
################################################################################
################################################################################