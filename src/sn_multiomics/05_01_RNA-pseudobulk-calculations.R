################################################################################
################################################################################
rm(list = ls()); gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
              "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
library(Signac);library(Seurat);library(patchwork);library(future)
options(future.globals.maxSize = 500000 * 1024^2) # setting the maximum memory for 50 GB
plan("multicore", workers = 30)
plan()
registerDoMC(30)
################################################################################
################################################################################
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/neurosurgery")
setwd(project.dir)
################################################################################
################################################################################
int.metadata.filtered2 <- read_rds("data/derivatives/cells-metadata-filtered.rds") 
celltypes.meta <- read_rds(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/data/Allen/human-m1-10x/celltypes-metadata.rds"))
samples.meta <- read_csv("data/samples-metadata_ME.csv")
################################################################################
################################################################################
################################################################################
################################################################################
# pseudobulk is done by basically summing the gene count per cluster/celltype
# nice approach, but I think you need to balance the number cells first

## find the 10th percentile of the number of cells per cluster
## Endo and VLMC dropped because of their terrible quality and being undetected in some participants
bs.cell.count <- int.metadata.filtered2 %>%
  group_by(sampleID_3, RNA_predicted_cluster) %>%
  dplyr::summarise(count = n()) %>% 
  group_by(RNA_predicted_cluster) %>% 
  dplyr::summarise(bs_count = round(quantile(count, 0.1))) %>% ungroup() %>%
  filter(!RNA_predicted_cluster %in% c("Endo", "VLMC"))
write_csv(bs.cell.count, "data/derivatives/manual-boostrapped-pseudobulk/counts.csv")

# read the DCA imputed gene expression data
pload("data/derivatives/rna-imputation/dca/dca-output-all-long.RData.pxz")
gc()

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# for each cluster per sample, randomly select the chosen number from the table above
#   do that for n=100 iterations. 
#   In each iteration, those selected samples will then be aggregated
#   Then, average all iterations

# group by cluster
registerDoMC(cores = 6)
for(ss in 1:length(unique(int.metadata.filtered2$sampleID_3))) {
  sample <- unique(int.metadata.filtered2$sampleID_3)[ss]
  s.rna <- imp.rna.dca %>%
    inner_join(int.metadata.filtered2 %>%
                 dplyr::filter(sampleID_3 == sample) %>%
                 dplyr::select(cell, sampleID_3, nFeature_RNA, 
                               RNA_predicted_cluster))
  s.cl <- unique(bs.cell.count$RNA_predicted_cluster)
  gc()
  s.cl.res <- foreach(cc = 1:length(s.cl), .combine = inner_join) %dopar% {
    cl <- s.cl[cc]
    # get the number of samples to sample
    cl.c <- bs.cell.count$bs_count[which(bs.cell.count$RNA_predicted_cluster == cl)]
    s.cl.bs <- foreach(b = 1:100, .combine = inner_join) %dopar% {
      # extract random # cells
      s.rna.cl <- s.rna %>%
        dplyr::filter(RNA_predicted_cluster == cl) %>%
        group_by(gene) %>%
        # sample_n(size = cl, replace = T, weight = nFeature_RNA) %>%
        sample_n(size = cl.c, replace = T) %>%
        dplyr::summarise(aggregated_imputed_DCA = sum(imputed_DCA))
      colnames(s.rna.cl)[2] <- paste0(cl, "__", sample, "__", b)
      return(s.rna.cl)
    }
    s.cl.bs.2 <- s.cl.bs %>%
      pivot_longer(cols = -c("gene")) %>%
      group_by(gene) %>%
      dplyr::summarise(avg = mean(value, na.rm = T))
    colnames(s.cl.bs.2)[2] <- paste0(cl, "__", sample)
    return(s.cl.bs.2)
  }
  write_rds(s.cl.res, paste0("data/derivatives/manual-boostrapped-pseudobulk/RNA/", 
                             sample, ".rds"),
            compress = "gz")
}

# read and combine them in one dataframe
rna.pb.bs <- lapply(paste0(unique(int.metadata.filtered2$sampleID_3), ".rds"), 
                    function(x) {
                      read_rds(paste0("data/derivatives/manual-boostrapped-pseudobulk/RNA/", x))
                    })
rna.pb.bs <- inner_join(rna.pb.bs[[1]], rna.pb.bs[[2]]) %>%
  inner_join(rna.pb.bs[[3]]) %>%
  inner_join(rna.pb.bs[[4]]) %>%
  inner_join(rna.pb.bs[[5]]) %>%
  inner_join(rna.pb.bs[[6]]) %>%
  inner_join(rna.pb.bs[[7]]) %>%
  inner_join(rna.pb.bs[[8]]) %>%
  inner_join(rna.pb.bs[[9]]) %>%
  inner_join(rna.pb.bs[[10]]) %>%
  inner_join(rna.pb.bs[[11]]) %>%
  inner_join(rna.pb.bs[[12]])
write_rds(rna.pb.bs, "data/derivatives/manual-boostrapped-pseudobulk/RNA/all.rds",
          compress = "gz")
################################################################################
################################################################################
################################################################################
################################################################################
