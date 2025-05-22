################################################################################
################################################################################
rm(list = ls()); gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
              "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
library(Seurat);library(Signac)
################################################################################
################################################################################
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/neurosurgery")
setwd(project.dir)
################################################################################
################################################################################
int.metadata.filtered2 <- read_rds("data/derivatives/cells-metadata-filtered.rds")
bs.cell.count <- read_csv("data/derivatives/manual-boostrapped-pseudobulk/counts.csv")
################################################################################
################################################################################
################################################################################
################################################################################
######## using MACS3 peaks
integrated.samples.so <- read_rds("data/derivatives/ME-processing/processed-seurat-objects/integrated-multiome-samples.rds")
atac.peaks <- GetAssayData(integrated.samples.so, assay = "MACS3_peaks", slot = "counts")
atac.peaks.long <- data.table(cell = atac.peaks@Dimnames[[2]],
                              peak = atac.peaks@Dimnames[[1]],
                              count = atac.peaks@x) %>%
  as.data.frame()

registerDoMC(cores = 6)
for(ss in 1:length(unique(int.metadata.filtered2$sampleID_3))) {
  # ss=1
  sample <- unique(int.metadata.filtered2$sampleID_3)[ss]
  
  sample.peaks.long <- atac.peaks.long %>%
    inner_join(int.metadata.filtered2 %>% 
                 dplyr::filter(sampleID_3 == sample) %>%
                 dplyr::select(cell, sampleID_3, RNA_predicted_cluster))
  gc()
  
  s.cl <- unique(bs.cell.count$RNA_predicted_cluster)
  gc()
  s.cl.res <- foreach(cc = 1:length(s.cl), .combine = inner_join) %dopar% {
    cl <- s.cl[cc]
    # get the number of samples to sample
    cl.c <- bs.cell.count$bs_count[which(bs.cell.count$RNA_predicted_cluster == cl)]
    s.cl.bs <- foreach(b = 1:100, .combine = inner_join) %dopar% {
      # extract random # cells
      s.atac.cl <- sample.peaks.long %>%
        dplyr::filter(RNA_predicted_cluster == cl) %>%
        group_by(peak) %>%
        sample_n(size = cl.c, replace = T) %>%
        dplyr::summarise(aggregated_fragment = sum(count))
      colnames(s.atac.cl)[2] <- paste0(cl, "__", sample, "__", b)
      return(s.atac.cl)
    }
    s.cl.bs.2 <- data.frame(peak = s.cl.bs$peak,
                            avg = apply(s.cl.bs[,-1], MARGIN = 1, FUN = function(x) mean(x)))
    colnames(s.cl.bs.2)[2] <- paste0(cl, "__", sample)
    return(s.cl.bs.2)
  }
  system("mkdir -p data/derivatives/manual-boostrapped-pseudobulk/no-weighting-50-or-100-100it-1/atac-peaks")
  write_rds(s.cl.res, paste0("data/derivatives/manual-boostrapped-pseudobulk/no-weighting-50-or-100-100it-1/atac-peaks/", 
                             sample, ".rds"),
            compress = "gz")
}
## read, combine, and save
atac.pb.peaks <- lapply(unique(int.metadata.filtered2$sampleID_3), FUN = function(x){
  read_rds(paste0("data/derivatives/manual-boostrapped-pseudobulk/no-weighting-50-or-100-100it-1/atac-peaks/",x,".rds"))
})
atac.pb.peaks.2 <- full_join(atac.pb.peaks[[1]],atac.pb.peaks[[2]]) %>%
  full_join(atac.pb.peaks[[3]]) %>% full_join(atac.pb.peaks[[4]]) %>% 
  full_join(atac.pb.peaks[[5]]) %>% full_join(atac.pb.peaks[[6]]) %>% 
  full_join(atac.pb.peaks[[7]]) %>% full_join(atac.pb.peaks[[8]]) %>% 
  full_join(atac.pb.peaks[[9]]) %>% full_join(atac.pb.peaks[[10]]) %>% 
  full_join(atac.pb.peaks[[11]]) %>% full_join(atac.pb.peaks[[12]]) %>%
  mutate_at(.vars = vars(-c("peak")), .funs = function(x) replace_na(x, replace = 0))
write_rds(atac.pb.peaks.2, "data/derivatives/manual-boostrapped-pseudobulk/no-weighting-50-or-100-100it-1/atac-peaks/all.rds")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################