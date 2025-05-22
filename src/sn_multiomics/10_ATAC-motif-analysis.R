################################################################################
################################################################################
rm(list = ls()); gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
              "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
library(Seurat);library(Signac);library(EnsDb.Hsapiens.v86);library(BSgenome.Hsapiens.UCSC.hg38);library(AnnotationHub)
library(future);options(future.globals.maxSize = 512 * 1024^3) # setting the maximum memory for 400 GB
plan("multisession", workers = 1)
################################################################################
################################################################################
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/neurosurgery")
setwd(project.dir)
################################################################################
################################################################################
int.metadata.filtered2 <- read_rds("data/derivatives/cells-metadata-filtered.rds") 
samples.meta <- read_csv("data/samples-metadata_ME.csv") 
################################################################################
################################################################################
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
################################################################################
################################################################################
################################################################################
# read the pseudobulk peaks, make a seurat object, normalize
atac.pb.peaks <- read_rds("data/derivatives/manual-boostrapped-pseudobulk/no-weighting-50-or-100-100it-1/atac-peaks/all.rds")

atac.pb.so <- CreateSeuratObject(counts = atac.pb.peaks[,-1], assay = "pb_MACS3_peaks")
# create ATAC assay and add it to the object
atac.pb.so[["ATAC_USE"]] <- CreateChromatinAssay(counts = round(atac.pb.peaks %>% column_to_rownames("peak")),
                                                 sep = c("-", "-"),
                                                 annotation = annotation)
DefaultAssay(atac.pb.so) <- "ATAC_USE"
atac.pb.so
so.meta <- data.frame(cell = colnames(atac.pb.so)) %>%
  mutate(RNA_predicted_cluster = sub("__.*", "", cell),
         sampleID = sub(".*__", "", cell)) %>%
  left_join(int.metadata.filtered2 %>% distinct(sampleID, stimulation, sampleID_2)) %>%
  mutate(cell_stim = paste0(RNA_predicted_cluster, "__", stimulation)) %>%
  column_to_rownames("cell")
atac.pb.so <- AddMetaData(atac.pb.so, metadata = so.meta)
Idents(atac.pb.so) <- "cell_stim"

## normalize
DefaultAssay(atac.pb.so) <- "ATAC_USE"
atac.pb.so <- FindTopFeatures(atac.pb.so, min.cutoff = "q5")
atac.pb.so <- RunTFIDF(atac.pb.so)
atac.pb.so <- RunSVD(atac.pb.so)

################################################################################
################################################################################
# add motif 
library(JASPAR2020);library(TFBSTools);library(patchwork)

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(JASPAR2020, opts = list(species = 9606, all_versions = FALSE))

# add motif information
atac.pb.so <- AddMotifs(atac.pb.so, 
                        genome = BSgenome.Hsapiens.UCSC.hg38,
                        pfm = pfm)
atac.pb.so
write_rds(atac.pb.so, "data/derivatives/manual-boostrapped-pseudobulk/no-weighting-50-or-100-100it-1/atac-peaks/motif-so.rds", compress = "gz")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## differential accessibility on peaks using Seurat LR
registerDoMC(cores = 3)
global.motifs <- foreach(cc = 1:length(unique(atac.pb.so$RNA_predicted_cluster)), .combine = rbind) %dopar% {
  cl <- unique(atac.pb.so$RNA_predicted_cluster)[cc]
  da.peaks <- FindMarkers(atac.pb.so, 
                          ident.1 = paste0(cl, "__stim"), ident.2 = paste0(cl, "__nostim"),
                          only.pos = T, test.use = "LR",
                          latent.vars = "nCount_pb_MACS3_peaks")
  write_rds(da.peaks, paste0("data/derivatives/atac-DA/DA-peaks_", cl, ".rds"), compress = "gz")
  
  ## get top differentially accessible peaks
  top.da.peak <- rownames(da.peaks[da.peaks$p_val < 0.005 & da.peaks$pct.1 > 0.2, ])
  
  ## get overrepresented motifs in the differentially-accessible peaks
  enriched.motifs <- FindMotifs(atac.pb.so, features = top.da.peak) %>%
    mutate(RNA_predicted_cluster = cl) 
  rownames(enriched.motifs) <- NULL
  return(enriched.motifs)
}
system("mkdir -p data/derivatives/atac-DA")
write_rds(global.motifs, "data/derivatives/atac-DA/global-motifs.rds", compress = "gz")
global.motifs <- read_rds("data/derivatives/atac-DA/global-motifs.rds")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
motif.all <- global.motifs %>%
  mutate(source = "global")

library(ggseqlogo)
motif.all %>%
  group_by(source) %>% dplyr::mutate(BH = p.adjust(pvalue, method = "BH")) %>% ungroup() %>%
  dplyr::filter(BH < 0.001) %>%
  ggplot(aes(x = RNA_predicted_cluster, y = `motif.name`, fill = `fold.enrichment`)) +
  geom_tile() +
  scale_fill_gradient2(low = antique.colors[10], high = antique.colors[11], name = "fold enrichment") + 
  my.guides +
  ggh4x::facet_grid2(cols = vars(source), scales = "free", space = "free") +
  labs(x = "cluster", y = "motif name",
       caption = paste0("all motifs shown have an adjusted-pvalue (BH) < 0.001")) +
  bw.theme +
    theme(legend.title = element_text(vjust = 1))
ggsave2("figs/keep/motif-enrichment-global-vs-prm.png", height = 16, width = 8)

MotifPlot(atac.pb.so, 
          motifs = unique(motif.all$motif[which(motif.all$motif.name %in% c("HINFP", "NRF1", "ZBTB33", 
                                                                            "ZBTB14", "TCFL5", "KLF15"))])) +
  scale_fill_manual(values = antique.colors) +
  facet_wrap(~seq_group, scales = "free_x") +
  bw.theme
ggsave2("figs/keep/motif-enrichment-top-6.png", height = 5.5, width = 8)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
