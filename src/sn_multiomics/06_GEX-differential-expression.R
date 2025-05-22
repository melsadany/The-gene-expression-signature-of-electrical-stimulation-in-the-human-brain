################################################################################
################################################################################
rm(list = ls()); gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
              "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
library(patchwork)
registerDoMC(30)
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
################################################################################
# read the pseudobulk data
rna.pb.bs <- read_rds("data/derivatives/manual-boostrapped-pseudobulk/RNA/all.rds")

## drop genes with 0 count total
rna.pb.bs <- rna.pb.bs[rowSums(rna.pb.bs[,-1])>0,]

rna.pb.bs.meta <- data.frame(col = colnames(rna.pb.bs)[-1]) %>%
  mutate(RNA_predicted_cluster = sub("__.*", "", col),
         sampleID_3 = sub(".*__", "", col))
################################################################################
################################################################################
################################################################################
# I tried several methods to do the differential gene expression analysis
################################################################################
################################################################################
################################################################################
library(edgeR)
registerDoMC(cores = 6)
edger.res <- foreach(cc = 1:length(unique(rna.pb.bs.meta$RNA_predicted_cluster)), .combine = rbind) %dopar% {
  cl <- unique(rna.pb.bs.meta$RNA_predicted_cluster)[cc]
  c.names <- rna.pb.bs.meta %>%
    dplyr::filter(RNA_predicted_cluster == cl)
  c.df <- rna.pb.bs %>%
    column_to_rownames("gene") %>%
    dplyr::select(c.names$col)
  cl.meta <- c.names %>%
    left_join(samples.meta) %>%
    mutate(stimulation = ifelse(grepl("A", sampleID_3), 0, 1),
           sexM = ifelse(sex=="M",1,0),
           anes = as.numeric(anes))
  
  ## epilepsy
  cl.meta.ep <- cl.meta %>% dplyr::filter(condition =="epilepsy") %>% mutate(sampleID_4 = factor(sampleID_4, levels = unique(sampleID_4)))
  c.df.ep <- c.df %>% dplyr::select(cl.meta.ep$col)
  # design
  design.ep <- model.matrix(~ stimulation, data = cl.meta.ep)
  rownames(design.ep) <- colnames(c.df.ep)
  y.ep <- DGEList(counts = c.df.ep, genes = rownames(c.df.ep))
  y.ep <- calcNormFactors(y.ep)
  y.ep <- estimateDisp(y.ep, design.ep, robust=TRUE)
  fit.ep <- glmFit(y.ep, design.ep)
  lrt.ep <- glmLRT(fit.ep)
  de.ep <- decideTestsDGE(lrt.ep)
  glm.res.ep <- as.data.frame(topTags(lrt.ep, n = nrow(c.df.ep))) %>%
    mutate(RNA_predicted_cluster = cl, dx = "epilepsy")
  
  ## glioma
  cl.meta.gl <- cl.meta %>% dplyr::filter(condition =="glioma") %>% mutate(sampleID_4 = factor(sampleID_4, levels = unique(sampleID_4)))
  c.df.gl <- c.df %>% dplyr::select(cl.meta.gl$col)
  # design
  design.gl <- model.matrix(~ stimulation, data = cl.meta.gl)
  rownames(design.gl) <- colnames(c.df.gl)
  y.gl <- DGEList(counts = c.df.gl, genes = rownames(c.df.gl))
  y.gl <- calcNormFactors(y.gl)
  y.gl <- estimateDisp(y.gl, design.gl, robust=TRUE)
  fit.gl <- glmFit(y.gl, design.gl)
  lrt.gl <- glmLRT(fit.gl)
  de.gl <- decideTestsDGE(lrt.gl)
  glm.res.gl <- as.data.frame(topTags(lrt.gl, n = nrow(c.df.gl))) %>%
    mutate(RNA_predicted_cluster = cl, dx = "glioma")
  
  ## combined
  cl.meta.co <- cl.meta %>% mutate(sampleID_4 = factor(sampleID_4, levels = unique(sampleID_4)))
  c.df.co <- c.df %>% dplyr::select(cl.meta.co$col)
  # design
  design.co <- model.matrix(~ condition + anes + sexM + age + stimulation, data = cl.meta.co)
  rownames(design.co) <- colnames(c.df.co)
  y.co <- DGEList(counts = c.df.co, genes = rownames(c.df.co))
  y.co <- calcNormFactors(y.co)
  y.co <- estimateDisp(y.co, design.co, robust=TRUE)
  fit.co <- glmFit(y.co, design.co)
  lrt.co <- glmLRT(fit.co)
  de.co <- decideTestsDGE(lrt.co)
  glm.res.co <- as.data.frame(topTags(lrt.co, n = nrow(c.df.co))) %>%
    mutate(RNA_predicted_cluster = cl, dx = "combined")
  
  rbind(glm.res.ep, glm.res.gl, glm.res.co)
}
# save results
system("mkdir -p data/derivatives/gex-DE/final")
write_rds(edger.res, "data/derivatives/gex-DE/final/edgeR-all-w-dx.rds",compress = "gz")
################################################################################
################################################################################
################################################################################
################################################################################
library(glmmSeq)
registerDoMC(cores = 10)
lmmseq.res <- foreach(cc = 1:length(unique(rna.pb.bs.meta$RNA_predicted_cluster)), .combine = rbind) %dopar% {
  cl <- unique(rna.pb.bs.meta$RNA_predicted_cluster)[cc]
  c.names <- rna.pb.bs.meta %>%
    dplyr::filter(RNA_predicted_cluster == cl)
  gc()
  c.df <- rna.pb.bs %>%
    column_to_rownames("gene") %>%
    dplyr::select(c.names$col)
  cl.meta <- c.names %>%
    left_join(samples.meta) %>%
    mutate(stimulation = ifelse(grepl("A", sampleID_3), 0, 1),
           sexM = ifelse(sex=="M",1,0),
           anes = as.numeric(anes))
  
  ## epilepsy
  cl.meta.ep <- cl.meta %>% filter(condition=="epilepsy") %>% mutate(sampleID_4 = factor(sampleID_4, levels = unique(sampleID_4)))
  c.df.ep <- c.df %>% select(cl.meta.ep$col)
  c.df.log.ep <- log2(c.df.ep[!rowSums(c.df.ep)==0,]+1)
  lmmres.ep <- lmmSeq(~ stimulation + (1 | sampleID_4),
                      maindata = c.df.log.ep, metadata = cl.meta.ep, progress = TRUE, cores = 10)
  lmmres.ep <- glmmQvals(lmmres.ep)
  tt.ep <- apply(lmmres.ep@predict, 1, function(x) {log2(mean(x[c(2,4,6)])/mean(x[c(1,3,5)]))}) %>%
    as.data.frame() %>%dplyr::rename("logFC" = 1) %>%rownames_to_column("gene")
  lmmres.df.ep <- lmmres.ep@stats %>% as.data.frame() %>% 
    select(`coef.stimulation`, pval=`stimulation.2`,qval=`stimulation.3`) %>%
    rownames_to_column("gene") %>%
    mutate(RNA_predicted_cluster = cl,FDR = p.adjust(pval, method = "fdr"),BF = p.adjust(pval, method = "bonferroni")) %>%
    left_join(tt.ep) %>%relocate(gene, coef.stimulation, logFC, pval, qval, FDR, BF) %>%
    mutate(dx = "epilepsy")
  
  
  ## glioma
  cl.meta.gl <- cl.meta %>% filter(condition=="glioma") %>% mutate(sampleID_4 = factor(sampleID_4, levels = unique(sampleID_4)))
  c.df.gl <- c.df %>% select(cl.meta.gl$col)
  c.df.log.gl <- log2(c.df.gl[!rowSums(c.df.gl)==0,]+1)
  lmmres.gl <- lmmSeq(~ stimulation + (1 | sampleID_4),
                      maindata = c.df.log.gl, metadata = cl.meta.gl, progress = TRUE, cores = 10)
  lmmres.gl <- glmmQvals(lmmres.gl)
  tt.gl <- apply(lmmres.gl@predict, 1, function(x) {log2(mean(x[c(2,4,6)])/mean(x[c(1,3,5)]))}) %>%
    as.data.frame() %>%dplyr::rename("logFC" = 1) %>%rownames_to_column("gene")
  lmmres.df.gl <- lmmres.gl@stats %>% as.data.frame() %>% 
    select(`coef.stimulation`, pval=`stimulation.2`,qval=`stimulation.3`) %>%
    rownames_to_column("gene") %>%
    mutate(RNA_predicted_cluster = cl,FDR = p.adjust(pval, method = "fdr"),BF = p.adjust(pval, method = "bonferroni")) %>%
    left_join(tt.gl) %>%relocate(gene, coef.stimulation, logFC, pval, qval, FDR, BF) %>%
    mutate(dx = "glioma")
  
  
  ## combined
  cl.meta.co <- cl.meta %>% mutate(sampleID_4 = factor(sampleID_4, levels = unique(sampleID_4)))
  c.df.co <- c.df %>% select(cl.meta.co$col)
  c.df.log.co <- log2(c.df.co[!rowSums(c.df.co)==0,]+1)
  lmmres.co <- lmmSeq(~ stimulation + (1 | sampleID_4),
                      maindata = c.df.log.co, metadata = cl.meta.co, progress = TRUE, cores = 10)
  lmmres.co <- glmmQvals(lmmres.co)
  tt.co <- apply(lmmres.co@predict, 1, function(x) {log2(mean(x[c(2,4,6)])/mean(x[c(1,3,5)]))}) %>%
    as.data.frame() %>%dplyr::rename("logFC" = 1) %>%rownames_to_column("gene")
  lmmres.df.co <- lmmres.co@stats %>% as.data.frame() %>% 
    select(`coef.stimulation`, pval=`stimulation.2`,qval=`stimulation.3`) %>%
    rownames_to_column("gene") %>%
    mutate(RNA_predicted_cluster = cl,FDR = p.adjust(pval, method = "fdr"),BF = p.adjust(pval, method = "bonferroni")) %>%
    left_join(tt.co) %>%relocate(gene, coef.stimulation, logFC, pval, qval, FDR, BF) %>%
    mutate(dx = "combined")
  
  rbind(lmmres.df.ep, lmmres.df.gl, lmmres.df.co)
}
# save results
write_rds(lmmseq.res, "data/derivatives/gex-DE/final/lmmseq-all-w-dx.rds",compress = "gz")
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
# read results and find the agreement of edgeR and lmmSeq
lmmseq.cl <- read_rds("data/derivatives/gex-DE/final/lmmseq-all-w-dx.rds")
edger.cl <- read_rds("data/derivatives/gex-DE/final/edgeR-all-w-dx.rds")

all.cl <- inner_join(lmmseq.cl %>% 
                       dplyr::select(dx, RNA_predicted_cluster, gene, lmmSeq_coef_stimulation = coef.stimulation, lmmSeq_FDR = FDR,lmmSeq_pval = pval),
                     edger.cl %>% 
                       dplyr::select(dx, RNA_predicted_cluster, gene = genes, edgeR_logFC = logFC, edgeR_FDR = FDR, edgeR_pval = PValue)) %>%
  dplyr::select(dx, RNA_predicted_cluster, gene, lmmSeq_coef_stimulation, edgeR_logFC, lmmSeq_FDR, edgeR_FDR,lmmSeq_pval, edgeR_pval)

################################################################################
################################################################################
################################################################################
################################################################################
## define thresholds for significance here
all.cl.2 <- all.cl %>%
  mutate(lmmSeq_DEG_cons = ifelse(lmmSeq_FDR<= 0.05&abs(lmmSeq_coef_stimulation)>1, T, F), 
         edgeR_DEG_cons = ifelse(edgeR_FDR<=0.05&abs(edgeR_logFC)>1, T, F),
         sig_both_cons = ifelse(lmmSeq_DEG_cons&edgeR_DEG_cons, T, F))
################################################################################
################################################################################
write_csv(all.cl.2, "data/derivatives/gex-DE/final/all-results.csv")
all.cl.2 <- read_csv("data/derivatives/gex-DE/final/all-results.csv")
################################################################################
################################################################################
################################################################################
#### you have different criterial for thresholds
thresholds <- data.frame(name = c("cons"),
                         pval_thr = c(0.05),
                         lfc_thr = c(1),
                         coef_thr = c(1)) 
foreach(c = 1:nrow(thresholds)) %dopar% {
  thr.name <- thresholds$name[c]
  p.thr <- thresholds$pval_thr[c]
  co.thr <- thresholds$coef_thr[c]
  lfc.thr <- thresholds$lfc_thr[c]
  
  df <- all.cl.2 %>% 
    select(dx, RNA_predicted_cluster, gene, edgeR_logFC, lmmSeq_coef_stimulation,ends_with(thr.name)) %>%
    rename(edgeR_DEG = paste0("edgeR_DEG_", thr.name),
           lmmSeq_DEG = paste0("lmmSeq_DEG_", thr.name),
           sig_both_DEG = paste0("sig_both_", thr.name)) %>%
    mutate(dx = factor(dx, levels = c("epilepsy", "glioma", "combined")),
           group = ifelse(sig_both_DEG==T, "DEG in both",ifelse(edgeR_DEG == T, "edgeR DEG ONLY",ifelse(lmmSeq_DEG == T, "lmmSeq DEG ONLY", "ns"))))
  
  df %>% 
    ggplot(aes(x = lmmSeq_coef_stimulation, y = edgeR_logFC,  label = ifelse(group != "ns", gene,""),color = group)) +
    geom_point(aes(alpha = ifelse(sig_both_DEG==F, F, T))) +
    ggrepel::geom_text_repel(max.overlaps = 100, size = 2) +
    geom_hline(yintercept = c(-lfc.thr,lfc.thr), linetype = 2, color = "red", alpha = 0.5) +
    geom_vline(xintercept = c(-co.thr,co.thr), linetype = 2, color = "red", alpha = 0.5) +
    geom_hline(yintercept = 0, linetype=2, color="black",alpha =0.2)+geom_vline(xintercept = 0, linetype=2, color="black",alpha =0.2)+
    scale_color_manual(values = c(abstract.colors.2[1:3], "grey"),name = "") +
    scale_alpha_manual(values = c(0.3,1), guide = "none") +
    ggh4x::facet_grid2(cols = vars(dx), rows = vars(RNA_predicted_cluster), scales = "fixed", space = "free") +
    labs(x = "lmmSeq coefficient for stimulation",
         y = expression("log"[2]~"FC from edgeR"),
         caption = paste0("n(pairs): 6 (epilepsy:4, glioma:2)\n",
                          "edgeR model:\n",
                          "      combined: ~ dx + anesthesia + sexM + age + stimulation\n",
                          "      epilepsy/glioma: ~ stimulation\n",
                          "lmmSeq model:\n",
                          "      combined: ~ stimulation + anesthesia + sexM + age + dx + (1|participant)\n",
                          "      epilepsy/glioma: ~ stimulation + (1|participant)\n",
                          "DEGs criteria:      (abs(logFC) >",lfc.thr," | abs(coefficient) >", co.thr,  ") & FDR < ", p.thr)) +
    bw.theme
  ggsave2(paste0("figs/keep/lmmseq-vs-edgeR-DGE-results_", thr.name, "-thr.pdf"),width = 8, height = 13)
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
