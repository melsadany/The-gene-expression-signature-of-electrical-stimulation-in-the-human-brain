################################################################################
################################################################################
rm(list = ls());gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
              "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
library(granulator)
################################################################################
################################################################################
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/neurosurgery")
setwd(project.dir)
source("src/98_manuscript-color-palettes.R")
################################################################################
################################################################################
genes.meta <- read_rds("../../data/genomics/bioMart-genes.rds") %>%
  group_by(ensembl_gene_id, hgnc_symbol) %>%
  dplyr::summarise(transcript_length = max(transcript_length)) %>% ungroup()
################################################################################
################################################################################
# get the bulk data
mouse.bulk.r <- fread("data/raw/GEO/GSE224952_Mouse_StimulatedVSBaseline.counts.txt.gz") %>%
  select(-id) %>% relocate(symbol)
human.nostim.bulk.r <- fread("data/raw/GEO/GSE224952_combined_UnstimulatedVSBaseline.counts.txt.gz")
human.stim.bulk.r <- fread("data/raw/GEO/GSE224952_combined_StimulatedVSBaseline.counts.txt.gz")

# reference data
library(Seurat)
allen.human <- read_rds("../../data/Allen/human-m1-10x/subclass-aggregated-sig.rds")
allen.mouse <- read_rds("../../data/Allen/mouse/subclass-aggregated-sig.rds")

celltypes.meta.human <- read_rds("../../data/Allen/human-m1-10x/celltypes-metadata.rds")
celltypes.meta.mouse <- read_rds("../../data/Allen/mouse/celltypes-metadata.rds") %>%
  mutate(cluster_label_v2 = sub("-PVM", "", cluster_label_2))

################################################################################
################################################################################
# convert to gene symbol

human.nostim.bulk <- inner_join(genes.meta[,-3],
                                human.nostim.bulk.r %>% rename(ensembl_gene_id = id)) %>%
  filter(nchar(hgnc_symbol)>0) %>% select(-ensembl_gene_id) %>% 
  distinct(hgnc_symbol, .keep_all = T) %>% column_to_rownames("hgnc_symbol")
human.stim.bulk <- inner_join(genes.meta[-3],
                              human.stim.bulk.r %>% rename(ensembl_gene_id = id)) %>%
  filter(nchar(hgnc_symbol)>0) %>% select(-ensembl_gene_id) %>% 
  distinct(hgnc_symbol, .keep_all = T) %>% column_to_rownames("hgnc_symbol")

mouse.bulk <- mouse.bulk.r %>% 
  distinct(symbol, .keep_all = T) %>%
  column_to_rownames("symbol")

rm(human.nostim.bulk.r);rm(human.stim.bulk.r);rm(mouse.bulk.r);gc()
################################################################################
################################################################################
# normalize data
library(preprocessCore)
human.nostim.bulk.norm <- normalize.quantiles(human.nostim.bulk %>% as.matrix()) %>% as.data.frame()
rownames(human.nostim.bulk.norm) <- rownames(human.nostim.bulk)
colnames(human.nostim.bulk.norm) <- colnames(human.nostim.bulk)

human.stim.bulk.norm <- normalize.quantiles(human.stim.bulk %>% as.matrix()) %>% as.data.frame()
rownames(human.stim.bulk.norm) <- rownames(human.stim.bulk)
colnames(human.stim.bulk.norm) <- colnames(human.stim.bulk)

mouse.bulk.norm <- normalize.quantiles(mouse.bulk %>% as.matrix()) %>% as.data.frame()
rownames(mouse.bulk.norm) <- rownames(mouse.bulk)
colnames(mouse.bulk.norm) <- colnames(mouse.bulk)

# ref

## sum by celltype_cluster?
allen.human.cl <- allen.human %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -c("gene"), names_to = "subclass_label") %>%
  mutate(subclass_label = str_replace_all(string = subclass_label, pattern = "\\.", replacement = " ")) %>%
  left_join(celltypes.meta.human %>% mutate(subclass_label = str_replace_all(string = subclass_label,
                                                                             pattern = "\\/", replacement = " "),
                                            subclass_label = str_replace_all(string = subclass_label,
                                                                             pattern = "-", replacement = " "))) %>%
  group_by(gene, cluster_label_v2) %>%
  dplyr::summarise(value = sum(value)) %>% 
  pivot_wider(names_from = cluster_label_v2, id_cols = gene) %>%
  column_to_rownames("gene")
allen.mouse.cl <- allen.mouse %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -c("gene"), names_to = "subclass_label") %>%
  mutate(subclass_label = str_replace_all(string = subclass_label, pattern = "\\.", replacement = " ")) %>%
  left_join(celltypes.meta.mouse %>% mutate(subclass_label = str_replace_all(string = subclass_label,
                                                                             pattern = "\\/", replacement = " "),
                                            subclass_label = str_replace_all(string = subclass_label,
                                                                             pattern = "-", replacement = " "))) %>%
  group_by(gene, cluster_label_v2) %>%
  dplyr::summarise(value = sum(value)) %>% 
  pivot_wider(names_from = cluster_label_v2, id_cols = gene) %>%
  column_to_rownames("gene")



allen.human.norm <- normalize.quantiles(allen.human.cl %>% as.matrix()) %>% as.data.frame()
rownames(allen.human.norm) <- rownames(allen.human.cl)
colnames(allen.human.norm) <- colnames(allen.human.cl)

allen.mouse.norm <- normalize.quantiles(allen.mouse.cl %>% as.matrix()) %>% as.data.frame()
rownames(allen.mouse.norm) <- rownames(allen.mouse.cl)
colnames(allen.mouse.norm) <- colnames(allen.mouse.cl)

rm(human.nostim.bulk);rm(human.stim.bulk);rm(mouse.bulk);rm(allen.human);rm(allen.mouse);gc()
################################################################################
################################################################################
################################################################################
# run deconvolution
human.stim.deco <- deconvolute(human.stim.bulk.norm %>% as.matrix(), 
                               sigMatrix = allen.human.norm %>% as.matrix(), 
                               methods = c("dtangle", "nnls", "rls"), 
                               use_cores = 20)
human.nostim.deco <- deconvolute(human.nostim.bulk.norm %>% as.matrix(), 
                                 sigMatrix = allen.human.norm %>% as.matrix(), 
                                 methods = c("dtangle", "nnls", "rls"), 
                                 use_cores = 20)
mouse.deco <- deconvolute(mouse.bulk.norm %>% as.matrix(), 
                          sigMatrix = allen.mouse.norm %>% as.matrix(), 
                          methods = c("dtangle", "nnls", "rls"), 
                          use_cores = 20)

save(list = c("human.nostim.deco", "human.stim.deco", "mouse.deco"), file = "data/derivatives/deconvolution-data.rda")
load("data/derivatives/deconvolution-data.rda")
p1 <- human.stim.deco$proportions$dtangle_sig1 %>%
  as.data.frame() %>%
  rownames_to_column("name") %>%
  pivot_longer(cols = -c("name"), names_to = "cluster_label_v2", values_to = "proportion") %>%
  mutate(stimulation = ifelse(grepl("Baseline", name), "baseline", "stim"),
         sample = sub("BulkRNA", "", sub("Baseline|Stimulated", "", name)),
         sample = sub("Male|Female", "", sample),
         sample = sub("subject", "sample", sample)) %>%
  ggplot(aes(x = stimulation, y = proportion, fill = cluster_label_v2)) +
  geom_bar(stat = "identity") +
  ggh4x::facet_grid2(cols = vars(sample)) +
  scale_fill_manual(values = cell.colors, name = "cluster label") +
  labs(x = "", title = "bulk RNA-Seq human samples (stimulation)") +
  bw.theme
p1

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# try deconvoluting the bulk DE results
human.stim.bulk.r.degs <- readxl::read_xlsx("data/YV/Supplementary data 2_Bulk_Stimulated_human_10.01.2024.xlsx")%>%
  select(gene = 2, logFC = 3) %>% distinct(gene,.keep_all = T)
mouse.stim.bulk.r.degs <- readxl::read_xlsx("data/YV/Supplementary data 4_Bulk_Stimulated_mouse_10.01.2024.xlsx")%>%
  select(gene = 2, logFC = 3)

human.stim.sn.degs <- read_csv("data/derivatives/gex-DE/final/all-results.csv") %>%
  filter(dx == "combined") %>%
  select(gene, coef = lmmSeq_coef_stimulation, RNA_predicted_cluster) %>%
  pivot_wider(names_from = RNA_predicted_cluster, values_from = coef, id_cols = gene)

# intersecting genes
int.g <- inner_join(human.stim.bulk.r.degs %>% select(gene),
                    human.stim.sn.degs %>% select(gene))
human.stim.sn.degs <- human.stim.sn.degs %>% filter(gene %in% int.g$gene)
human.stim.bulk.r.degs <- human.stim.bulk.r.degs %>% filter(gene %in% int.g$gene)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
#### try a glm to predict the logFC of bulk using sn
df <- human.stim.bulk.r.degs %>% inner_join(human.stim.sn.degs)
df1 <- jtools::summ(glm(logFC ~ ., data = df %>% select(-gene)), confin = T, pval = T)$coeftable %>%
  as.data.frame() %>% rownames_to_column("x") %>%
  rename(Estimate = 2, confin_min = 3, confin_max = 4, pval = 6)
p004_2<- df1 %>%
  ggplot(aes(x=Estimate, y = x, color = x)) +
  geom_point(show.legend = F) +
  geom_errorbarh(aes(xmin=confin_min, xmax=confin_max), show.legend = F,
                 height = 0.1) +
  scale_color_manual(values = c(cell.colors, "grey")) +
  geom_vline(xintercept = 0, linetype = 2, color = "pink") +
  bw.theme + labs(y="", caption = paste0("glm(bulk-log2FC ~ stimulation coefficients per cell cluster)"))
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
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# sum the baseline samples

stim.cl.diff <- p1$data %>%
  group_by(cluster_label_v2, stimulation) %>%
  dplyr::summarise(prop_sum = sum(proportion)) %>%
  ungroup() %>% group_by(stimulation) %>%
  mutate(proportion = (prop_sum/sum(prop_sum))*100) %>%
  mutate(stimulation = case_when(stimulation == "stim" ~ "stimulation", stimulation == "baseline" ~ "baseline")) %>%
  pivot_wider(names_from = stimulation, values_from = proportion, id_cols = cluster_label_v2) %>%
  mutate(diff = stimulation - baseline,
         avg = (stimulation+baseline)/2) 

tt <- p1$data %>%
  mutate(stimulation = case_when(stimulation == "stim" ~ "stimulation", stimulation == "baseline" ~ "baseline")) %>%
  pivot_wider(names_from = stimulation, values_from = proportion, id_cols = c(cluster_label_v2, sample)) %>%
  mutate(diff = stimulation - baseline)
tt.summ <- tt %>%
  group_by(cluster_label_v2) %>%
  summarise(mean_diff = mean(diff), 
            sem_diff = sd(diff) / sqrt(n()))
p003 <- tt.summ %>%
  ggplot(aes(x = cluster_label_v2, y = mean_diff, fill = cluster_label_v2)) +
  geom_bar(stat = "identity", show.legend = F) +
  geom_errorbar(aes(ymin = mean_diff-sem_diff, ymax = mean_diff+sem_diff),
                width = 0.2, alpha = 0.4) +
  ylim(c(-50,50)) +
  geom_text(aes(x = "Micro", y = 30, 
                label = "more presence in stimulation samples")) +
  geom_text(aes(x = "Micro", y = -30, 
                label = "more presence in baseline samples")) +
  geom_hline(yintercept = c(-5,0,5), color = redblu.col[1], linetype = 2)+
  # geom_text(aes(x = 0, y = 5, label = "5%"), size = 3, nudge_x = 0.7) +
  # geom_text(aes(x = 0, y = -5, label = "-5%"), size = 3, nudge_x = 0.7) +
  scale_fill_manual(values = cell.colors, name = "cluster label") +
  labs(x = "", y = "proportion difference\n(stimulation - baseline)",
       title = "stimulation-based deconvolution differences",
       subtitle = "human bulk RNA-Seq") +
  bw.theme

patchwork::wrap_plots(p003,p004_2,nrow = 1)
ggsave2("figs/keep/deconv-bulk_v3.pdf", width = 10, height = 4.5)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
