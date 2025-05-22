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
source("src/98_manuscript-color-palettes.R")
################################################################################
################################################################################
################################################################################
# read the DEGs from the analysis 
ns.res <- read_csv("data/derivatives/gex-DE/final/all-results.csv") %>% filter(dx=="combined")
################################################################################
################################################################################
################################################################################
################################################################################
# data from epilepsy neuronal sub-types paper
# https://www.nature.com/articles/s41467-020-18752-7#Sec40
# Pfisterer et al.

lit.epi.genes <- readxl::read_xlsx("data/raw/downloaded/epi-neuronal-subtypes-supp.xlsx",
                                   sheet = "S7", skip = 2) %>%
  select(Gene) %>% drop_na()
## the sc analysis
epi.neuro.inh <- readxl::read_xlsx("data/raw/downloaded/epi-neuronal-subtypes-supp.xlsx",
                                   sheet = "S6", range = "A4:I3759") %>%
  mutate_at(.vars = vars(c(3:8)), .funs = function(x) as.numeric(x)) %>%
  dplyr::filter(abs(log2FoldChange)>1, padj<0.05) %>% distinct(Gene)
epi.neuro.exc <- readxl::read_xlsx("data/raw/downloaded/epi-neuronal-subtypes-supp.xlsx",
                                   sheet = "S6", range = "L4:T32092") %>%
  mutate_at(.vars = vars(c(3:8)), .funs = function(x) as.numeric(x)) %>%
  dplyr::filter(abs(log2FoldChange)>1, padj<0.05) %>% distinct(Gene)
################################################################################
# epilepsy gwas
# Genome-wide mega-analysis identifies 16 loci and highlights diverse biological mechanisms in the common epilepsies. Nat. Commun. (2018)
# https://www.nature.com/articles/s41467-018-07524-z
epi.gwas.genes <- c("ATXN1","BCL11A","BRD7","C3orf33","FANCL","GABRA2",
                    "GJA1","GRIK1","HEATR3","KCNAB1","KCNN2","PCDH7",
                    "PNPO","SCN1A","SCN2A","SCN3A","STAT4","STX1B",
                    "TTC21B","ZEB2")
################################################################################
# data from the Allen institute for neurosurgery activated or PM genes
# https://pmc.ncbi.nlm.nih.gov/articles/PMC6919571/#ABS1
ns.pm.genes <- readxl::read_xlsx("data/raw/downloaded/allen-ns-vs-pm.xlsx") %>%
  rename(logFC = 2, adj_pval = 6)
################################################################################
# data from the epilepsy paper from Leo and Jake
epi.spiking.genes <- readxl::read_xlsx("data/raw/downloaded/Leo-apilepsy-paper.xlsx")
################################################################################
# konopka papers
# mem biomarkers
mem.genes <- readxl::read_xlsx("data/raw/downloaded/berto-konopka-mem-biomarkers-nn-2021.xlsx",
                               sheet = "SME Multivariate Stats")
mem.genes.2 <- readxl::read_xlsx("data/raw/downloaded/berto-konopka-mem-biomarkers-cc-2018.xlsx",
                                 skip = 3) %>%
  filter(Gene %in% ns.res$gene)
################################################################################
################################################################################
## Glioblastoma
## Science Advances paper
# https://www.science.org/doi/10.1126/sciadv.adn4306#sec-4
wang.gbm <- readxl::read_xlsx("data/raw/downloaded/Wang-et-al-GBM-DEGs-SA.xlsx", sheet = "Table S4") %>%
  dplyr::filter(abs(avg_log2FC)>1, Cluster %in% c("GBM cell", "Astrocyte", "Oligodendrocyte", "Inhibitory neuron",
                                                  "Excitatory neuron"))
################################################################################
################################################################################
################################################################################
################################################################################
all.lit.genes <- data.frame(gene = unique(c(lit.epi.genes$Gene,
                                            epi.gwas.genes, epi.neuro.exc$Gene, epi.neuro.inh$Gene,
                                            ns.pm.genes$gene,
                                            wang.gbm$Gene, 
                                            epi.spiking.genes$symbol,
                                            mem.genes$Gene, mem.genes.2$Gene,
                                            ns.res$gene))) %>%
  mutate(EPI_lit = (gene %in% lit.epi.genes$Gene),
         EPI_GWAS = (gene %in% epi.gwas.genes),
         EPI_Brueggeman = (gene %in% epi.spiking.genes$symbol),
         EPI_Pfisterer_inh = (gene %in% epi.neuro.inh$Gene),
         EPI_Pfisterer_exc = (gene %in% epi.neuro.exc$Gene),
         GBM_Wang_GBM = (gene %in% wang.gbm$Gene[wang.gbm$Cluster=="GBM cell"]),
         GBM_Wang_Astro = (gene %in% wang.gbm$Gene[wang.gbm$Cluster=="Astrocyte"]),
         GBM_Wang_Oligo = (gene %in% wang.gbm$Gene[wang.gbm$Cluster=="Oligodendrocyte"]),
         GBM_Wang_Inh = (gene %in% wang.gbm$Gene[wang.gbm$Cluster=="Inhibitory neuron"]),
         GBM_Wang_Exc = (gene %in% wang.gbm$Gene[wang.gbm$Cluster=="Excitatory neuron"]),
         MEM_Berto_NN_2021 = (gene %in% mem.genes$Gene),
         MEM_Berto_CC_2018 = (gene %in% mem.genes.2$Gene),
         Allen_neurosurgery_up = (gene %in% ns.pm.genes$gene[ns.pm.genes$logFC>0]),
         Allen_PM_up = (gene %in% ns.pm.genes$gene[ns.pm.genes$logFC<0]))
################################################################################
################################################################################
################################################################################
## plot
intersection.genes <- ns.res %>% 
  select(RNA_predicted_cluster, gene, starts_with("lmmSeq_DEG")) %>%
  left_join(all.lit.genes) %>%
  pivot_longer(cols = colnames(all.lit.genes)[-1]) %>%
  filter(value==T) %>%
  mutate(category = case_when(grepl("MEM", name) ~ "Memory",
                              grepl("EPI", name) ~ "Epilepsy",
                              grepl("GBM", name) ~ "Glioblastoma",
                              grepl("Allen", name) ~ "Allen"),
         name = case_when(name == "EPI_Brueggeman" ~ "Brueggeman et al.",
                          name == "EPI_lit" ~ "literature",
                          name == "EPI_GWAS" ~ "GWAS",
                          name == "EPI_Pfisterer_inh" ~ "Pfisterer et al. Inhibitory neurons",
                          name == "EPI_Pfisterer_exc" ~ "Pfisterer et al. Excitatory neurons",
                          name == "GBM_Wang_GBM" ~ "Wang et al. GBM cells",
                          name == "GBM_Wang_Astro" ~ "Wang et al. Astrocytes",
                          name == "GBM_Wang_Oligo" ~ "Wang et al. Oligodendrocytes",
                          name == "GBM_Wang_Inh" ~ "Wang et al. Inhibitory neurons",
                          name == "GBM_Wang_Exc" ~ "Wang et al. Excitatory neurons",
                          name == "MEM_Berto_NN_2021" ~ "Berto et al. NN",
                          name == "MEM_Berto_CC_2018" ~ "Berto et al. CC",
                          name == "Allen_neurosurgery_up" ~ "neurosurgery upregulated",
                          name == "Allen_PM_up" ~ "postmortem upregulated"))
write_rds(intersection.genes, "data/derivatives/all-lit-genes.rds")

thresholds <- c("cons")
foreach(t = 1:length(thresholds)) %dopar% {
  
  p01 <- intersection.genes %>% 
    rename(sig_both = paste0("lmmSeq_DEG_", thresholds[t])) %>%
    filter(sig_both==T) %>%
    ggplot(aes(x = name, y = gene, fill = sig_both)) +
    geom_tile(show.legend = F) +
    scale_fill_manual(values = abstract.colors.2[1], name = "stimulation") +
    ggh4x::facet_grid2(rows = vars(RNA_predicted_cluster), cols = vars(category),
                       scales = "free",
                       space = "free") +
    labs(y = "", x = "",
         caption = paste0("only showing intersecting DEGs.")) +
    bw.theme +
    theme(axis.text.x.bottom = element_text(angle = 65, hjust = 1))
  
  
  p02 <- intersection.genes %>% 
    rename(sig_both = paste0("lmmSeq_DEG_", thresholds[t])) %>%
    filter(sig_both==T) %>%
    ggplot(aes(x = name)) +
    geom_bar(fill = antique.colors[5]) +
    geom_text(stat = "count", aes(label = ..count..), 
              position = position_stack(vjust = 0.5),
              color = "white") +
    ggh4x::facet_grid2(rows = vars(RNA_predicted_cluster), cols = vars(category),
                       scales = "free",
                       space = "free") +
    bw.theme +
    theme(axis.text.x.bottom = element_text(angle = 65, hjust = 1)) +
    labs(x = "Literature Gene Set", y = "count", 
         title = "Intersecting Genes")
  
  patchwork::wrap_plots(p02,p01)
  ggsave2(paste0("figs/keep/stimultion-DEGs-vs-literature_", thresholds[t], "-thr.pdf"),
          width = 10, height = ifelse(t==1,22,34)) 
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
registerDoMC(cores = 4)
fisher.res <- foreach(c = 1:length(unique(ns.res$RNA_predicted_cluster)), .combine = rbind) %dopar% {
  cl <- unique(ns.res$RNA_predicted_cluster)[c]
  df <- ns.res %>% filter(RNA_predicted_cluster==cl) %>%
    select(gene, eSTIM = lmmSeq_DEG_cons) %>%
    inner_join(all.lit.genes)
  l.r <- do.call(rbind,
                 lapply(df[,-c(1:2)], function(x) {
                   ft <- fisher.test(table(df$eSTIM,x))
                   data.frame(RNA_predicted_cluster = cl,
                              # dataset = colnames(x),
                              OR = as.numeric(ft$estimate),
                              pval = ft$p.value,
                              confin_min = ft$conf.int[1],
                              confin_max = ft$conf.int[2])
                 })) %>%
    rownames_to_column("dataset")
  l.r
}
fisher.res.2 <- fisher.res %>% 
  mutate(name = dataset,
         name = case_when(name == "EPI_Brueggeman" ~ "Brueggeman et al.",
                          name == "EPI_lit" ~ "literature",
                          name == "EPI_GWAS" ~ "GWAS",
                          name == "EPI_Pfisterer_inh" ~ "Pfisterer et al. Inhibitory neurons",
                          name == "EPI_Pfisterer_exc" ~ "Pfisterer et al. Excitatory neurons",
                          name == "GBM_Wang_GBM" ~ "Wang et al. GBM cells",
                          name == "GBM_Wang_Astro" ~ "Wang et al. Astrocytes",
                          name == "GBM_Wang_Oligo" ~ "Wang et al. Oligodendrocytes",
                          name == "GBM_Wang_Inh" ~ "Wang et al. Inhibitory neurons",
                          name == "GBM_Wang_Exc" ~ "Wang et al. Excitatory neurons",
                          name == "MEM_Berto_NN_2021" ~ "Berto et al. NN",
                          name == "MEM_Berto_CC_2018" ~ "Berto et al. CC",
                          name == "Allen_neurosurgery_up" ~ "neurosurgery upregulated",
                          name == "Allen_PM_up" ~ "postmortem upregulated")) %>%
  left_join(intersection.genes %>% distinct(name, category))
fisher.res.2 %>%
  ggplot(aes(x = name, y = RNA_predicted_cluster, fill = OR, 
             label = ifelse(pval < 0.05, "*",""))) +
  geom_tile() + 
  scale_fill_gradient2(low = antique.colors[5],high = antique.colors[5], name = "OR") + my.guides + 
  geom_text() +
  ggh4x::facet_grid2(cols = vars(category),scales = "free",space = "free") +
  bw.theme + theme(axis.text.x.bottom = element_text(angle = 65, hjust = 1)) +
  labs(x = "Literature Gene Set")

fisher.res.2 %>% 
  ggplot(aes(x = log2(OR), y = -log10(pval), color = RNA_predicted_cluster,
             label = ifelse(pval<0.05,name,""))) +
  geom_hline(yintercept = -log10(0.05),color="red",linetype=2) +
  geom_vline(xintercept = log2(1),color="red",linetype=2) +
  geom_point() + ggrepel::geom_text_repel(size = 3, show.legend = F) +
  scale_color_manual(values = cell.colors, name="cluster") +
  ggh4x::facet_grid2(cols = vars(category),scales = "free",space = "free") +
  bw.theme + theme(axis.text.x.bottom = element_text(angle = 65, hjust = 1)) +
  labs(x = expression('log'[2]~'(OR)'),
       y = expression('-log'[10]~'(p-value)'))
ggsave2("figs/keep/stimultion-DEGs-vs-literature_cons-thr_fisher.pdf",width = 12, height = 6)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################