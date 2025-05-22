################################################################################
################################################################################
rm(list = ls()); gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
              "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
################################################################################
################################################################################
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/neurosurgery")
setwd(project.dir)
source("src/98_manuscript-color-palettes.R")
################################################################################
################################################################################
################################################################################
# read the DEGs from the analysis 
ns.res <- read_csv("data/derivatives/gex-DE/final/all-results.csv") %>% dplyr::filter(dx=="combined")
################################################################################
################################################################################
################################################################################
library(fgsea)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(goseq)
library(clusterProfiler)

genes.meta <- data.frame(gene = unique(ns.res$gene))
genes.meta <- inner_join(genes.meta,
                         mapIds(org.Hs.eg.db, keys = genes.meta$gene,
                                column = "ENTREZID", keytype = "SYMBOL") %>%
                           as.data.frame() %>%
                           rownames_to_column("gene") %>%
                           dplyr::rename("entrezid" = 2))
################################################################################
################################################################################
registerDoMC(cores = 3)
## select method and threshold type
thresholds <- c("cons")
methods <- c("lmmSeq_DEG")
go.res <- foreach(cc = 1:length(unique(ns.res$RNA_predicted_cluster)), .combine = rbind) %dopar% {
  cl <- unique(ns.res$RNA_predicted_cluster)[cc]
  df <- ns.res %>%
    dplyr::filter(RNA_predicted_cluster == cl) %>%
    left_join(genes.meta) 
  
  rres <- foreach(t = 1:length(thresholds), .combine = rbind) %dopar% {
    thr.name <- thresholds[t]
    df2 <- df %>% dplyr::select(gene, paste0(methods, "_", thr.name))
    
    mres <- foreach(m = 2:4, .combine = rbind) %dopar% {
      df3 <- df2[,c(1,m)] %>% dplyr::rename(DEG = 2)
      # GO 
      DEGs.g <- as.integer(df3$DEG)
      names(DEGs.g) <- df3$gene
      pwf <- nullp(DEGs.g, genome = "hg19", id = "geneSymbol")
      goRes <- goseq(pwf, genome = "hg19", id = "geneSymbol") %>%
        group_by(ontology) %>%mutate(FDR = p.adjust(over_represented_pvalue, method = "fdr")) %>%ungroup()
      
      return(goRes %>%
               mutate(RNA_predicted_cluster = cl,
                      type = colnames(df2)[m], thr = thr.name))
    }
    return(mres)
  }
  return(rres)
}
write_csv(go.res, "data/derivatives/enrichment-analysis/goseq-GO-all.csv")
library(ggbreak)


go.res %>%
  filter(FDR < 0.05) %>%
  mutate(type = sub("_.*", "", type),
         type = ifelse(type=="sig", "both methods", type)) %>% 
  filter(type == "lmmSeq",thr=="cons") %>%
  group_by(ontology, RNA_predicted_cluster, thr) %>%
  top_n(n = 10, wt = -FDR) %>% ungroup() %>%
  mutate(hitsPerc=numDEInCat*100/numInCat, 
         thr = factor(thr, levels = c("cons", "rel5", "rel2", "emp"))) %>%
  ggplot(aes(x = numDEInCat, y = -log10(FDR),
             label = term, color = RNA_predicted_cluster)) +
  geom_point(size = 0.6,position = "jitter") +
  ggrepel::geom_text_repel(direction = "y", nudge_y = -0.01, 
                           show.legend = F, size = 1.5, max.overlaps = 10) +
  scale_color_manual(values = cell.colors) +
  ggh4x::facet_grid2(rows=  vars(ontology), 
                     scales = "free", independent = T) +
  labs(x="count of DEGs in category",  y=expression("-log"[10]~"FDR"),
       color = "cluster") +
  bw.theme
ggsave2("figs/keep/GO-analysis.pdf", height = 8, width = 5)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
library(ActivePathways)
scores <- ns.res %>% dplyr::select(RNA_predicted_cluster, gene, lmmSeq_FDR) %>%
  pivot_wider(names_from = RNA_predicted_cluster, id_cols = gene, values_from = lmmSeq_FDR) %>%
  column_to_rownames("gene")
scores <- as.matrix(scores); scores[is.na(scores)] <- 1
registerDoMC(cores = 4)
thr <- 0.1
enriched_pathways <- foreach(cc =1:ncol(scores), .combine = rbind) %dopar% {
  ranks <- scores[,cc]%>%as.data.frame()%>%rownames_to_column("gene")%>%
    dplyr::rename(rank=2)%>%dplyr::arrange(rank) %>% mutate(rank = -log10(rank))
  write_tsv(ranks,
              paste0("data/derivatives/enrichment-analysis/ActivePathways/", colnames(scores)[cc], "_gene-ranks.txt"),
              col_names = F)
  df <- ActivePathways(scores[scores[,cc]<thr,cc]%>%as.matrix(), 
                       gmt = correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/data/genomics/gene-sets/Human_Reactome_April_01_2025_symbol.gmt"), 
                       cytoscape_file_tag = paste0(project.dir, "/data/derivatives/enrichment-analysis/ActivePathways/", colnames(scores)[cc], "_"))
  ifelse(nrow(df)>0, return(df %>% mutate(RNA_predicted_cluster = colnames(scores)[cc])), return(NULL))
}
write_csv(enriched_pathways, "data/derivatives/enrichment-analysis/ActivePathways/combined-results.csv")

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################