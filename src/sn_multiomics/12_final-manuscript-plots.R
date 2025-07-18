################################################################################
################################################################################
rm(list = ls()); gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
              "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
library(patchwork)
################################################################################
################################################################################
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/neurosurgery")
setwd(project.dir)
source("src/98_manuscript-color-palettes.R")
################################################################################
################################################################################
int.samples.so <- read_rds("data/derivatives/ME-processing/processed-seurat-objects/integrated-multiome-samples.rds")
samples.meta <- read_csv("data/samples-metadata_ME.csv")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## UMAP
int.metadat <- read_rds("data/derivatives/cells-metadata-filtered.rds") %>%
  filter(!RNA_predicted_cluster %in% c("Endo", "VLMC"), sampleID_2 != "P19")

p.umap <- int.metadat %>%
  ggplot(aes(x=umap_1, y = umap_2, color = RNA_predicted_cluster)) +
  geom_point(show.legend = F,size = 0.5) + scale_color_manual(values = cell.colors.2) +
  annotate(geom = "text", x=-8,y=10, label = "Exc", color = cell.colors.2[["Exc"]], fontface="bold")+
  annotate(geom = "text", x=-4,y=1, label = "Inh", color = cell.colors.2[["Inh"]], fontface="bold")+
  annotate(geom = "text", x=10,y=-3, label = "Micro", color = cell.colors.2[["Micro"]], fontface="bold")+
  annotate(geom = "text", x=5,y=10, label = "Oligo", color = cell.colors.2[["Oligo"]], fontface="bold") +
  annotate(geom = "text", x=3,y=-9, label = "OPC", color = cell.colors.2[["OPC"]], fontface="bold")+
  annotate(geom = "text", x=-4,y=-6, label = "Astro", color = cell.colors.2[["Astro"]], fontface="bold")+
  bw.theme +
  labs(x = "UMAP-1", y="UMAP-2")
pdf("figs/keep/aligned/01_umap__V2.pdf",height = 5,width = 5)
p.umap
dev.off()
################################################################################
################################################################################
################################################################################
## lmmSeq DE volcano plot
## there's an outlier gene in the Exc
ns.de <- read_csv("data/derivatives/gex-DE/final/all-results.csv") %>% filter(dx=="combined", -log10(lmmSeq_pval) < 100) %>%
  mutate(final_DEG = lmmSeq_DEG_cons)
conts <- ns.de %>% mutate(dir = case_when(final_DEG & lmmSeq_coef_stimulation>0~ "Up", 
                                          final_DEG & lmmSeq_coef_stimulation<0 ~"Down")) %>%
  group_by(RNA_predicted_cluster, dir) %>% 
  dplyr::summarise(c=n(),y = -log10(min(lmmSeq_pval))) %>% ungroup() %>%
  rbind(data.frame(RNA_predicted_cluster="Exc",dir="Down",c=0,y=0)) %>%
  group_by(RNA_predicted_cluster) %>% dplyr::mutate(y=max(y)) %>%
  filter(!is.na(dir)) %>% mutate(x = case_when(dir=="Up"&RNA_predicted_cluster%in%c("Exc","Inh") ~ 1.5, 
                                               dir=="Up" ~ 3.3, 
                                               dir=="Down"&RNA_predicted_cluster%in%c("Exc","Inh") ~ -1.7,
                                               dir=="Down" ~ -2.2), 
                                 y=y-2,
                                 y2 = ifelse(RNA_predicted_cluster%in%c("Exc","Inh"), 30,40))



p.volc.lmmseq.n <- ns.de %>% 
  filter(RNA_predicted_cluster %in% c("Exc","Inh")) %>%
  mutate(siggg = abs(lmmSeq_coef_stimulation)>1&lmmSeq_FDR<0.05) %>%
  ggplot(aes(x = lmmSeq_coef_stimulation, y = -log10(lmmSeq_pval))) +
  geom_point(aes(color = ifelse(siggg, RNA_predicted_cluster, "grey"),
                 alpha = ifelse(siggg, 1, 0.1),
                 size = siggg), 
             show.legend = F) + 
  scale_color_manual(values = c(cell.colors.2, "grey")) +
  ggnewscale::new_scale_color() + 
  ylim(c(0,32)) + xlim(c(-2,2)) +
  geom_vline(xintercept = c(1,-1), linetype=2, color = "pink") +geom_vline(xintercept = c(0), linetype=1, color = "black",alpha=0.2) +
  geom_hline(yintercept = -log10(0.01), linetype=2, color = "pink") +
  ggrepel::geom_text_repel(aes(label = ifelse(siggg,gene,"")),
                           size=2, max.time = 0.1, max.overlaps = 30) +
  geom_text(aes(x = x, y = y2,label = paste0(dir,": ", c),color=dir),
            data = conts %>% filter(RNA_predicted_cluster %in% c("Exc","Inh")), 
            size = 3,show.legend = F,fontface = "bold") +
  ggh4x::facet_grid2(rows = vars(RNA_predicted_cluster), axes = "all",scales = "free_y",
                     strip = ggh4x::strip_themed(text_y = list(element_text(color = cell.colors.2[["Inh"]],size=12,face = "bold"),
                                                               element_text(color = cell.colors.2[["Exc"]],size=12,face = "bold")))) + 
  scale_color_manual(values = c(redblu.col[c(2,1)]))+
  scale_size_manual(values = c(0.3,1.2)) +
  labs(x = "", y = expression("-log"[10]~"(p-value)")) +
  bw.theme

p.volc.lmmseq.nn <- ns.de %>% 
  filter(!RNA_predicted_cluster %in% c("Exc","Inh")) %>%
  mutate(siggg = abs(lmmSeq_coef_stimulation)>1&lmmSeq_FDR<0.05) %>%
  ggplot(aes(x = lmmSeq_coef_stimulation, y = -log10(lmmSeq_pval))) +
  geom_point(aes(color = ifelse(siggg, RNA_predicted_cluster, "grey"),
                 alpha = ifelse(siggg, 1, 0.1),
                 size = siggg), 
             show.legend = F) + 
  scale_color_manual(values = c(cell.colors.2, "grey")) +
  ggnewscale::new_scale_color() + 
  # ylim(c(0,32)) + xlim(c(-3.5,4)) +
  geom_vline(xintercept = c(1,-1), linetype=2, color = "pink") +geom_vline(xintercept = c(0), linetype=1, color = "black",alpha=0.2) +
  geom_hline(yintercept = -log10(0.01), linetype=2, color = "pink") +
  ggrepel::geom_text_repel(aes(label = ifelse(siggg,gene,"")),
                           size=2, max.time = 0.1, max.overlaps = 30) +
  geom_text(aes(x = x, y = y,label = paste0(dir,": ", c),color=dir),
            data = conts %>% filter(!RNA_predicted_cluster %in% c("Exc","Inh")), 
            size = 3,show.legend = F,fontface = "bold") +
  ggh4x::facet_grid2(rows = vars(RNA_predicted_cluster), axes = "all",scales = "free_y",
                     strip = ggh4x::strip_themed(text_y = list(element_text(color = cell.colors.2[["Astro"]],size=12,face = "bold"),
                                                               element_text(color = cell.colors.2[["Micro"]],size=12,face = "bold"),
                                                               element_text(color = cell.colors.2[["Oligo"]],size=12,face = "bold"),
                                                               element_text(color = cell.colors.2[["OPC"]],size=12,face = "bold")))) +

  scale_color_manual(values = c(redblu.col[c(2,1)]))+
  scale_size_manual(values = c(0.3,1.2)) +
  labs(x = "stimulation coefficient", y = expression("-log"[10]~"(p-value)")) +
  bw.theme


# p.volc.lmmseq <- ns.de %>% mutate(siggg = abs(lmmSeq_coef_stimulation)>1&lmmSeq_FDR<0.05,
#                                   dir = case_when(final_DEG & lmmSeq_coef_stimulation>0~ "upregulated", 
#                                                   final_DEG & lmmSeq_coef_stimulation<0 ~"downregulated")) %>%
#   ggplot(aes(x = lmmSeq_coef_stimulation, y = -log10(lmmSeq_pval))) +
#   geom_point(aes(color = ifelse(siggg, RNA_predicted_cluster, "grey"),
#                  alpha = ifelse(siggg, 1, 0.1),
#                  size = siggg), 
#              show.legend = F) + 
#   scale_color_manual(values = c(cell.colors.2, "grey")) +
#   ggnewscale::new_scale_color() + 
#   geom_vline(xintercept = c(1,-1), linetype=2, color = "pink") +geom_vline(xintercept = c(0), linetype=1, color = "black",alpha=0.2) +
#   geom_hline(yintercept = -log10(0.01), linetype=2, color = "pink") +
#   # ggrepel::geom_text_repel(aes(label = ifelse(gene%in%genes.of.int& !is.na(dir),gene,"")), max.overlaps = 10,size=2) +
#   ggrepel::geom_text_repel(aes(label = ifelse(siggg,gene,"")),
#                            size=2, max.time = 0.1, max.overlaps = 30) +
#   geom_text(aes(x = x, y = y,label = paste0(dir,": ", c),color=dir),
#             data = conts, size = 3,show.legend = F,fontface = "bold") +
#   ggh4x::facet_grid2(rows = vars(RNA_predicted_cluster), scales = "free",
#                      independent = T,
#                      strip = ggh4x::strip_themed(text_y = list(element_text(color = cell.colors.2[["Astro"]],size=12,face = "bold"),
#                                                                element_text(color = cell.colors.2[["Inh"]],size=12,face = "bold"),
#                                                                element_text(color = cell.colors.2[["Exc"]],size=12,face = "bold"),
#                                                                element_text(color = cell.colors.2[["Micro"]],size=12,face = "bold"),
#                                                                element_text(color = cell.colors.2[["Oligo"]],size=12,face = "bold"),
#                                                                element_text(color = cell.colors.2[["OPC"]],size=12,face = "bold")))) + 
#   scale_color_manual(values = c(redblu.col[c(2,1)]))+
#   scale_size_manual(values = c(0.3,1.2)) +
#   labs(x = "stimulation coefficient", y = expression("-log"[10]~"(p-value)")) +
#   bw.theme

p.volc.lmmseq <- patchwork::wrap_plots(p.volc.lmmseq.n, p.volc.lmmseq.nn, heights = c(1,2))

pdf("figs/keep/aligned/02_lmmseq-GEX-volcano_0626.pdf", width = 5, height = 12)
p.volc.lmmseq
dev.off()
################################################################################
################################################################################
################################################################################
## GO results
go.res <- read_csv("data/derivatives/enrichment-analysis/goseq-GO-all.csv")
library(ggbreak)
p.go <- go.res %>%
  filter(FDR < 0.05, thr == "cons") %>%
  mutate(type = sub("_.*", "", type),
         type = ifelse(type=="sig", "both methods", type)) %>% filter(type == "lmmSeq") %>%
  group_by(ontology, RNA_predicted_cluster, thr) %>%
  top_n(n = 10, wt = -FDR) %>% ungroup() %>%
  mutate(hitsPerc=numDEInCat*100/numInCat, 
         thr = factor(thr, levels = c("cons", "rel5", "rel2", "emp"))) %>%
  ggplot(aes(x = hitsPerc, y = -log10(FDR),
             label = term, color = RNA_predicted_cluster)) +
  geom_point(size = 0.6,position = "jitter") +
  ggrepel::geom_text_repel(direction = "y", nudge_y = -0.01, 
                           show.legend = F, size = 2, max.overlaps = 20) +
  scale_color_manual(values = cell.colors.2) +
  # ylim(c(1,10)) + scale_y_break(c(7,10)) +
  ggh4x::facet_grid2(rows=  vars(ontology), scales = "free", independent = T) +
  labs(x="percentage of DEGs in category",  y=expression("-log"[10]~"(FDR)"),color = "cluster") +
  bw.theme
pdf("figs/keep/aligned/09_GO-scatter.pdf",height = 14,width = 8)
p.go
dev.off()
################################################################################
################################################################################
################################################################################
## neuroestimator
## find it in the neuroestimator script
################################################################################
################################################################################
################################################################################
################################################################################
#### DEGs intersection
library(UpSetR);library(ComplexUpset);library(ggupset)
# bulk.mouse <- readxl::read_xlsx("data/YV/Supplementary data 4_Bulk_Stimulated_mouse_10.01.2024.xlsx", range = "A1:G15399") %>%
#   rename("gene" = 2) %>% mutate(gene = toupper(gene)) %>% rename_at(.vars = -c(1:2), .funs = function(x) paste0(x, "_mouse")) %>%
#   mutate(sig = FDR_mouse<0.05&abs(logFC_mouse)>0.2)
# bulk.human <- readxl::read_xlsx("data/YV/Supplementary data 2_Bulk_Stimulated_human_10.01.2024.xlsx") %>%
#   rename("gene" = 2)  %>% rename_at(.vars = -c(1:2), .funs = function(x) paste0(x, "_human")) %>%
#   mutate(sig = FDR_human<0.05&abs(logFC_human)>0.2)

# inter.df <- data.frame(gene = unique(c(ns.de$gene, bulk.human$gene, bulk.mouse$gene))) %>%
#   mutate(Micro_human = as.numeric(gene %in% ns.de$gene[(ns.de$lmmSeq_DEG_cons)&ns.de$RNA_predicted_cluster=="Micro"]),
#          Astro_human = as.numeric(gene %in% ns.de$gene[(ns.de$lmmSeq_DEG_cons)&ns.de$RNA_predicted_cluster=="Astro"]),
#          Oligo_human = as.numeric(gene %in% ns.de$gene[(ns.de$lmmSeq_DEG_cons)&ns.de$RNA_predicted_cluster=="Oligo"]),
#          OPC_human = as.numeric(gene %in% ns.de$gene[(ns.de$lmmSeq_DEG_cons)&ns.de$RNA_predicted_cluster=="OPC"]),
#          Exc_human = as.numeric(gene %in% ns.de$gene[(ns.de$lmmSeq_DEG_cons)&ns.de$RNA_predicted_cluster=="Exc"]),
#          Inh_human = as.numeric(gene %in% ns.de$gene[(ns.de$lmmSeq_DEG_cons)&ns.de$RNA_predicted_cluster=="Inh"]),
#          Bulk_human = as.numeric(gene %in% bulk.human$gene[bulk.human$sig]),
#          Bulk_mouse = as.numeric(gene %in% bulk.mouse$gene[bulk.mouse$sig]))
inter.df <- data.frame(gene = unique(c(ns.de$gene))) %>%
  mutate(Micro = as.numeric(gene %in% ns.de$gene[(ns.de$lmmSeq_DEG_cons)&ns.de$RNA_predicted_cluster=="Micro"]),
         Astro = as.numeric(gene %in% ns.de$gene[(ns.de$lmmSeq_DEG_cons)&ns.de$RNA_predicted_cluster=="Astro"]),
         Oligo = as.numeric(gene %in% ns.de$gene[(ns.de$lmmSeq_DEG_cons)&ns.de$RNA_predicted_cluster=="Oligo"]),
         OPC = as.numeric(gene %in% ns.de$gene[(ns.de$lmmSeq_DEG_cons)&ns.de$RNA_predicted_cluster=="OPC"]),
         Exc = as.numeric(gene %in% ns.de$gene[(ns.de$lmmSeq_DEG_cons)&ns.de$RNA_predicted_cluster=="Exc"]),
         Inh = as.numeric(gene %in% ns.de$gene[(ns.de$lmmSeq_DEG_cons)&ns.de$RNA_predicted_cluster=="Inh"]))

inter.df.long <- inter.df %>%
  pivot_longer(cols = -1) %>% filter(value == 1) %>%
  pivot_wider(names_from = value, values_from = name) %>% rename("cl_list" = 2)
p.upset <- inter.df.long %>%
  ggplot(aes(x = cl_list)) +
  geom_bar(width = 0.7) +
  scale_x_upset() + 
  scale_color_manual(values = cell.colors.2) +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=0, size = 2) +
  geom_text(data = inter.df.long %>% group_by(cl_list) %>% summarize(gene_names = paste(gene, collapse = "\n"), gene_count = n()) %>% filter(gene_count<20),
            aes(x = cl_list, y = gene_count +4,
                label = gene_names), nudge_y = 3, angle = 0, hjust = 0.5, size = 1,vjust=0)+
  labs(x = "", y="count",
       caption = paste0("DEGs criteria:\n",
                        "    Bulk_human:    abs(logFC)>0.2 & FDR < 0.05\n",
                        "    Bulk_mouse:    abs(logFC)>0.2 & FDR < 0.05\n",
                        "    any sn-RNA-Seq_human:    abs(stimulation coefficient) > 1 & FDR < 0.05")) +
  bw.theme
upset(inter.df[rowSums(inter.df[,-1])>1,], name = "",names(cell.colors.2),
      queries=list(upset_query(set=names(cell.colors.2)[1], color=cell.colors.2[1], fill=cell.colors.2[1]),
                   upset_query(set=names(cell.colors.2)[2], color=cell.colors.2[2], fill=cell.colors.2[2]),
                   upset_query(set=names(cell.colors.2)[3], color=cell.colors.2[3], fill=cell.colors.2[3]),
                   upset_query(set=names(cell.colors.2)[4], color=cell.colors.2[4], fill=cell.colors.2[4]),
                   upset_query(set=names(cell.colors.2)[5], color=cell.colors.2[5], fill=cell.colors.2[5]),
                   upset_query(set=names(cell.colors.2)[6], color=cell.colors.2[6], fill=cell.colors.2[6]))) +
  bw.theme


inter.df.long %>%
  ggplot(aes(x = cl_list)) +
  geom_bar() +
  scale_x_upset() +
  axis_combmatrix(override_plotting_function = function(df) {
    df %>%
      mutate(# Create a factor indicating if a cluster is in this intersection, or "not observed"
        cluster_fill = case_when(map_lgl(labels_split, ~ "Astro" %in% .x) ~ "Astro",
                                 map_lgl(labels_split, ~ "Micro" %in% .x) ~ "Micro",
                                 map_lgl(labels_split, ~ "Oligo" %in% .x) ~ "Oligo",
                                 map_lgl(labels_split, ~ "OPC" %in% .x) ~ "OPC",
                                 map_lgl(labels_split, ~ "Exc" %in% .x) ~ "Exc",
                                 map_lgl(labels_split, ~ "Inh" %in% .x) ~ "Inh")) %>%
      ggplot(aes(x = at, y = single_label)) +
      geom_rect(aes(fill = index %% 2 == 0),
                ymin = df$index - 0.5, ymax = df$index + 0.5,
                xmin = 0, xmax = 1, show.legend = F) +
      geom_point(data = function(d) d[d$observed, , drop = FALSE],
                 aes(fill = cluster_fill), shape = 21, size = 4, 
                 color = "black", show.legend = F) +  # shape 21 for fill + outline
      geom_line(data = function(d) d[d$observed, , drop = FALSE],
                aes(group = labels), linewidth = 1.2, color = "grey30") +
      scale_fill_manual(values = c(
        `TRUE` = "white",
        `FALSE` = "#F7F7F7",
        cell.colors.2
        # "Other" = "grey80",
        # "not observed" = "lightgrey"
        )) +
      theme_void()
  })


pdf("figs/keep/aligned/05_upset.pdf", width = 9, height = 5)
p.upset
dev.off()
################################################################################
################################################################################
################################################################################
## TF
global.motifs <- read_rds("data/derivatives/atac-DA/global-motifs.rds")
prm.motifs <- read_rds("data/derivatives/atac-DA/prm-tes-motifs.rds")
motif.all <- global.motifs %>%
  mutate(source = "global") %>%
  rbind(prm.motifs %>% mutate(source = "PRM_TES"))
atac.pb.so <- read_rds("data/derivatives/manual-boostrapped-pseudobulk/archive_050625/no-weighting-50-or-100-100it-1/atac-peaks/motif-so.rds")
library(ggseqlogo)
p.motif <- Signac::MotifPlot(atac.pb.so, 
          motifs = unique(motif.all$motif[which(motif.all$motif.name %in% c("HINFP", "NRF1", "ZBTB33", 
                                                                            "ZBTB14", "TCFL5", "KLF15"))])) +
  scale_fill_manual(values = antique.colors) +
  facet_wrap(~seq_group, scales = "free_x") +
  bw.theme
pdf("figs/keep/aligned/06_TF-motifs.pdf", width = 7, height = 5)
p.motif;dev.off()


inter.df.motif <- motif.all %>% 
  filter(source=="global",p.adjust<0.05) %>% 
  pivot_wider(names_from = RNA_predicted_cluster, 
              values_from = p.adjust, id_cols = motif.name) %>% 
  mutate_at(.vars = vars(-c(motif.name)), .funs = function(x) as.numeric(ifelse(x<0.05,T,F))) 
inter.df.motif.long <- inter.df.motif %>%
  pivot_longer(cols = -1) %>% filter(value == 1) %>%
  pivot_wider(names_from = value, values_from = name) %>% rename("cl_list" = 2)
p.upset.motif <- inter.df.motif.long %>%
  ggplot(aes(x = cl_list)) +
  geom_bar(width = 0.7) +
  scale_x_upset() + 
  geom_text(stat='count', aes(label=after_stat(count)), vjust=0, size = 2) +
  geom_text(data = inter.df.motif.long %>% group_by(cl_list) %>% 
              summarize(motif_names = paste(motif.name, collapse = "\n"), 
                        motif_count = n()) %>% filter(motif_count<20),
            aes(x = cl_list, y = motif_count +4,
                label = motif_names), nudge_y = 3, angle = 0, hjust = 0.5, size = 1,vjust=0)+
  labs(x = "", y="count") +
  bw.theme
pdf("figs/keep/aligned/05_upset-TF.pdf", width = 9, height = 5)
p.upset.motif
dev.off()
################################################################################
################################################################################
################################################################################
## ATAC peak for selected gene
# plot gene activity
int.samples.so.filt <- subset(int.samples.so, cells = int.metadat$cell[int.metadat$RNA_predicted_cluster %in% c("Astro", "Micro", "Oligo", "OPC", "Exc", "Inh")])
int.samples.so.filt$cluster_stim <- factor(paste0(int.samples.so.filt$RNA_predicted_cluster, "__", ifelse(int.samples.so.filt$stimulation == "stim", "stimulated", "baseline")))
library(Seurat);library(Signac)
Idents(int.samples.so.filt) <- "cluster_stim"
bs.cell.count <- read_csv("data/derivatives/manual-boostrapped-pseudobulk/counts.csv")
source("src/utils/CoveragePlot_ME.R")
p.ccl4.peaks <- CoveragePlot_ME(obj = CoveragePlot(int.samples.so.filt, region = "CCL4", 
                                                   assay = "MACS3_peaks",
                                                   extend.upstream = 2000, extend.downstream = 1000))
p.npas4.peaks <- CoveragePlot_ME(obj = CoveragePlot(int.samples.so.filt, region = "NPAS4", 
                                                    assay = "MACS3_peaks",
                                                    extend.upstream = 2000, extend.downstream = 1000))
pdf("figs/keep/aligned/04_atac-peaks.pdf", width = 10, height = 10)
wrap_plots(p.ccl4.peaks,p.npas4.peaks, nrow = 1)
dev.off()
################################################################################
################################################################################
################################################################################
## gene expression for selected ones
genes <- c("CCL4", "NPAS4")
rna.pb.bs <- read_rds("data/derivatives/manual-boostrapped-pseudobulk/RNA/all.rds")
rna.pb.bs.meta <- data.frame(cc = colnames(rna.pb.bs)[-1]) %>%
  mutate(RNA_predicted_cluster = sub("__.*", "", cc),sampleID = sub(".*__", "", cc)) %>%
  left_join(int.metadat %>% distinct(sampleID, sampleID_2, stimulation)) %>% filter(sampleID_2 %in% paste0("P",c(4,5,6,18,20,22)))
rna.pb.bs.00 <- rna.pb.bs %>% 
  pivot_longer(cols = -c("gene"), names_to = "cc") %>% left_join(rna.pb.bs.meta) %>%
  mutate(sampleID_2 = factor(sampleID_2, levels = paste0("P", c(4:6,22,18:20))))
## radar plot
library(ggradar)
tt <- rna.pb.bs.00 %>% dplyr::filter(gene %in% genes) %>%
  mutate(stimulation = case_when(stimulation == "stim" ~ "stimulation", stimulation == "nostim" ~ "baseline")) %>%
  left_join(samples.meta) %>% 
  pivot_wider(names_from = sampleID_4, values_from = value, id_cols = c(stimulation,RNA_predicted_cluster,gene)) %>%
  select(stimulation,RNA_predicted_cluster,gene,paste0("E",c(1:4)),"G1","G2")

source("src/utils/gene_radar_plot.R")
p.a <- gene_radar_plot(data = tt, cluster = "Astro", gene.n = "CCL4",title = T)
p.ex <- gene_radar_plot(data = tt, cluster = "Exc", gene.n = "CCL4")
p.in <- gene_radar_plot(data = tt, cluster = "Inh", gene.n = "CCL4")
p.mi <- gene_radar_plot(data = tt, cluster = "Micro", gene.n = "CCL4")
p.ol <- gene_radar_plot(data = tt, cluster = "Oligo", gene.n = "CCL4")
p.op <- gene_radar_plot(data = tt, cluster = "OPC", gene.n = "CCL4", legend = T)
radar.ccl4 <- wrap_plots(p.a,p.ex,p.in,p.mi,p.ol,p.op,ncol = 1)
p.a <- gene_radar_plot(data = tt, cluster = "Astro", gene.n = "NPAS4",title = T)
p.ex <- gene_radar_plot(data = tt, cluster = "Exc", gene.n = "NPAS4")
p.in <- gene_radar_plot(data = tt, cluster = "Inh", gene.n = "NPAS4")
p.mi <- gene_radar_plot(data = tt, cluster = "Micro", gene.n = "NPAS4")
p.ol <- gene_radar_plot(data = tt, cluster = "Oligo", gene.n = "NPAS4")
p.op <- gene_radar_plot(data = tt, cluster = "OPC", gene.n = "NPAS4", legend = T)
radar.npas4 <- wrap_plots(p.a,p.ex,p.in,p.mi,p.ol,p.op,ncol = 1)
pdf("figs/keep/aligned/07_gene-exp-radar-plot.pdf", height = 20,width = 10)
wrap_plots(radar.ccl4, radar.npas4, nrow = 1)
dev.off()
################################################################################
################################################################################
################################################################################
################################################################################
## heatmap for genes from bulk and their expression in sn-multiome
p.heat <- rna.pb.bs.00 %>% left_join(samples.meta) %>%
  group_by(RNA_predicted_cluster, gene, condition, stimulation) %>% summarise(value = log2(mean(value))) %>% ungroup() %>%
  dplyr::filter(gene %in% bulk.human$gene[bulk.human$sig]) %>% mutate(stimulation = case_when(stimulation=="stim" ~ "stimulated",stimulation=="nostim" ~ "baseline")) %>%
  ggplot(aes(x = stimulation, y = gene, fill = value))+
  geom_tile() + scale_fill_gradient(low = "white", high = abstract.colors.2[3], name = expression("log"[2]~"(average expression)")) + my.guides +
  ggh4x::facet_nested(~RNA_predicted_cluster+condition, scales = "free")+
  # ggh4x::facet_grid2(rows = vars(RNA_predicted_cluster), cols = vars(sampleID_2), scales = "free") +
  bw.theme + theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "", y="")
pdf("figs/keep/bulk-DEGs-exp-heatmap-in-sn-multiomics-human.pdf", height = 16, width = 8)
p.heat;dev.off()

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
