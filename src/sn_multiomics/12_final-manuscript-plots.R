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
ns.de <- read_csv("data/derivatives/gex-DE/final/all-results_082025.csv") %>% filter(dx=="combined", -log10(lmmSeq_pval) < 100) %>%
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
                     strip = ggh4x::strip_themed(text_y = list(element_text(color = cell.colors.2[["Exc"]],size=12,face = "bold"),
                                                               element_text(color = cell.colors.2[["Inh"]],size=12,face = "bold")))) + 
  scale_color_manual(values = c(redblu.col.2[c(2,1)]))+
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

  scale_color_manual(values = c(redblu.col.2[c(2,1)]))+
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

pdf("figs/keep/aligned/02_lmmseq-GEX-volcano_082025.pdf", width = 5, height = 12)
p.volc.lmmseq
dev.off()

###### lmmSeq vs edgeR per dx
ns.de.2 <- read_csv("data/derivatives/gex-DE/final/all-results_082025.csv")


pdf("figs/keep/aligned/02_lmmseq-vs-edgeR-per-dx-volcano_082025.pdf", width = 8, height = 14)
ns.de.2 %>% mutate(cat = case_when(sig_both_cons ~ "DEG in both",
                                   lmmSeq_DEG_cons ~ "lmmSeq DEG only",
                                   edgeR_DEG_cons ~ "edgeR DEG only",
                                   T~"ns")) %>%
  ggplot(aes(lmmSeq_coef_stimulation, edgeR_logFC, color = cat)) +
  geom_vline(xintercept = c(1,-1), linetype=2, color = "pink") +geom_vline(xintercept = c(0), linetype=1, color = "black",alpha=0.2) +
  geom_hline(yintercept = c(1,-1), linetype=2, color = "pink") +geom_hline(yintercept = c(0), linetype=1, color = "black",alpha=0.2) +
  geom_point(aes(alpha = lmmSeq_DEG_cons|edgeR_DEG_cons, size = lmmSeq_DEG_cons|edgeR_DEG_cons)) +
  ggrepel::geom_text_repel(aes(label = ifelse(lmmSeq_DEG_cons|edgeR_DEG_cons,gene,"")),
                           show.legend = F,size=2, max.time = 0.1, max.overlaps = 30) +
  ggh4x::facet_grid2(rows = vars(RNA_predicted_cluster), cols = vars(dx),
                     scales = "free", space = "free") +
  scale_size_manual(values = c(0.1,1)) + scale_alpha_manual(values = c(0.1,1)) +
  scale_color_manual(values = c(palette.1[1:3],"grey"), name="") +
  labs(x = "lmmSeq coefficient for stimulation", y = expression("edgeR log"[2]~"(fold change)"),
       # caption = paste0("n(pairs): 6 (epilepsy:4, glioma:2)\n",
       #                  "edgeR model:\n",
       #                  "      combined: ~ dx + anesthesia + sexM + age + stimulation\n",
       #                  "      epilepsy/glioma: ~ stimulation\n",
       #                  "lmmSeq model:\n",
       #                  "      combined: ~ stimulation + anesthesia + sexM + age + dx + (1|participant)\n",
       #                  "      epilepsy/glioma: ~ stimulation + (1|participant)\n",
       #                  "DEGs criteria:      (abs(logFC) > 1 | abs(coefficient) > 1) & FDR < 0.05")
       ) +
  bw.theme + guides(size="none",alpha="none")+
  theme(strip.text = element_text(size=14),
        axis.text=element_text(size=10),axis.title = element_text(size=17))
dev.off()
ggsave2("figs/keep/aligned/02_lmmseq-vs-edgeR-per-dx-volcano_082025.png", width = 8, height = 14)
################################################################################
################################################################################
################################################################################
## GO results
go.res <- read_csv("data/derivatives/enrichment-analysis/goseq-GO-all.csv")
library(ggbreak)
go.res.2 <- go.res %>%
  filter(FDR < 0.05, thr == "cons") %>%
  mutate(type = sub("_.*", "", type),
         type = ifelse(type=="sig", "both methods", type)) %>% filter(type == "lmmSeq") %>%
  group_by(ontology, RNA_predicted_cluster, thr) %>%
  top_n(n = 10, wt = -FDR) %>% ungroup() %>%
  mutate(hitsPerc=numDEInCat*100/numInCat)
p.go <- go.res.2 %>%
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

## V2
p.go2 <- go.res.2 %>%
  ggplot(aes(RNA_predicted_cluster, term, color = -log10(FDR))) +
  geom_point(aes(size = numDEInCat), alpha=0.6) +
  # scale_color_gradient2(low = redblu.col.2[2], high = redblu.col.2[1]) +
  # scale_color_steps2(low = redblu.col.2[2], high = redblu.col.2[1], n.breaks = 6)+
  scale_color_gradientn(colors = colorRampPalette(redblu.col.2[c(2,1)])(10),
                        values = scales::rescale(quantile(-log10(go.res.2$FDR), 
                                                          probs = seq(0, 1, length.out = 10))))+
  my.guides +
  ggh4x::facet_grid2(rows = vars(ontology), scales = "free", space = "free") +
  bw.theme + labs(x="", y="") +
  theme(legend.box = "vertical", 
        axis.text.x.bottom = element_text(angle=90,hjust=1,vjust = 0.5))
pdf("figs/keep/aligned/09_GO-scatter_V2.pdf",height = 10,width = 7.3)
p.go2
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
bulk.mouse <- readxl::read_xlsx("data/YV/Supplementary data 4_Bulk_Stimulated_mouse_10.01.2024.xlsx", range = "A1:G15399") %>%
  rename("gene" = 2) %>% mutate(gene = toupper(gene)) %>% rename_at(.vars = -c(1:2), .funs = function(x) paste0(x, "_mouse")) %>%
  mutate(sig = FDR_mouse<0.05&abs(logFC_mouse)>0.2)
bulk.human <- readxl::read_xlsx("data/YV/Supplementary data 2_Bulk_Stimulated_human_10.01.2024.xlsx") %>%
  rename("gene" = 2)  %>% rename_at(.vars = -c(1:2), .funs = function(x) paste0(x, "_human")) %>%
  mutate(sig = FDR_human<0.05&abs(logFC_human)>0.2)

inter.df <- data.frame(gene = unique(c(ns.de$gene, bulk.human$gene, bulk.mouse$gene))) %>%
  mutate(Micro_human = as.numeric(gene %in% ns.de$gene[(ns.de$lmmSeq_DEG_cons)&ns.de$RNA_predicted_cluster=="Micro"]),
         Astro_human = as.numeric(gene %in% ns.de$gene[(ns.de$lmmSeq_DEG_cons)&ns.de$RNA_predicted_cluster=="Astro"]),
         Oligo_human = as.numeric(gene %in% ns.de$gene[(ns.de$lmmSeq_DEG_cons)&ns.de$RNA_predicted_cluster=="Oligo"]),
         OPC_human = as.numeric(gene %in% ns.de$gene[(ns.de$lmmSeq_DEG_cons)&ns.de$RNA_predicted_cluster=="OPC"]),
         Exc_human = as.numeric(gene %in% ns.de$gene[(ns.de$lmmSeq_DEG_cons)&ns.de$RNA_predicted_cluster=="Exc"]),
         Inh_human = as.numeric(gene %in% ns.de$gene[(ns.de$lmmSeq_DEG_cons)&ns.de$RNA_predicted_cluster=="Inh"]),
         Bulk_human = as.numeric(gene %in% bulk.human$gene[bulk.human$sig]),
         Bulk_mouse = as.numeric(gene %in% bulk.mouse$gene[bulk.mouse$sig]))
inter.df.long <- inter.df %>%
  pivot_longer(cols = -1) %>% filter(value == 1) %>%
  pivot_wider(names_from = value, values_from = name) %>% rename("cl_list" = 2)
## full supplementary version
p.upset <- inter.df.long %>%
  ggplot(aes(x = cl_list)) +
  geom_bar(width = 0.7) +
  scale_x_upset() + 
  geom_text(stat='count', aes(label=after_stat(count)), vjust=0, size = 2) +
  geom_text(data = inter.df.long %>% group_by(cl_list) %>% summarize(gene_names = paste(gene, collapse = ", "), gene_count = n()) %>% 
              filter(gene_count<20),
            aes(x = cl_list, y = gene_count +4,label = gene_names), 
            nudge_y = 0, hjust = 0,vjust=0, size = 1.8, angle=90)+
  labs(x = "", y="count",
       caption = paste0("DEGs criteria:\n",
                        "    Bulk_human:    abs(logFC)>0.2 & FDR < 0.05\n",
                        "    Bulk_mouse:    abs(logFC)>0.2 & FDR < 0.05\n",
                        "    any sn-RNA-Seq_human:    abs(stimulation coefficient) > 1 & FDR < 0.05")) +
  bw.theme
pdf("figs/keep/aligned/05_upset-supp.pdf", width = 9, height = 6.5)
p.upset
dev.off()

## main panel version
inter.df.long2 <- inter.df %>% select(-contains("Bulk")) %>%
  pivot_longer(cols = -1) %>% filter(value == 1) %>%
  pivot_wider(names_from = value, values_from = name) %>% rename("cl_list" = 2)


p.upset2 <- inter.df.long2 %>%
  ggplot(aes(x = cl_list)) +
  geom_bar(width = 0.7) +
  scale_x_upset() + 
  scale_color_manual(values = cell.colors.2) +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=0, size = 2) +
  geom_text(data = inter.df.long2 %>% group_by(cl_list) %>% summarize(gene_names =str_replace_all(paste0(paste(gene[1:min(length(gene),7)], collapse = ", "),
                                                                                                         ifelse(length(gene)>7,paste0("\n",paste(gene[8:length(gene)], collapse = ", ")),"")),
                                                                                                  ",NA",""), gene_count = n()) %>% filter(gene_count<20),
            aes(x = cl_list, y = gene_count +4,label = gene_names), 
            nudge_y = 0, angle = 90, hjust = 0, size = 2.3,vjust=0.5)+
  labs(x = "", y="count") +
  bw.theme
p.upset2

micro.t <- as.character(inter.df.long2 %>% filter(cl_list %in% list("Micro_human")) %>% select(gene) %>% unlist())
# Create the table for micro genes
library(gridExtra);library(grid)
table_grob <- tableGrob(cbind(micro.t[1:38],micro.t[39:76],c(micro.t[77:112],"","")),
                        theme = ttheme_minimal(
                          base_size = 6,base_colour = "black",
                          base_family = "sans",parse = FALSE,
                          padding = unit(c(1, 1), "mm"),
                          core = list(fg_params = list(hjust = 0, x = 0.1)),
                          colhead = list(fg_params = list(fontface = "bold"))))
table_grob <- gtable::gtable_add_grob(table_grob,
  grobs = rectGrob(gp = gpar(fill = NA, lwd = 1.5)), 
  t = 1, b = nrow(table_grob), l = 1, r = ncol(table_grob))
table_grob <- gtable::gtable_add_rows(table_grob,heights = unit(0.2, "in"), pos = 0)
table_grob <- gtable::gtable_add_grob(table_grob,textGrob("Microglia-specific genes", gp = gpar(fontsize = 8)),
                                      t = 1, l = 1, r = ncol(table_grob))

wrap_plots(table_grob,plot_spacer(),ncol=1,heights = c(1,0.2))
wrap_plots(plot_spacer(),
           # table_grob,
           wrap_plots(table_grob,plot_spacer(),ncol=1,heights = c(1,0.2)),
           plot_spacer(),p.upset2,nrow = 1,
           widths = c(0.5,1,0.5,6))
wrap_elements(full = table_grob) + p.upset2 +
  plot_layout(widths = c(1, 6))
p.upset2.full
pdf("figs/keep/aligned/05_upset.pdf", width = 8, height = 5.5)
p.upset2
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
          motifs = unique(motif.all$motif[which(motif.all$motif.name %in% 
                                                  # c("HINFP", "NRF1", "ZBTB33", "ZBTB14", "TCFL5", "KLF15"))])) +
                                                  c("EGR1", "EGR2", "EGR3", "EGR4"))])) +
  scale_fill_manual(values = palette.1) +
  facet_wrap(~seq_group, scales = "free_x", nrow = 1) +
  bw.theme
pdf("figs/keep/aligned/06_TF-motifs_V3.pdf", width = 10, height = 3)
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
              summarize(motif_names = paste(motif.name, collapse = ", "), 
                        motif_count = n()) %>% filter(motif_count<20),
            aes(x = cl_list, y = motif_count +3,label = motif_names), nudge_y = 0, 
            angle = 90, hjust = 0, size = 1.7,vjust=0.5)+
  labs(x = "", y="count") +
  bw.theme
pdf("figs/keep/aligned/05_upset-TF_V2.pdf", width = 9, height = 5.5)
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
  geom_tile() + scale_fill_gradient(low = "white", 
                                    # high = abstract.colors.2[3], 
                                    high = palette.1[1],
                                    name = expression("log"[2]~"(average expression)")) + my.guides +
  ggh4x::facet_nested(~RNA_predicted_cluster+condition, scales = "free", 
                      nest_line = element_line(linewidth=0.7), solo_line = T)+
  # ggh4x::facet_grid2(rows = vars(RNA_predicted_cluster), cols = vars(sampleID_2), scales = "free") +
  bw.theme + theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "", y="")+theme(axis.text = element_text(size=10),
                           strip.text = element_text(size=14))
p.heat
pdf("figs/keep/bulk-DEGs-exp-heatmap-in-sn-multiomics-human.pdf", height = 16, width = 11)
p.heat;dev.off()

################################################################################
################################################################################
################################################################################
### additional figures requested by Ted and Jyoti

### ATAC UMAP
library(Seurat);library(Signac)
int.samples.so.v2 <- int.samples.so
DefaultAssay(int.samples.so.v2) <- "Gene.Activity"
int.samples.so.v2 <- ScaleData(int.samples.so.v2, verbose = FALSE)
int.samples.so.v2 <- FindVariableFeatures(int.samples.so.v2)
int.samples.so.v2 <- RunPCA(int.samples.so.v2, npcs = 30, verbose = FALSE)
int.samples.so.v2 <- RunUMAP(int.samples.so.v2, reduction = "pca", dims = 1:30)
int.samples.so.v2 <- FindNeighbors(int.samples.so.v2, reduction = "pca", dims = 1:30)
int.samples.so.v2 <- FindClusters(int.samples.so.v2, resolution = 0.5)


int.metadat.v2 <- int.metadat %>%
  left_join(int.samples.so.v2@reductions$umap@cell.embeddings %>%
              as.data.frame() %>% rownames_to_column("cell") %>%
              rename(atac_umap_1=2, atac_umap_2=3))
p.umap.2 <- int.metadat.v2 %>%
  ggplot(aes(x=atac_umap_1, y = atac_umap_2, color = RNA_predicted_cluster)) +
  geom_point(show.legend = F,size = 0.5) + scale_color_manual(values = cell.colors.2) +
  annotate(geom = "text", x=7,y=5, label = "Exc", color = cell.colors.2[["Exc"]], fontface="bold")+
  annotate(geom = "text", x=8,y=-2, label = "Inh", color = cell.colors.2[["Inh"]], fontface="bold")+
  annotate(geom = "text", x=7,y=10, label = "Micro", color = cell.colors.2[["Micro"]], fontface="bold")+
  annotate(geom = "text", x=-7,y=0, label = "Oligo", color = cell.colors.2[["Oligo"]], fontface="bold") +
  annotate(geom = "text", x=3,y=-7, label = "OPC", color = cell.colors.2[["OPC"]], fontface="bold")+
  annotate(geom = "text", x=0,y=6, label = "Astro", color = cell.colors.2[["Astro"]], fontface="bold")+
  bw.theme +
  labs(x = "UMAP-1", y="UMAP-2")
pdf("figs/keep/aligned/01_4A-atac-umap.pdf",height = 5,width = 5)
p.umap.2
dev.off()
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
### cell gene markers
# DefaultAssay(int.samples.so) <- "RNA"
# cell.markers <- data.frame(celltype = c(rep("Exc",7), rep("Inh",10),rep("Astro",3),
#                                         rep("Micro",6),rep("OPC",3),rep("Oligo",2)),
#                            gene = c("SLC17A7","SATB2","RORB","CUX2","TLE4","NR4A2","SEMA3C",
#                                     "GAD1","GAD2","SOX6","PVALB","SST","VIP","LHX6","NDNF","CALB2","SULF1",
#                                     "AQP4","GJB6","FGFR3",
#                                     "CTSS","C1QB","CD11b","CD45","CX3CR1","TMEM119",
#                                     "CSPG4","PDGFRA","VCAN",
#                                     "MOG","MAG"))
# rna.data <- int.samples.so@assays$RNA$data[rownames(int.samples.so)%in%cell.markers$gene,] %>%as.matrix()%>%
#   as.data.frame() %>% rownames_to_column("gene") %>%
#   pivot_longer(cols = -c(gene), names_to = "cell", values_to = "RNA") %>%
#   inner_join(int.metadat %>% select(cell, RNA_predicted_cluster, sampleID)) %>%
#   inner_join(int.samples.so@assays$Gene.Activity$data[rownames(int.samples.so@assays$Gene.Activity$data)%in%cell.markers$gene,] %>%
#                as.matrix()%>%as.data.frame() %>% rownames_to_column("gene") %>%
#                pivot_longer(cols = -c(gene), names_to = "cell", values_to = "ATAC"))
# rna.data %>% left_join(cell.markers) %>% 
#   # filter(celltype == "Oligo") %>%
#   filter(gene %in% c("AQP4","SLC17A7", "GAD1", "C1QB","VCAN","MOG")) %>%
#   mutate(gene = factor(gene, levels = c("AQP4","SLC17A7","GAD1","C1QB","MOG" ,"VCAN"))) %>%
#   pivot_longer(cols = c(RNA, ATAC), names_to = "assay", values_to = "count") %>%
#   ggplot(aes(x = RNA_predicted_cluster, y = count, fill = RNA_predicted_cluster)) +
#   geom_violin(show.legend = F) + scale_fill_manual(values = cell.colors.2) +
#   ggh4x::facet_grid2(rows = vars(gene), cols = vars(assay),
#                      scales = "free") +
#   bw.theme + labs(x = "cell cluster")

pb <- read_rds("data/derivatives/manual-boostrapped-pseudobulk/RNA/all.rds")
pb.atac <- read_rds("data/derivatives/manual-boostrapped-pseudobulk/MACS-GA/all.rds")
pb %>%
  pivot_longer(cols = -c(gene)) %>%
  mutate(celltype = sub("__.*", "", name)) %>%
  filter(gene %in% c("AQP4","SLC17A7", "GAD1", "C1QB","VCAN","MOG")) %>%
  mutate(gene = factor(gene, levels = c("AQP4","SLC17A7","GAD1","C1QB","MOG" ,"VCAN"))) %>%
  group_by(gene, celltype) %>% summarise(value = mean(value)) %>%
  ungroup() %>% group_by(gene) %>% mutate(value = scale(value,T,T)[,1]) %>% ungroup() %>% 
  ggplot(aes(x = celltype, y = gene, fill = value, size = value)) +
  geom_dotplot() + redblu.col.gradient.2(label = "scaled expression") +
  bw.theme
## lollipop
pb %>%
  pivot_longer(cols = -c(gene)) %>%
  mutate(celltype = sub("__.*", "", name)) %>%
  filter(gene %in% c("AQP4", "SLC17A7", "GAD1", "C1QB", "VCAN", "MOG")) %>%
  mutate(gene = factor(gene, levels = c("VCAN","MOG","C1QB","AQP4","GAD1","SLC17A7"))) %>%
  group_by(gene, celltype) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  group_by(gene) %>%
  mutate(value = scale(value, center = TRUE, scale = TRUE)[,1],
         y_pos = as.numeric(gene) + (as.numeric(factor(celltype)) - 3.5) * 0.1) %>%  # spread out by celltype
  ungroup() -> plot_data
# Get top celltype per gene
label_data <- plot_data %>%
  group_by(gene) %>% slice_max(order_by = value, n = 1) %>% ungroup()
pb.atac %>%
  pivot_longer(cols = -c(gene)) %>%
  mutate(celltype = sub("__.*", "", name)) %>%
  # filter(gene %in% cell.markers$gene) %>% mutate(gene = factor(gene,levels=cell.markers$gene)) %>%
  filter(gene %in% c("GJB6", "SLC17A7", "GAD2", "C1QB", "PDGFRA", "MOG")) %>%
  mutate(gene = factor(gene, levels = c("SLC17A7", "GAD2", "GJB6", "C1QB", "MOG","PDGFRA"))) %>%
  group_by(gene, celltype) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  group_by(gene) %>%
  mutate(value = scale(value, center = TRUE, scale = TRUE)[,1],
         y_pos = as.numeric(gene) + (as.numeric(factor(celltype)) - 3.5) * 0.1) %>%  # spread out by celltype
  ungroup() -> plot_data.atac
label_data.atac <- plot_data.atac %>%
  group_by(gene) %>% slice_max(order_by = value, n = 1) %>% ungroup()

loli.rna <- plot_data %>%
  ggplot(aes(x = value, y = y_pos, color = celltype)) +
  geom_segment(aes(x = 0, xend = value, y = y_pos, yend = y_pos), 
               color = "grey70", size = 0.4) +
  geom_point(size = 3, show.legend = F) + geom_vline(xintercept = 0, linetype=1,color="grey")+
  geom_text(data = label_data, aes(label = celltype), 
            hjust = 1, color = "black", size = 3.5, vjust = -1) +
  scale_y_continuous(breaks = 1:6, labels = levels(plot_data$gene), expand = expansion(add = 0.3)) +
  scale_color_manual(values = cell.colors.2)+
  labs(y = "Gene", x = "Scaled Average RNA Expression") +
  bw.theme +
  theme(panel.grid.major.y = element_blank())
loli.atac <- plot_data.atac %>%
  ggplot(aes(x = value, y = y_pos, color = celltype)) +
  geom_segment(aes(x = 0, xend = value, y = y_pos, yend = y_pos), 
               color = "grey70", size = 0.4) +
  geom_point(size = 3, show.legend = F) + geom_vline(xintercept = 0, linetype=1,color="grey")+
  geom_text(data = label_data.atac, aes(label = celltype), 
            hjust = 1, color = "black", size = 3.5, vjust = -0.5) +
  scale_y_continuous(breaks = 1:6, labels = levels(plot_data.atac$gene), 
                     expand = expansion(add = 0.3)) +
  scale_color_manual(values = cell.colors.2)+
  labs(y = "Gene", x = "Scaled Average ATAC Gene Activity") +
  bw.theme +
  theme(panel.grid.major.y = element_blank())
pdf("figs/keep/aligned/10_4B-marker-genes-per-celltype_V3.pdf", height = 7, width = 7)
patchwork::wrap_plots(loli.rna,loli.atac,nrow = 1)
dev.off()
################################################################################
################################################################################
################################################################################
### comprehensive TF motif

mm.data <- motif.all %>%
  filter(source=="global") %>% mutate(sigg = p.adjust<0.05) %>%
  left_join(ns.de %>% dplyr::filter(dx=="combined") %>% 
              select(motif.name = gene, RNA_predicted_cluster, lmmSeq_DEG_cons)) %>%
  mutate(lmmSeq_DEG_cons = ifelse(is.na(lmmSeq_DEG_cons),F,lmmSeq_DEG_cons),
         fontface_label = ifelse(sigg & lmmSeq_DEG_cons, "bold", "plain"),
         font_size = ifelse(sigg & lmmSeq_DEG_cons, 2.5, 2))
p.volc.motif <- mm.data %>%
  ggplot(aes(x = log2(fold.enrichment), y = -log10(pvalue),
             label= ifelse(p.adjust<0.05, motif.name, ""))) +
  geom_point(aes(color = ifelse(sigg, RNA_predicted_cluster, "grey"),
                 alpha = ifelse(sigg, 1, 0.1),
                 size = sigg), 
             show.legend = F) + 
  scale_color_manual(values = c(cell.colors.2, "grey")) +
  ggnewscale::new_scale_color() + 
  geom_vline(xintercept = c(0), linetype=1, color = "black",alpha=0.2) +
  geom_hline(yintercept = -log10(0.01), linetype=2, color = "pink") +
  ggrepel::geom_text_repel(aes(label = ifelse(sigg&lmmSeq_DEG_cons,motif.name,"")),
                           fontface = mm.data$fontface_label,
                           size=mm.data$font_size, max.time = 0.1, max.overlaps = 500) +
  # scale_color_manual(values = c("black", "#00BFC4"),
  #                    name = "DEG in RNA differential expression") +
  ggnewscale::new_scale_color() + 
  ggh4x::facet_grid2(rows = vars(RNA_predicted_cluster), scales = "free",
                     strip = ggh4x::strip_themed(text_y = list(element_text(color = cell.colors.2[["Astro"]],size=12,face = "bold"),
                                                               element_text(color = cell.colors.2[["Exc"]],size=12,face = "bold"),
                                                               element_text(color = cell.colors.2[["Inh"]],size=12,face = "bold"),
                                                               element_text(color = cell.colors.2[["Micro"]],size=12,face = "bold"),
                                                               element_text(color = cell.colors.2[["Oligo"]],size=12,face = "bold"),
                                                               element_text(color = cell.colors.2[["OPC"]],size=12,face = "bold")))) + 
  scale_color_manual(values = c(redblu.col[c(2,1)]))+
  scale_size_manual(values = c(0.3,1.2)) +
  labs(x = expression("-log"[2]~"fold nrichment"), y = expression("-log"[10]~"(p-value)")) +
  bw.theme
pdf("figs/keep/aligned/11_4D-motif-volcano.pdf", width = 5, height = 12)
p.volc.motif
dev.off()

library(ggseqlogo)
p.motif <- Signac::MotifPlot(atac.pb.so, 
                             motifs = c("MA0162.4","MA0472.2","MA0732.1","MA0733.1")) +
  scale_fill_manual(values = palette.1) +
  facet_wrap(~seq_group, scales = "free_x") +
  bw.theme
pdf("figs/keep/aligned/12_TF-motifs_V2.pdf", width = 5, height = 5)
p.motif;dev.off()


motif.all %>%
  filter(source=="global",
         motif.name %in% paste0("EGR",c(1:4))) %>%
  ggplot(aes(x = RNA_predicted_cluster, y = motif.name, size = fold.enrichment,
             color = fold.enrichment)) +
  geom_point() + scale_color_gradient2(low = redblu.col.2[2],high = redblu.col.2[1],name = "fold enrichment") +
  bw.theme
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
## wrap plots
# run the neuroestimator plot and get it here
ne.plot <- read_rds("figs/keep/aligned/03_neuroestimator_0919.rds")
pdf("figs/keep/aligned/fig-03_V2.pdf", height = 14, width = 12)
wrap_plots(wrap_plots(wrap_plots(plot_spacer(), p.umap, loli.rna,widths = c(1,1.4,0.7)), 
                      p.upset2, ne.plot, ncol = 1,heights = c(1,1.3,1.3))
           ,p.volc.lmmseq, nrow = 1, widths = c(1.8,1))
dev.off()

pdf("figs/keep/aligned/fig-03_V3.pdf", height = 14, width = 12)
wrap_plots(wrap_plots(wrap_plots(plot_spacer(), p.umap, loli.rna,widths = c(1,1.4,0.7)), 
                      wrap_plots(wrap_plots(plot_spacer(),table_grob,plot_spacer(),ncol=1,
                                            heights = c(0.2,1,0.2)),
                                 p.upset2,nrow = 1,
                                 widths = c(1,6)), 
                      ne.plot, ncol = 1,heights = c(1,1.3,1.3))
           ,p.volc.lmmseq, nrow = 1, widths = c(1.8,1))
dev.off()

pdf("figs/keep/aligned/fig-04.pdf", height = 10, width = 9)
wrap_plots(wrap_plots(wrap_plots(p.umap.2, loli.atac, nrow = 1),
                      p.motif, ncol = 1, heights = c(1,1.5)), 
           p.volc.motif, widths = c(3,1))
dev.off()
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
