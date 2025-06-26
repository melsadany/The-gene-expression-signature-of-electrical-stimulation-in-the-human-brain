################################################################################
################################################################################
rm(list = ls()); gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
              "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
library(patchwork)
library(Seurat);library(Signac)
registerDoMC(30)
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/neurosurgery/src/98_manuscript-color-palettes.R"))
################################################################################
################################################################################
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/neurosurgery")
setwd(project.dir)
################################################################################
################################################################################
integrated.samples.so <- read_rds("data/derivatives/ME-processing/processed-seurat-objects/integrated-multiome-samples.rds")
int.metadata.filtered2 <- read_rds("data/derivatives/cells-metadata-filtered.rds")
celltypes.meta <- read_rds(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/data/Allen/human-m1-10x/celltypes-metadata.rds"))
samples.meta <- read_csv("data/samples-metadata_ME.csv") 
bs.cell.count <- read_csv("data/derivatives/manual-boostrapped-pseudobulk/counts.csv")
cl.keep <- c("OPC","Oligo","Exc","Micro","Inh","Astro")
################################################################################
################################################################################
################################################################################
################################################################################
## RNA
rna.counts <- GetAssayData(integrated.samples.so, 
                           assay = "RNA", layer = "counts") %>%
  as.data.frame()
rownames(rna.counts) <- rownames(integrated.samples.so@assays$RNA@features)
colnames(rna.counts) <- rownames(integrated.samples.so@assays$RNA@cells)
gc()
## ATAC
atac.counts <- GetAssayData(integrated.samples.so, 
                            assay = "Gene.Activity", layer = "counts") %>%
  as.data.frame()
rownames(atac.counts) <- integrated.samples.so@assays$Gene.Activity@counts@Dimnames[[1]]
colnames(atac.counts) <- integrated.samples.so@assays$Gene.Activity@counts@Dimnames[[2]]
gc()
## ATAC frags
atac.frags.r <- read_rds("data/derivatives/calc-ATAC-frags/genes/all-PRM-TES.rds")
atac.frags <- atac.frags.r %>% 
  column_to_rownames("cell") %>% t() %>% 
  as.data.frame()

gc()
################################################################################
################################################################################
library(neuroestimator)

rna.act <- neuroestimator(rna.counts, species = "hsapiens") 
atac.act <- neuroestimator(atac.counts, species = "hsapiens") 
atac.frag.act <- neuroestimator(atac.frags, species = "hsapiens") 

ne.act <- inner_join(rna.act %>%
                       rename(RNA_predicted_activity = 1) %>%
                       rownames_to_column("cell"),
                     atac.act %>%
                       rename(ATAC_predicted_activity = 1) %>%
                       rownames_to_column("cell")) %>%
  inner_join(atac.frag.act %>%
               rename(ATAC_frags_predicted_activity = 1) %>%
               rownames_to_column("cell")) %>%
  inner_join(int.metadata.filtered2)
write_rds(ne.act, "data/derivatives/neuroestimator-results.rds", compress = "gz")
ne.act <- read_rds("data/derivatives/neuroestimator-results.rds")
################################################################################
################################################################################
################################################################################
## bootstrapping by participant
registerDoMC(cores = 4)
ne.bs <- foreach(jk = 1:length(cl.keep), .combine = rbind) %dopar% {
  c <- cl.keep[jk]
  c.ne <- ne.act %>% 
    dplyr::filter(RNA_predicted_cluster == c, cell %in% int.metadata.filtered2$cell) %>%
    dplyr::select(sampleID_3, sampleID_4, stimulation,
                  RNA_predicted_activity, ATAC_predicted_activity, ATAC_frags_predicted_activity)
  set.seed(123)
  s.count <- bs.cell.count$bs_count[bs.cell.count$RNA_predicted_cluster==c]
  
  bs.res <- foreach(i = 1:1000, .combine = rbind) %dopar% {
    bs.data <- c.ne %>%
      group_by(sampleID_3) %>%
      sample_n(s.count, replace = TRUE) %>%
      ungroup()
    
    ttdf <- bs.data %>%
      group_by(sampleID_3) %>%
      dplyr::summarise(RNA_predicted_activity = mean(RNA_predicted_activity),
                       ATAC_predicted_activity = mean(ATAC_predicted_activity),
                       ATAC_frags_predicted_activity = mean(ATAC_frags_predicted_activity)) %>%
      left_join(int.metadata.filtered2 %>% ungroup() %>% 
                  distinct(sampleID_3, sampleID_4, stimulation))
    
    ## get activity difference for each participant
    activity_diff <- inner_join(ttdf %>%
                                  pivot_wider(names_from = stimulation, values_from = RNA_predicted_activity, id_cols = sampleID_4) %>%
                                  mutate(RNA_activity_diff = stim - nostim) %>% 
                                  dplyr::select(sampleID_4, RNA_activity_diff),
                                ttdf %>%
                                  pivot_wider(names_from = stimulation, values_from = ATAC_predicted_activity, id_cols = sampleID_4) %>%
                                  mutate(ATAC_activity_diff = stim - nostim) %>% 
                                  dplyr::select(sampleID_4, ATAC_activity_diff)) %>%
      inner_join(ttdf %>%
                   pivot_wider(names_from = stimulation, values_from = ATAC_frags_predicted_activity, id_cols = sampleID_4) %>%
                   mutate(ATAC_frags_activity_diff = stim - nostim) %>% 
                   dplyr::select(sampleID_4, ATAC_frags_activity_diff))
    
    ### get the null distribution
    null.bs.data <- bs.data %>%
      dplyr::select(4:6) %>%
      mutate(stimulation = sample(c(rep("stimulated", nrow(bs.data)/2), rep("unstimulated", nrow(bs.data)/2))))
    null.rna <- mean(null.bs.data$RNA_predicted_activity[which(null.bs.data$stimulation=="stimulated")])-
      mean(null.bs.data$RNA_predicted_activity[which(null.bs.data$stimulation=="unstimulated")])
    null.atac <- mean(null.bs.data$ATAC_predicted_activity[which(null.bs.data$stimulation=="stimulated")])-
      mean(null.bs.data$ATAC_predicted_activity[which(null.bs.data$stimulation=="unstimulated")])
    null.atac.frags <- mean(null.bs.data$ATAC_frags_predicted_activity[which(null.bs.data$stimulation=="stimulated")])-
      mean(null.bs.data$ATAC_frags_predicted_activity[which(null.bs.data$stimulation=="unstimulated")])
    
    ## combine and return
    wide.df <- activity_diff %>%
      pivot_longer(cols = -sampleID_4) %>%
      mutate(name2 = paste0(sub("_activity_diff", "", name), "__", sampleID_4)) %>%
      select(name2,value) %>%
      column_to_rownames("name2") %>% t()
    ## return
    data.frame(iteration = i,
               wide.df,
               RNA__null = null.rna,
               ATAC__null = null.atac,
               ATAC_frags__null = null.atac.frags)
  }
  
  ## combine from all iterations, get mean, and CI
  df33 <- do.call(rbind, lapply(bs.res[,-1], FUN = function(x) mean(x))) %>%
    as.data.frame() %>%
    rownames_to_column("source_1") %>%
    rename("mean_difference" = 2) %>%
    mutate(sampleID_4 = as.character(sub(".*__", "", source_1)),
           type = sub("__.*", "", source_1),
           celltype = c,
           confin_min = do.call(rbind, lapply(bs.res[,-1], FUN = function(x) quantile(x, 0.025))) %>% as.numeric(),
           confin_max = do.call(rbind, lapply(bs.res[,-1], FUN = function(x) quantile(x, 0.975))) %>% as.numeric(),
           confin_min_90 = do.call(rbind, lapply(bs.res[,-1], FUN = function(x) quantile(x, 0.05))) %>% as.numeric(),
           confin_max_90 = do.call(rbind, lapply(bs.res[,-1], FUN = function(x) quantile(x, 0.95))) %>% as.numeric(),
           confin_min_85 = do.call(rbind, lapply(bs.res[,-1], FUN = function(x) quantile(x, 0.075))) %>% as.numeric(),
           confin_max_85 = do.call(rbind, lapply(bs.res[,-1], FUN = function(x) quantile(x, 0.925))) %>% as.numeric(),
           sampleID_4 = ifelse(grepl("null", source_1), "all", sampleID_4),
           type = ifelse(grepl("null", source_1), sub("null_", "", type), type))
  
  return(df33)
}
write_rds(ne.bs, "data/derivatives/neuroestimator-bootstrapping-results.rds")
ne.bs <- read_rds("data/derivatives/neuroestimator-bootstrapping-results.rds")
################################################################################
################################################################################
################################################################################

## get participants data for plotting
ne.bs2 <- ne.bs %>%
  dplyr::filter(sampleID_4 != "all") %>%
  pivot_wider(names_from = type, values_from = mean_difference, id_cols = c(celltype, sampleID_4),
              names_prefix = "mean__") %>%
  inner_join(inner_join(ne.bs %>% dplyr::filter(sampleID_4 != "all") %>%
                          pivot_wider(names_from = type, values_from = confin_min, id_cols = c(celltype, sampleID_4),
                                      names_prefix = "confin_min__"),
                        ne.bs %>% dplyr::filter(sampleID_4 != "all") %>%
                          pivot_wider(names_from = type, values_from = confin_max, id_cols = c(celltype, sampleID_4),
                                      names_prefix = "confin_max__"))) %>%
  mutate(RNA_sig = (sign(confin_min__RNA)==sign(confin_max__RNA)),
         ATAC_sig = (sign(confin_min__ATAC)==sign(confin_max__ATAC)),
         ATAC_frags_sig = (sign(confin_min__ATAC_frags)==sign(confin_max__ATAC_frags)),
         sig = RNA_sig & ATAC_sig,
         sig_frags = RNA_sig & ATAC_frags_sig)

## get the null distribution eclipse data
ne.bs.null <- ne.bs %>%
  filter(sampleID_4 == "all") %>%
  pivot_wider(names_from = type, values_from = mean_difference, id_cols = c(celltype, sampleID_4),
              names_prefix = "center__") %>%
  inner_join(ne.bs %>%
               dplyr::filter(sampleID_4 == "all") %>%
               mutate(diameter = confin_max -confin_min,
                      diameter_90 = confin_max_90 -confin_min_90,
                      diameter_85 = confin_max_85 -confin_min_85) %>%
               pivot_longer(cols = starts_with("diameter")) %>%
               mutate(type_2 = paste0(name, "__", type)) %>%
               pivot_wider(names_from = type_2, values_from = value, 
                           id_cols = c(celltype, sampleID_4))) %>%
  select(celltype, starts_with(c("center", "diameter")))


#### plot 
ne.bs2 %>%
  dplyr::filter(sampleID_4 != "all") %>%
  left_join(ne.bs.null) %>%
  left_join(samples.meta %>% mutate(sampleID_4 = paste0("P", parse_number(sampleID_3)))) %>%
  mutate(sampleID_4 = factor(sampleID_4, levels = names(sample.colors.2))) %>%
  mutate(sig_3 = case_when(RNA_sig&ATAC_sig ~ "RNA and ATAC",
                           RNA_sig&!ATAC_sig ~ "RNA only",
                           !RNA_sig&ATAC_sig ~ "ATAC only",
                           !RNA_sig&!ATAC_sig ~ "neither"),
         sig_3 = factor(sig_3, levels = c("RNA and ATAC", "RNA only", "ATAC only", "neither")),
         sampleID_4 = factor(sampleID_4, levels = c(paste0("E",c(1:4)),"G1","G2")),
         celltype = factor(celltype, levels = c("Exc","Inh","Astro","Micro","Oligo","OPC"))) %>%
  ggplot(aes(x = mean__ATAC, y = mean__RNA, 
             # color = sampleID_4, 
             color = celltype,
             alpha = sig_3)) +
  geom_point(aes(shape = sig_3),position = position_dodge(width = 0.6), size =2.5) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2, color = "red") +
  geom_errorbarh(aes(xmin = confin_min__ATAC, xmax = confin_max__ATAC),
                 linewidth = 0.8, height = 0, show.legend = F, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = confin_min__RNA, ymax = confin_max__RNA),
                linewidth = 0.8, width = 0, show.legend = F, position = position_dodge(width = 0.6)) +
  xlim(c(-0.075,0.2)) +
  # ylim(c(-0.15,0.15)) +
  ggforce::geom_ellipse(aes(x0 = center__ATAC, y0 = center__RNA,
                            a = diameter__ATAC, b = diameter__RNA,
                            angle = 0), 
                        linetype = 2, alpha = 0.3, size = 0.3, color = "grey") +
  # ggforce::geom_ellipse(aes(x0 = center__ATAC, y0 = center__RNA,
  #                           a = diameter_90__ATAC, b = diameter_90__RNA,
  #                           angle = 0), 
  #                       linetype = 2, alpha = 0.3, size = 0.3, color = "grey") +
  # ggforce::geom_ellipse(aes(x0 = center__ATAC, y0 = center__RNA,
  #                           a = diameter_85__ATAC, b = diameter_85__RNA,
  #                           angle = 0), 
  #                       linetype = 2, alpha = 0.3, size = 0.3, color = "grey") +
  # scale_alpha_manual(values = c(0.3,1), name = "significant in both RNA and ATAC predicted activity") +
  scale_alpha_manual(values = c("RNA and ATAC" = 1,
                                "RNA only" = 0.4,
                                "ATAC only" = 0.4,
                                "neither" = 0.1), 
                     guide = "legend",
                     name = "Activity difference significance") +
  scale_shape_manual(values = c("RNA and ATAC" = 0,
                                "RNA only" = 2,
                                "ATAC only" = 8,
                                "neither" = 13),
                     name = "Activity difference significance",
                     guide = "legend") +
  # scale_color_manual(values = sample.colors.2, name = "participant") +
  scale_color_manual(values = cell.colors.2, name = "cell type", guide = "none") +
  guides(alpha = guide_legend(override.aes = list(shape = c(0, 2, 8, 13)), nrow = 1)) +
  ggh4x::facet_wrap2(~celltype, nrow = 2, scales = "free",
                     strip = ggh4x::strip_themed(text_x = list(element_text(color = cell.colors.2[["Exc"]],size=12,face = "bold"),
                                                               element_text(color = cell.colors.2[["Inh"]],size=12,face = "bold"),
                                                               element_text(color = cell.colors.2[["Astro"]],size=12,face = "bold"),
                                                               element_text(color = cell.colors.2[["Micro"]],size=12,face = "bold"),
                                                               element_text(color = cell.colors.2[["Oligo"]],size=12,face = "bold"),
                                                               element_text(color = cell.colors.2[["OPC"]],size=12,face = "bold")))) +
  labs(x = "mean difference in NEUROeSTIMator predicted activity (ATAC)",
       y = "mean difference in NEUROeSTIMator predicted activity (RNA)",
       caption = paste0("data points have error bars for the 95% CI\n",
                        "the dashed ellipse show the null distribution with 95% CIs per cell type")) +
  bw.theme +
  theme(panel.grid.major = element_line(linewidth = 0.05),
        panel.grid.minor = element_line(linewidth = 0.05),
        legend.box = "vertical", legend.spacing.y = unit(0.0005, units = "in"))
ggsave2("figs/keep/aligned/03_neuroestimator_0626.pdf", width = 9.5, height = 7)


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
#### do it by stimulation state
registerDoMC(cores = 4)
ne.st.bs <- foreach(jk = 1:length(cl.keep), .combine = rbind) %dopar% {
  c <- cl.keep[jk]
  c.ne <- ne.act %>% 
    dplyr::filter(RNA_predicted_cluster == c, cell %in% int.metadata.filtered2$cell) %>%
    left_join(samples.meta) %>%
    dplyr::select(sampleID_3, sampleID_4, stimulation,condition,
                  RNA_predicted_activity, ATAC_predicted_activity, ATAC_frags_predicted_activity)
  set.seed(123)
  s.count <- bs.cell.count$bs_count[bs.cell.count$RNA_predicted_cluster==c]
  
  bs.res <- foreach(i = 1:1000, .combine = rbind) %dopar% {
    bs.data <- c.ne %>%
      group_by(sampleID_3) %>%
      sample_n(s.count, replace = TRUE) %>%
      ungroup()
    
    ttdf.dx <- bs.data %>%
      group_by(condition,stimulation) %>%
      dplyr::summarise(RNA_predicted_activity = mean(RNA_predicted_activity),
                       ATAC_predicted_activity = mean(ATAC_predicted_activity),
                       ATAC_frags_predicted_activity = mean(ATAC_frags_predicted_activity)) %>%
      left_join(int.metadata.filtered2 %>% ungroup() %>% 
                  left_join(samples.meta) %>%
                  distinct(condition, stimulation))
    ttdf.st <- bs.data %>%
      group_by(stimulation) %>%
      dplyr::summarise(RNA_predicted_activity = mean(RNA_predicted_activity),
                       ATAC_predicted_activity = mean(ATAC_predicted_activity),
                       ATAC_frags_predicted_activity = mean(ATAC_frags_predicted_activity)) %>%
      left_join(int.metadata.filtered2 %>% ungroup() %>% 
                  left_join(samples.meta) %>%
                  distinct(stimulation)) %>%
      mutate(condition = "combined")
    ttdf <- rbind(ttdf.dx, ttdf.st)
    
    ### get the null distribution
    null.bs.data <- bs.data %>%
      dplyr::select(5:7) %>%
      mutate(stimulation = sample(c(rep("stimulated", nrow(bs.data)/2), rep("unstimulated", nrow(bs.data)/2))))
    
    null.rna <- c(mean(null.bs.data$RNA_predicted_activity[which(null.bs.data$stimulation=="stimulated")]),
                  mean(null.bs.data$RNA_predicted_activity[which(null.bs.data$stimulation=="unstimulated")]))
    null.atac <- c(mean(null.bs.data$ATAC_predicted_activity[which(null.bs.data$stimulation=="stimulated")]),
                   mean(null.bs.data$ATAC_predicted_activity[which(null.bs.data$stimulation=="unstimulated")]))
    null.atac.frags <- c(mean(null.bs.data$ATAC_frags_predicted_activity[which(null.bs.data$stimulation=="stimulated")]),
                         mean(null.bs.data$ATAC_frags_predicted_activity[which(null.bs.data$stimulation=="unstimulated")]))
    
    ## combine and return
    wide.df <- ttdf %>%ungroup()%>%
      pivot_longer(cols = -c(condition,stimulation)) %>%
      mutate(name2 = paste0(sub("_predicted_activity", "", name), "__", condition,"___",stimulation)) %>%
      select(name2,value) %>%
      column_to_rownames("name2") %>% t()
    ## return
    data.frame(iteration = i,
               wide.df,
               RNA__null___stim = null.rna[1],RNA__null___nostim = null.rna[2],
               ATAC__null___stim = null.atac[1],ATAC__null___nostim = null.atac[2],
               ATAC_frags__null___stim = null.atac.frags[1],ATAC_frags__null___nostim = null.atac.frags[2])
  }
  
  ## combine from all iterations, get mean, and CI
  df33 <- do.call(rbind, lapply(bs.res[,-1], FUN = function(x) mean(x))) %>%
    as.data.frame() %>%
    rownames_to_column("source_1") %>%
    rename("activity" = 2) %>%
    mutate(stimulation = sub(".*___", "", source_1),
           condition = as.character(sub(".*__", "", sub("___.*", "", source_1))),
           type = sub("__.*", "", sub("___.*", "", source_1)),
           celltype = c,
           confin_min = do.call(rbind, lapply(bs.res[,-1], FUN = function(x) quantile(x, 0.025))) %>% as.numeric(),
           confin_max = do.call(rbind, lapply(bs.res[,-1], FUN = function(x) quantile(x, 0.975))) %>% as.numeric(),
           confin_min_90 = do.call(rbind, lapply(bs.res[,-1], FUN = function(x) quantile(x, 0.05))) %>% as.numeric(),
           confin_max_90 = do.call(rbind, lapply(bs.res[,-1], FUN = function(x) quantile(x, 0.95))) %>% as.numeric(),
           confin_min_85 = do.call(rbind, lapply(bs.res[,-1], FUN = function(x) quantile(x, 0.075))) %>% as.numeric(),
           confin_max_85 = do.call(rbind, lapply(bs.res[,-1], FUN = function(x) quantile(x, 0.925))) %>% as.numeric(),
           condition = ifelse(grepl("null", source_1), "all", condition),
           type = ifelse(grepl("null", source_1), sub("null_", "", type), type))
  
  return(df33)
}
write_rds(ne.st.bs, "data/derivatives/neuroestimator-bootstrapping-results-by-st.rds")
ne.st.bs <- read_rds("data/derivatives/neuroestimator-bootstrapping-results-by-st.rds")
################################################################################
## get st data for plotting
ne.st.bs2 <- ne.st.bs %>%
  dplyr::filter(condition != "all") %>%
  pivot_wider(names_from = type, values_from = activity, id_cols = c(celltype, condition,stimulation),
              names_prefix = "mean__") %>%
  inner_join(inner_join(ne.st.bs %>% 
                          dplyr::filter(condition != "all") %>%
                          pivot_wider(names_from = type, values_from = confin_min, id_cols = c(celltype, condition,stimulation),
                                      names_prefix = "confin_min__"),
                        ne.st.bs %>% 
                          dplyr::filter(condition != "all") %>%
                          pivot_wider(names_from = type, values_from = confin_max, id_cols = c(celltype, condition,stimulation),
                                      names_prefix = "confin_max__"))) %>%
  mutate(RNA_sig = (sign(confin_min__RNA)==sign(confin_max__RNA)),
         ATAC_sig = (sign(confin_min__ATAC)==sign(confin_max__ATAC)),
         ATAC_frags_sig = (sign(confin_min__ATAC_frags)==sign(confin_max__ATAC_frags)),
         sig = RNA_sig & ATAC_sig,
         sig_frags = RNA_sig & ATAC_frags_sig)

## get the null distribution eclipse data
ne.st.bs.null <- ne.st.bs %>%
  filter(condition == "all") %>%
  pivot_wider(names_from = type, values_from = activity, id_cols = c(celltype, condition,stimulation),
              names_prefix = "center__") %>%
  inner_join(ne.st.bs %>%
               dplyr::filter(condition == "all") %>%
               mutate(diameter = confin_max -confin_min,
                      diameter_90 = confin_max_90 -confin_min_90,
                      diameter_85 = confin_max_85 -confin_min_85) %>%
               pivot_longer(cols = starts_with("diameter")) %>%
               mutate(type_2 = paste0(name, "__", type)) %>%
               pivot_wider(names_from = type_2, values_from = value, 
                           id_cols = c(celltype, condition,stimulation))) %>%
  select(celltype,stimulation, starts_with(c("center", "diameter")))



#### plot 
# dx.stim.colors <- c()
# names(dx.stim.colors) <- c("epilepsy_baseline", "epilepsy_stimulated", "glioma_baseline", "glioma_stimulated")
ne.st.bs2 %>%
  dplyr::filter(condition != "all") %>%
  left_join(ne.st.bs.null) %>%
  mutate(group = paste0(condition,"_", stimulation),
         condition = factor(condition, levels = c(names(dx.colors),"combined")),
         stimulation = case_when(stimulation =="stim" ~ "stimulated", stimulation == "nostim" ~ "baseline")) %>%
  ggplot(aes(x = mean__ATAC, y = mean__RNA, group = group, color = condition, alpha = stimulation)) +
  geom_point(position = position_dodge(width = 0.6), size =2.5) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.2, color = "red") +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.2, color = "red") +
  geom_errorbarh(aes(xmin = confin_min__ATAC, xmax = confin_max__ATAC),linewidth = 0.8, height = 0, show.legend = F, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = confin_min__RNA, ymax = confin_max__RNA),linewidth = 0.8, width = 0, show.legend = F, position = position_dodge(width = 0.6)) +
  xlim(c(0.42,0.76)) +
  scale_color_manual(values = c(dx.colors,"combined" = stim.colors[1]%>%as.character()), name = "") +
  scale_alpha_manual(values = c("stimulated" = 1,"baseline" = 0.2), name = "") +
  facet_wrap(~celltype, nrow = 2) +
  labs(x = "NEUROeSTIMator predicted activity (ATAC)",
       y = "NEUROeSTIMator predicted activity (RNA)",
       caption = paste0("data points have error bars for the 95% CI\n")) +
  bw.theme + 
  theme(panel.grid.major = element_line(linewidth = 0.05),
        panel.grid.minor = element_line(linewidth = 0.05),
        legend.box = "vertical", legend.spacing.y = unit(0.0005, units = "in")) +
  guides(alpha = guide_legend(nrow = 1), color = guide_legend(nrow = 1))
ggsave2("figs/manuscript_fig4-5-by-st.pdf", width = 9.5, height = 7)

################################################################################
################################################################################
################################################################################


