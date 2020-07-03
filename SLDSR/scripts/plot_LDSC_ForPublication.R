
setwd("/Volumes/Transcend/eGTEx/SLDSR/")
# Plot LDSC results for 'adjusting for baseline' analyses
# Peter Hickey
# 2019-10-22
# Redone with just DMRs/VMRs of interest by Lindsay Rizzardi 
# after P.Hickey redid analysis with hg38 (so as not to use liftOver)
# 2020-05-27
# Setup ------------------------------------------------------------------------
setwd("/Volumes/Transcend/eGTEx/")
library(GenomicRanges)
library(readr)
library(dplyr)
library(ggplot2)
library(scales)
library(gplots)
library(cowplot)
library(here)
library(bsseq)
library(BSgenome.Hsapiens.UCSC.hg38)

categories <- readRDS(here("SLDSR", "objects", "eGTEx_features.rds"))
# NOTE: Focusing on CG-DMRs, CH-DMRs, and CG-VMRs
categories <- categories[grepl("CG|CH", names(categories)) &
                           grepl("DMR|VMR",  names(categories))]

x_stratified=as_tibble(read.csv(here("SLDSR", "tables", "LDSC_results.baseline_adjustments.csv")))
n=grep("New|DMR|subset",x_stratified$Feature)
x_stratified2=x_stratified[-n,]
x_stratified2=x_stratified2[-grep("Hippocampus",x_stratified2$Feature),]

# Plot limits ------------------------------------------------------------------

# TODO: Choose these empirically
ylim_coefficient_score <- c(-4.5, 9)
ylim_enrichment <- c(-60, 40)

# Plot results for DMRs --------------------------------------------------------
ct <- names(categories)
x_stratified3=x_stratified2 %>% mutate(sig_coef = Coefficient_holm < 0.05)
Feature=x_stratified3$Feature

x_stratified4=x_stratified3 %>% mutate(strata = ifelse(
  any(Coefficient_holm < 0.05),
  "Brain-linked (sig)",
  "Brain-linked (non-sig)"))

# Coefficient Z-score
g <- x_stratified3 %>%
  arrange(Feature) %>%
  ggplot(
    data = .,
    aes(
      x = Feature,
      y = Coefficient_z.score,
      col = Feature,
      shape = sig_coef,
      size = sig_coef)) +
  geom_point() +
  facet_wrap( ~ Trait, ncol = 5) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  scale_color_brewer(palette = "Set1") +
  coord_cartesian(ylim = ylim_coefficient_score)
ggsave(
  here(
    "SLDSR",
    "figures",
    "LR_Coefficient_Z-score.baseline_adjustments.VMRs.hg38.pdf"),
  g,
  height = 6,
  width = 7)

g <- x_stratified3 %>%
  ggplot(
    data = .,
    aes(
      x = Feature,
      y = Coefficient_z.score,
      col = Feature,
      shape = sig_coef,
      size = sig_coef)) +
  geom_jitter(width = 0.3) +
  facet_grid(. ~ Stratum) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  scale_color_brewer(palette = "Set1") +
  coord_cartesian(ylim = ylim_coefficient_score)
ggsave(
  here(
    "SLDSR",
    "figures",
    "LR_Coefficient_Z-score.baseline_adjustments.stratified.VMRs.hg38.pdf"),
  g,
  height = 4,
  width = 5)

# Enrichment
g <- x_stratified3 %>%
  filter(Stratum == "Brain-linked (sig)") %>%
  arrange(Feature) %>%
  ggplot(
    data = .,
    aes(
      x = Feature,
      y = Enrichment,
      col = Feature,
      shape = sig_coef)) +
  geom_point() +
  geom_pointrange(aes(
    ymin = Enrichment - 2 * Enrichment_std_error,
    ymax = Enrichment + 2 * Enrichment_std_error)) +
  facet_wrap( ~ Trait, ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_shape_manual(values = c(1, 16)) +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  scale_color_brewer(palette = "Set1") +
  geom_hline(yintercept = 0, lty = 2) +
  coord_cartesian(ylim = ylim_enrichment)
ggsave(
  here(
    "SLDSR",
    "figures",
    "LR_Enrichment.baseline_adjustments.VMRs.hg38.pdf"),
  g,
  height = 6,
  width = 7)

g <- x_stratified3 %>%
  ggplot(
    data = .,
    aes(
      x = Feature,
      y = Enrichment,
      col = Feature,
      shape = sig_coef,
      size = sig_coef)) +
  geom_jitter(width = 0.2) +
  facet_grid(. ~ Stratum, labeller = labeller(sig = label_both)) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  scale_color_brewer(palette = "Set1")+
                     coord_cartesian(ylim = ylim_enrichment)
                     ggsave(
                       here(
                         "SLDSR",
                         "figures",
                         "LR_Enrichment.baseline_adjustments.sig_stratified.VMRs.hg38.pdf"),
                       g,
                       height = 4,
                       width = 5)
                     
                     # Create legend used in all plot_ldsc.* output
                     g <- x_stratified3 %>%
                       arrange(Feature) %>%
                       ggplot(
                         data = .,
                         aes(
                           x = Feature,
                           y = -log10(Coefficient_p),
                           col = Feature,
                           shape = sig_coef,
                           size = sig_coef)) +
                       geom_point() +
                       scale_color_brewer(palette = "Set1") +
                       scale_shape_manual(values = c(1, 16)) +
                       scale_size_manual(values = c(2, 3))
                     legend_plot <- ggdraw(plot_grid(NULL, get_legend(g)))
                     ggsave(
                       here("SLDSR", "figures", "LR_Legend2.VMRs.hg38.pdf"),
                       legend_plot,
                       height = 6,
                       width = 6)

# Enrichment
ylim_enrichment=c(-5,10)

g <- x_stratified4 %>%
  filter(Stratum == "Brain-linked (sig)") %>%
  arrange(Feature) %>%
  ggplot(
    data = .,
    aes(
      x = Feature,
      y = Enrichment,
      col = Feature,
      shape = sig_coef)) +
  geom_point() +
  geom_pointrange(aes(
    ymin = Enrichment - 2 * Enrichment_std_error,
    ymax = Enrichment + 2 * Enrichment_std_error)) +
  facet_wrap( ~ Trait, ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_shape_manual(values = c(1, 16)) +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  scale_color_brewer(palette = "Set1") +
  geom_hline(yintercept = 0, lty = 2) +
  coord_cartesian(ylim = ylim_enrichment)
ggsave(
  here(
    "SLDSR",
    "figures",
    "LR_Enrichment.baseline_adjustments.VMRs.hg38.pdf"),
  g,
  height = 6,
  width = 7)

g <- x_stratified4 %>%
  ggplot(
    data = .,
    aes(
      x = Feature,
      y = Enrichment,
      col = Feature,
      shape = sig_coef,
      size = sig_coef)) +
  geom_jitter(width = 0.2) +
  facet_grid(. ~ Stratum, labeller = labeller(sig = label_both)) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  scale_color_brewer(palette = "Set1")+
                     coord_cartesian(ylim = c(-50,50))
                     ggsave(
                       here(
                         "SLDSR",
                         "figures",
                         "LR_Enrichment.baseline_adjustments.sig_stratified.VMRs.hg38.pdf"),
                       g,
                       height = 4,
                       width = 5)
                     
                     # Create legend used in all plot_ldsc.* output
                     g <- x_stratified4 %>%
                       arrange(Feature) %>%
                       ggplot(
                         data = .,
                         aes(
                           x = Feature,
                           y = -log10(Coefficient_p),
                           col = Feature,
                           shape = sig_coef,
                           size = sig_coef)) +
                       geom_point() +
                       scale_color_brewer(palette = "Set1") +
                       scale_shape_manual(values = c(1, 16)) +
                       scale_size_manual(values = c(2, 3))
                     legend_plot <- ggdraw(plot_grid(NULL, get_legend(g)))
                     ggsave(
                       here("SLDSR", "figures", "LR_Legend3.VMRs.pdf"),
                       legend_plot,
                       height = 6,
                       width = 6)

# now plot Z.score x Enrichment

ldsc=read.csv("/Users/lrizzardi/Desktop/Transcend/eGTEx/SLDSR/tables/LDSC_results.baseline_adjustments.csv")
library(ggplot2)
ldsc_sig=ldsc[which(ldsc$Coefficient_holm<0.05),]
ldsc_sig=ldsc_sig[-grep("New",ldsc_sig$Feature),]
ldsc_sig_VMR=ldsc_sig[grep("VMR",ldsc_sig$Feature),]
ldsc_sig_VMR=ldsc_sig_VMR[-grep("Hippocampus",ldsc_sig_VMR$Feature),] #because this is all hippo samples before splitting into HC1 and HC2
ldsc_sig_DMR=ldsc_sig[grep("DMR",ldsc_sig$Feature),]
pdf("LDSC_VMR_plot.hg38_2.pdf")
ggplot(ldsc_sig_VMR, aes(x=Coefficient_z.score, y=Enrichment,colour=Feature)) +
  geom_point() + theme_bw() +coord_cartesian(ylim = c(1,14),xlim=c(3,7))+
  geom_text(label=ldsc_sig_VMR$Trait)
  dev.off()
pdf("LDSC_DMR_plot.hg38_2.pdf") 
ggplot(ldsc_sig_DMR, aes(x=Coefficient_z.score, y=Enrichment,colour=Feature)) +
  geom_point() + theme_bw()+coord_cartesian(ylim = c(1,14),xlim=c(3,7))+
  geom_text(label=ldsc_sig_DMR$Trait)
  dev.off()

################
## Now with DMRs
################

x_stratified=as_tibble(read.csv(here("SLDSR", "tables", "LDSC_results.baseline_adjustments.csv")))
n=grep("New|VMR|subset|overall",x_stratified$Feature)
x_stratified2=x_stratified[-n,]
#x_stratified2=x_stratified2[-grep("Hippocampus",x_stratified2$Feature),]

# Plot limits ------------------------------------------------------------------

# TODO: Choose these empirically? (currently using BrainEpigenome values)
ylim_coefficient_score <- c(-4.5, 9)
ylim_enrichment <- c(-5, 15)

# Plot results for DMRs --------------------------------------------------------
ct <- names(categories)
x_stratified3=x_stratified2 %>% mutate(sig_coef = Coefficient_holm < 0.05)




# Coefficient Z-score
g <- x_stratified3 %>%
  arrange(Feature) %>%
  ggplot(
    data = .,
    aes(
      x = Feature,
      y = Coefficient_z.score,
      col = Feature,
      shape = sig_coef,
      size = sig_coef)) +
  geom_point() +
  facet_wrap( ~ Trait, ncol = 5) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(ylim = ylim_coefficient_score)
ggsave(
  here(
    "SLDSR",
    "figures",
    "LR_Coefficient_Z-score.baseline_adjustments.DMRs.hg38.pdf"),
  g,
  height = 6,
  width = 7)

g <- x_stratified3 %>%
  ggplot(
    data = .,
    aes(
      x = Feature,
      y = Coefficient_z.score,
      col = Feature,
      shape = sig_coef,
      size = sig_coef)) +
  geom_jitter(width = 0.3) +
  facet_grid(. ~ Stratum) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(ylim = ylim_coefficient_score)
ggsave(
  here(
    "SLDSR",
    "figures",
    "LR_Coefficient_Z-score.baseline_adjustments.stratified.DMRs.hg38.pdf"),
  g,
  height = 4,
  width = 5)

# Enrichment
g <- x_stratified3 %>%
  filter(Stratum == "Brain-linked (sig)") %>%
  arrange(Feature) %>%
  ggplot(
    data = .,
    aes(
      x = Feature,
      y = Enrichment,
      col = Feature,
      shape = sig_coef)) +
  geom_point() +
  geom_pointrange(aes(
    ymin = Enrichment - 2 * Enrichment_std_error,
    ymax = Enrichment + 2 * Enrichment_std_error)) +
  facet_wrap( ~ Trait, ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_shape_manual(values = c(1, 16)) +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept = 0, lty = 2) +
  coord_cartesian(ylim = c(-5,15))
ggsave(
  here(
    "SLDSR",
    "figures",
    "LR_Enrichment.baseline_adjustments.DMRs.hg38.pdf"),
  g,
  height = 6,
  width = 7)

g <- x_stratified3 %>%
  ggplot(
    data = .,
    aes(
      x = Feature,
      y = Enrichment,
      col = Feature,
      shape = sig_coef,
      size = sig_coef)) +
  geom_jitter(width = 0.2) +
  facet_grid(. ~ Stratum, labeller = labeller(sig = label_both)) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  scale_color_brewer(palette = "Dark2")+
                     coord_cartesian(ylim = ylim_enrichment)
                     ggsave(
                       here(
                         "SLDSR",
                         "figures",
                         "LR_Enrichment.baseline_adjustments.sig_stratified.DMRs.hg38.pdf"),
                       g,
                       height = 4,
                       width = 5)
                     
                     # Create legend used in all plot_ldsc.* output
                     g <- x_stratified3 %>%
                       arrange(Feature) %>%
                       ggplot(
                         data = .,
                         aes(
                           x = Feature,
                           y = -log10(Coefficient_p),
                           col = Feature,
                           shape = sig_coef,
                           size = sig_coef)) +
                       geom_point() +
                       scale_color_brewer(palette = "Set1") +
                       scale_shape_manual(values = c(1, 16)) +
                       scale_size_manual(values = c(2, 3))
                     legend_plot <- ggdraw(plot_grid(NULL, get_legend(g)))
                     ggsave(
                       here("SLDSR", "figures", "LR_Legend2.VMRs.pdf"),
                       legend_plot,
                       height = 6,
                       width = 6)
























### try heatmap...didn't use this
## Heatmap form....
# Zscore
install.packages('pheatmap') # if not installed already
library(pheatmap)
pheatmap(u, display_numbers = T)

meth_col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
t=x_stratified3[,c(1,3,8)] %>% spread(Feature,Coefficient_z.score) 
u=as.matrix(t[,c(2:10)])
rownames(u)=t$Trait


pdf("CG-DMR_OR_neweATAC.pdf")
heatmap.2(u,
          trace = "none",
          col = meth_col_fun,
          breaks=breaks,
          margins = c(10, 12),
          srtCol = 30,
          density.info = "none")

dev.off()

## Heatmap form.... (not used)
# Zscore
# pdf("CG-DMR_OR_neweATAC.pdf")
# heatmap.2(u,
#           trace = "none",
#           col = meth_col_fun,
#           breaks=breaks,
#           margins = c(10, 12),
#           srtCol = 30,
#           density.info = "none")
# 
# dev.off()
# 
# install.packages('pheatmap') # if not installed already
# library(pheatmap)
# pheatmap(u, display_numbers = T)
# 
# meth_col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
# t=x_stratified3[,c(1,3,8)] %>% spread(Feature,Coefficient_z.score) 
# u=as.matrix(t[,c(2:10)])
# rownames(u)=t$Trait
# 
# Feature=x_stratified3$Feature
# tissue=pData(bsseq)$Tissue
# group=pData(bsseq)$group
# ha = HeatmapAnnotation(Trait = Trait, Feature=Feature, 
#                        col = list(type = c("Brain - Amygdala"="#e7298a","Brain - Anterior cingulate cortex (BA24)"="#a6761d","Brain - Caudate (basal ganglia)"="#666666",
#                                            "Brain - Frontal Cortex (BA9)"="#e6ab02","Brain - Hippocampus"="#1b9e77","Brain - Hypothalamus"="#66a61e","Brain - Nucleus accumbens (basal ganglia)"="#d95f02","Brain - Putamen (basal ganglia)"="#7570b3"),
#                                   sex=c("Male"="blue","Female"="red"),group=c("Amygdala"="#e7298a","Cortical"="#a6761d","Hippocampus"="#1b9e77","Basal_ganglia"="#7570b3","Hypothalamus"="#66a61e")))
# 
# 
# pdf("Random20k_general_DMRs_heatmap_dendrogram.pdf")
# Heatmap(as.matrix(gen_DMR_meth), name = "methylation", col = meth_col_fun, 
#         show_row_names = FALSE, show_column_names = FALSE,use_raster=TRUE,raster_device="png") 
# dev.off()


