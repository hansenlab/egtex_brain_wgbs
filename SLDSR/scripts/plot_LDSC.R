# Plot LDSC results for 'adjusting for baseline' analyses
# Peter Hickey
# 2020-05-22

# NOTE: Run with module load conda_R/3.6.x

# Setup ------------------------------------------------------------------------

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

# Brain features ---------------------------------------------------------------

categories <- readRDS(here("SLDSR", "objects", "eGTEx_features.rds"))
# NOTE: Focusing on CG-DMRs, CH-DMRs, and CG-VMRs
categories <- categories[grepl("CG|CH", names(categories)) &
                           grepl("DMR|VMR",  names(categories))]

# Record proportion of SNPs and CpGs in each category --------------------------

categories_df <- bind_rows(
  lapply(
    X = names(categories),
    FUN = function(x) {
      tibble(
        Category = x,
        `Total width (bp)` = sum(width(categories[[x]])),
        `Mean width (bp)` = mean(width(categories[[x]])),
        `Median width (bp)` = median(width(categories[[x]])),
        n = length(categories[[x]]))
    })) %>%
  arrange(`Total width (bp)`)
stopifnot(nrow(categories_df) == length(categories))

snp_prop_table <- bind_rows(
  lapply(
    X = names(categories),
    FUN = function(cn) {
      x <- read_tsv(
        here(
          "SLDSR",
          "output",
          "ldsc",
          paste0(cn, ".Height.results"))) %>%
        filter(Category == "L2_0" | Category == "CNS_0")
      tibble(
        Category = cn,
        `Prop._SNPs` = unlist(x[, "Prop._SNPs"]))
    }))
colnames(snp_prop_table)[2] <- "Proportion of SNPs"

# TODO: Is this actually used anywhere? Should I do the same for CpHs?
cpgs <- findLoci(
  pattern = "CG",
  subject = BSgenome.Hsapiens.UCSC.hg38,
  include = paste0("chr", c(1:22, "X", "Y")),
  strand = "*")
cpgs <- unstrand(cpgs[strand(cpgs) == "+"])
cpg_prop_table <- tibble(
  Category = names(categories),
  `Proportion of CpGs` = sapply(
    X = categories,
    FUN = function(x) {
      sum(overlapsAny(cpgs, x)) / length(cpgs)
    }))

categories_df <- categories_df %>%
  inner_join(snp_prop_table, c("Category" = "Category")) %>%
  inner_join(cpg_prop_table, c("Category" = "Category")) %>%
  arrange(`Total width (bp)`) %>%
  mutate(Category = factor(Category, ordered = TRUE)) %>%
  arrange(Category)

traits_df <- read_csv(here("SLDSR", "tables", "traits_df.csv")) %>%
  mutate(
    N = ifelse(is.na(N_cases), N, N_cases + N_controls),
    TraitType2 = ifelse(Trait == "BMI", "BMI", TraitType)) %>%
  mutate(
    TraitColour = case_when(
      .$TraitType2 == "Additional_phenotype" ~ brewer_pal("qual")(5)[1],
      .$TraitType2 == "Behavioural-cognitive" ~ brewer_pal("qual")(5)[2],
      .$TraitType2 == "Neurological" ~ brewer_pal("qual")(5)[3],
      .$TraitType2 == "Psychiatric" ~ brewer_pal("qual")(5)[4]),
    Trait2Colour = case_when(
      .$TraitType2 == "Additional_phenotype" ~ brewer_pal("qual")(5)[1],
      .$TraitType2 == "Behavioural-cognitive" ~ brewer_pal("qual")(5)[2],
      .$TraitType2 == "Neurological" ~ brewer_pal("qual")(5)[3],
      .$TraitType2 == "Psychiatric" ~ brewer_pal("qual")(5)[4],
      .$TraitType2 == "BMI" ~ brewer_pal("qual")(5)[5]))

# Load data and construct objects ----------------------------------------------

fls <- unlist(
  lapply(
    X = names(categories),
    FUN = function(cn) {
      list.files(
        here("SLDSR", "output", "ldsc"),
        pattern = glob2rx(
          paste0(cn, "*.results")),
        full.names = TRUE)
    }))

# Read in files, tidy up, and rbind
x <- bind_rows(
  lapply(
    X = fls,
    FUN = function(fl) {
      suppressMessages(read_tsv(fl)) %>%
        filter(Category == "L2_0" | Category == "CNS_0") %>%
        mutate(
          Category = sapply(strsplit(basename(fl), "\\."), "[[", 1),
          Trait = sapply(
            strsplit(sub("\\.Phase1\\.results", "", basename(fl)), "\\."), "[[", 2),
          lower = Enrichment - 1.96 * Enrichment_std_error,
          upper = Enrichment + 1.96 * Enrichment_std_error,
          file = fl)
    }))

# Join tidied LDSC output with categories_df
x <- x %>%
  mutate(
    Category = factor(
      x = Category,
      levels = levels(categories_df$Category),
      ordered = TRUE)) %>%
  inner_join(categories_df, by = c("Category" = "Category")) %>%
  inner_join(traits_df, by = c("Trait" = "Trait"))

stopifnot(length(fls) == nrow(x))

# NOTE: Anttila report these traits "had in sufficient evidence of additive
#       heritability for robust analysis" and excluded them from further
#       analysis
x <- x %>%
  filter(Trait != "Agreeableness",
         Trait != "Cardioembolic_stroke",
         Trait != "Large-vessel_disease",
         Trait != "Small-vessel_disease")

# TODO: Remove Conscientiousness and Intracarebral_hemorrhage? Both have huge
#       enrichment SEs.
# x <- x %>%
#   filter(!Trait %in% c("Intracarebral_hemorrhage", "Conscientiousness"))

# Add adjusted P-values
x <- x %>%
  group_by(Category, file) %>%
  filter(grepl(Category, file)) %>%
  ungroup() %>%
  mutate(Coefficient_p = pnorm(`Coefficient_z-score`, lower.tail = FALSE)) %>%
  group_by(Trait) %>%
  mutate(
    Coefficient_holm = p.adjust(Coefficient_p, method = "holm"),
    Enrichment_holm = p.adjust(Enrichment_p, method = "holm"),
    Coefficient_holm_cutoff =
      max(Coefficient_p[Coefficient_holm < 0.05], na.rm = TRUE),
    Enrichment_holm_cutoff =
      max(Enrichment_p[Enrichment_holm < 0.05], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(sig_coef = Coefficient_holm < 0.05) %>%
  arrange(`Total width (bp)`) %>%
  mutate(`Pretty Trait` = gsub("_", " ", Trait)) %>%
  arrange(Category)

# Stratify traits --------------------------------------------------------------

# NOTE: Traits stratified by 'Brain-linked (sig)', 'Brain-linked (non-sig)', or
#       'Non-brain-linked' (stratification defined in 'baseline' analysis)

strata_df <- x %>%
  group_by(Trait) %>%
  summarise(
    strata = ifelse(
      any(Coefficient_holm < 0.05),
      "Brain-linked (sig)",
      "Brain-linked (non-sig)"),
    TraitType2 = unique(TraitType2),
    strata = ifelse(
      TraitType2 == "Additional_phenotype",
      "Non-brain-linked",
      strata)) %>%
  dplyr::select(Trait, strata)
saveRDS(strata_df, here("SLDSR", "objects", "trait_strata_df.rds"))
x_stratified <- inner_join(x, strata_df)

# Tables of results ------------------------------------------------------------

x_stratified %>%
  dplyr::select(
    -lower, -upper, -n, -`Prop._SNPs`, -file, -TraitType2, -TraitColour,
    -Trait2Colour, -Enrichment_p, -Enrichment_holm, -Enrichment_holm_cutoff,
    -`Mean width (bp)`, -`Median width (bp)`, -`Proportion of CpGs`,
    -Coefficient_holm_cutoff, -sig_coef, -Trait, -N, -N_cases, -N_controls,
    -TraitType) %>%
  dplyr::select(
    `Pretty Trait`, strata, Category, `Total width (bp)`,
    `Proportion of SNPs`, starts_with("Coefficient"), everything()) %>%
  rename(
    Feature = Category,
    Trait = `Pretty Trait`,
    `Proportion of h2` = `Prop._h2`,
    `Proportion of h2 standard error` = `Prop._h2_std_error`,
    Stratum = `strata`) %>%
  write_csv(here("SLDSR", "tables", "LDSC_results.baseline_adjustments.csv"))

x_stratified %>%
  dplyr::select(`Pretty Trait`, TraitType, N, N_cases, N_controls) %>%
  distinct() %>%
  rename(Trait = `Pretty Trait`, `Type` = `TraitType`) %>%
  write_csv(here("SLDSR", "tables", "annotated_traits_df.csv"))

# Plot limits ------------------------------------------------------------------

# TODO: Choose these empirically? (currently using BrainEpigenome values)
ylim_coefficient_score <- c(-4.5, 9)
ylim_enrichment <- c(-7, 45)

# Plot results for DMRs --------------------------------------------------------

ct <- grep("VMR", names(categories), value = TRUE, invert = TRUE)

# Coefficient Z-score
g <- x_stratified %>%
  filter(Category %in% ct) %>%
  arrange(Category) %>%
  ggplot(
    data = .,
    aes(
      x = Category,
      y = `Coefficient_z-score`,
      col = Category,
      shape = sig_coef,
      size = sig_coef)) +
  geom_point() +
  facet_wrap( ~ `Pretty Trait`, ncol = 5) +
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
    "Coefficient_Z-score.baseline_adjustments.DMRs.pdf"),
  g,
  height = 6,
  width = 7,
  useDingbats = FALSE)

g <- x_stratified %>%
  filter(Category %in% ct) %>%
  ggplot(
    data = .,
    aes(
      x = Category,
      y = `Coefficient_z-score`,
      col = Category,
      shape = sig_coef,
      size = sig_coef)) +
  geom_jitter(width = 0.3) +
  facet_grid(. ~ strata) +
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
    "Coefficient_Z-score.baseline_adjustments.stratified.DMRs.pdf"),
  g,
  height = 4,
  width = 5,
  useDingbats = FALSE)

# Enrichment
g <- x_stratified %>%
  filter(Category %in% ct) %>%
  filter(strata == "Brain-linked (sig)") %>%
  arrange(Category) %>%
  ggplot(
    data = .,
    aes(
      x = Category,
      y = Enrichment,
      col = Category,
      shape = sig_coef)) +
  geom_point() +
  geom_pointrange(aes(
    ymin = Enrichment - 2 * Enrichment_std_error,
    ymax = Enrichment + 2 * Enrichment_std_error)) +
  facet_wrap( ~ `Pretty Trait`, ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_shape_manual(values = c(1, 16)) +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept = 0, lty = 2) +
  coord_cartesian(ylim = ylim_enrichment)
ggsave(
  here(
    "SLDSR",
    "figures",
    "Enrichment.baseline_adjustments.DMRs.pdf"),
  g,
  height = 6,
  width = 7,
  useDingbats = FALSE)

g <- x_stratified %>%
  filter(Category %in% ct) %>%
  ggplot(
    data = .,
    aes(
      x = Category,
      y = Enrichment,
      col = Category,
      shape = sig_coef,
      size = sig_coef)) +
  geom_jitter(width = 0.2) +
  facet_grid(. ~ strata, labeller = labeller(sig = label_both)) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(ylim = ylim_enrichment)
ggsave(
  here(
    "SLDSR",
    "figures",
    "Enrichment.baseline_adjustments.sig_stratified.DMRs.pdf"),
  g,
  height = 4,
  width = 5,
  useDingbats = FALSE)

# Create legend used in all plot_ldsc.* output
g <- x %>%
  filter(Category %in% ct) %>%
  arrange(Category) %>%
  ggplot(
    data = .,
    aes(
      x = Category,
      y = -log10(Coefficient_p),
      col = Category,
      shape = sig_coef,
      size = sig_coef)) +
  geom_point() +
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3))
legend_plot <- ggdraw(plot_grid(NULL, get_legend(g)))
ggsave(
  here("SLDSR", "figures", "Legend.DMRs.pdf"),
  legend_plot,
  height = 6,
  width = 6,
  useDingbats = FALSE)

# Plot results for VMRs --------------------------------------------------------

ct <- grep("VMR", names(categories), value = TRUE)

# NOTE: Too many levels for a sane colour palette.
colours <- setNames(
  c(Polychrome::palette36.colors(min(length(ct), 36)), "dodgerBlue", "orange"),
  ct)

# Coefficient Z-score
g <- x_stratified %>%
  filter(Category %in% ct) %>%
  arrange(Category) %>%
  ggplot(
    data = .,
    aes(
      x = Category,
      y = `Coefficient_z-score`,
      col = Category,
      shape = sig_coef,
      size = sig_coef)) +
  geom_point() +
  facet_wrap( ~ `Pretty Trait`, ncol = 5) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  scale_colour_manual(values = colours) +
  coord_cartesian(ylim = ylim_coefficient_score)
ggsave(
  here(
    "SLDSR",
    "figures",
    "Coefficient_Z-score.baseline_adjustments.VMRs.pdf"),
  g,
  height = 12,
  width = 14,
  useDingbats = FALSE)

g <- x_stratified %>%
  filter(Category %in% ct) %>%
  ggplot(
    data = .,
    aes(
      x = Category,
      y = `Coefficient_z-score`,
      col = Category,
      shape = sig_coef,
      size = sig_coef)) +
  geom_jitter(width = 0.3) +
  facet_grid(. ~ strata) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  scale_colour_manual(values = colours) +
  coord_cartesian(ylim = ylim_coefficient_score)
ggsave(
  here(
    "SLDSR",
    "figures",
    "Coefficient_Z-score.baseline_adjustments.stratified.VMRs.pdf"),
  g,
  height = 8,
  width = 10,
  useDingbats = FALSE)

# Enrichment
g <- x_stratified %>%
  filter(Category %in% ct) %>%
  filter(strata == "Brain-linked (sig)") %>%
  arrange(Category) %>%
  ggplot(
    data = .,
    aes(
      x = Category,
      y = Enrichment,
      col = Category,
      shape = sig_coef)) +
  geom_point() +
  geom_pointrange(aes(
    ymin = Enrichment - 2 * Enrichment_std_error,
    ymax = Enrichment + 2 * Enrichment_std_error)) +
  facet_wrap( ~ `Pretty Trait`, ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_shape_manual(values = c(1, 16)) +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  scale_colour_manual(values = colours) +
  geom_hline(yintercept = 0, lty = 2) +
  coord_cartesian(ylim = ylim_enrichment)
ggsave(
  here(
    "SLDSR",
    "figures",
    "Enrichment.baseline_adjustments.VMRs.pdf"),
  g,
  height = 12,
  width = 14,
  useDingbats = FALSE)

g <- x_stratified %>%
  filter(Category %in% ct) %>%
  ggplot(
    data = .,
    aes(
      x = Category,
      y = Enrichment,
      col = Category,
      shape = sig_coef,
      size = sig_coef)) +
  geom_jitter(width = 0.2) +
  facet_grid(. ~ strata, labeller = labeller(sig = label_both)) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3)) +
  guides(col = FALSE, shape = FALSE, size = FALSE) +
  scale_colour_manual(values = colours) +
  coord_cartesian(ylim = ylim_enrichment)
ggsave(
  here(
    "SLDSR",
    "figures",
    "Enrichment.baseline_adjustments.sig_stratified.VMRs.pdf"),
  g,
  height = 8,
  width = 10,
  useDingbats = FALSE)

# Create legend used in all plot_ldsc.* output
g <- x %>%
  filter(Category %in% ct) %>%
  arrange(Category) %>%
  ggplot(
    data = .,
    aes(
      x = Category,
      y = -log10(Coefficient_p),
      col = Category,
      shape = sig_coef,
      size = sig_coef)) +
  geom_point() +
  scale_colour_manual(values = colours) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(2, 3))
legend_plot <- ggdraw(plot_grid(NULL, get_legend(g), ncol = 1))
ggsave(
  here("SLDSR", "figures", "Legend.VMRs.pdf"),
  legend_plot,
  height = 12,
  width = 12,
  useDingbats = FALSE)
