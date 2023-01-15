# Workflow for processing SIP data

library(tidyverse)
library(gt)
library(ggplot2)
library(ggtext)
library(patchwork)
library(qSIP)
library(jakR)

plot_theme = jak_theme(plot_title_size = 14, 
                    axis_title_size = 13,
                    base_size = 12)

table_theme = function(data) {
  data %>%
  gt_theme(background_color = "#037bcf35",
           font_color = "#000000",
           table_font_size = 14,
           column_label_size = 15,
           align = "left",
           table_padding = 2)
}

metadata_file <- "./processing/zotu_table.csv"
taxa_file <- "./processing/taxonomy_legend.csv"
sdat_file <- "./processing/copies_per_g_soil.csv"

t7_file <- "./processing/amplicon_Treatment_comparison-7.txt"
t30_file <- "./processing/amplicon_Treatment_comparison-30.txt"

qsip_df <- read_tsv(metadata_file)

qsip_df <- qsip_df %>%
  mutate(unique.tube = paste(Replicate, Time, Isotope, sep="_")) %>%
  relocate(unique.tube, .after = Trt.code)


qsip_df <- qsip_df %>%
  mutate(unique.tmt = paste(Time, Isotope, sep="_")) %>%
  relocate(unique.tmt, .after = unique.tube)


qsip_df <- qsip_df %>% 
  mutate(sum.abundance = rowSums(select(., starts_with("X")))) %>%
  relocate(sum.abundance, .after = unique.tmt)

qsip_df <- qsip_df %>% 
  pivot_longer(cols = starts_with("X"), names_to = "taxon", values_to = "rel.abundance")

qsip_df <- qsip_df %>%
  mutate(t.copies.ul = copies.ul * rel.abundance) %>%
  select(-sum.abundance)

#Taxonomy file

qsip_taxa <- read_tsv(taxa_file, col_names = c("code", "taxon"))

#setdiff(qsip_taxa$taxon, unique(qsip_df$taxon)) # should be 0
#setdiff(unique(qsip_df$taxon), qsip_taxa$taxon) # should be 0

qsip_df <- qsip_df %>%
  left_join(qsip_taxa, by = "taxon")

T7 = read_tsv(t7_file)
T30 = read_tsv(t30_file)

copies_tube = qsip_df %>%
  group_by(taxon, unique.tube) %>%
  summarize(Replicate = Replicate, Time = Time, Isotope = Isotope, tube = tube, unique.tmt = unique.tmt, copies.ul = sum(copies.ul), t.copies.ul = sum(t.copies.ul), code = code, .groups = "drop") %>%
  unique() %>%
  select(c("taxon", "tube", "Replicate", "Time", "Isotope", "unique.tube", 
           "unique.tmt", "copies.ul", "t.copies.ul", "code"))

#First round of filtering
qsip_t7 = filter_taxa(qsip_df,
                      treatment1 = NULL, 
                      treatment2 = T7$trt.code.2[1], 
                      treatment_refs = NULL, 
                      min_reps = 3,
                      min_fractions = 1, #
                      taxon_column = "taxon",
                      filter_column = "t.copies.ul",
                      treatment_column = "unique.tmt")

qsip_t7$plot + plot_theme + labs(title = "First pass filter")

#Second round of filtering
qsip_t7 = filter_taxa(qsip_t7$df, 
                         treatment1 = NULL, 
                         treatment2 = T7$trt.code.1[1], 
                         treatment_refs = T7$trt.refs[1], 
                         min_reps = 3,
                         min_fractions = 1,
                         filter_column = "t.copies.ul", 
                         treatment_column = "unique.tmt")

qsip_t7$plot + plot_theme +
  labs(title = "Second pass filter (log scale)") + scale_y_log10()

#T30 filtering
qsip_t30 <- filter_taxa(qsip_df,
                      treatment1 = NULL,
                      treatment2 = T30$trt.code.2[1],
                      treatment_refs = NULL,
                      min_reps = 3,
                      min_fractions = 1,
                      taxon_column = "taxon",
                      filter_column = "t.copies.ul",
                      treatment_column = "unique.tmt")

qsip_t30$plot + plot_theme + labs(title = "First pass filter")

qsip_t30 <- filter_taxa(qsip_t30$df,
                       treatment1 = NULL,
                       treatment2 = T30$trt.code.1[1],
                       treatment_refs = T30$trt.refs[1],
                       min_reps = 3,
                       min_fractions = 1,
                       filter_column = "t.copies.ul",
                       treatment_column = "unique.tmt")

qsip_t30$plot + plot_theme +
  labs(title = "Second pass filter") + scale_y_log10()


sdat <- read_tsv(sdat_file)

qsip_t7_time <- qsip_t7.comparisons <- all_taxa_calcs(X.all = qsip_t7,
                                output_directory = "./processing/rScript/t7/",
                                comparisons = T7,
                                M.soil = sdat,
                                taxon_column = "taxon",
                                density_column = "Density",
                                copies_ul_column = "t.copies.ul",
                                tube_column = "unique.tube",
                                treatment_column = "unique.tmt",
                                soil_g_column = "g.dry.soil.tube",
                                growth.model = "linear",
                                prop.O.from.water = 0.60,
                                v.frac = 1,
                                copies.cell = 6,
                                pgC.cell = 0.1,
                                CI = 0.90,
                                draws = 1000, #Lowers run time
                                tailed.test = 1)

qsip_t7.comparisons %>% as_tibble() %>%
  select(taxonID, wad.diff.p.value, starts_with("ape")) %>%
  #filter(ape.boot.CI.L > 0) %>%
  ggplot(aes(x = reorder(taxonID, ape.obs),
             y = ape.boot.mean,
             color = wad.diff.p.value)) +
    plot_theme +
    geom_point() +
    geom_errorbar(aes(ymin = ape.boot.CI.L, ymax = ape.boot.CI.U)) +
    labs(title = "T7", x = "Taxa (reordered by ape.obs)") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
          legend.position = "right")
#Save the figure
ggsave("./processing/rScript/wad_diff_T7.svg", width = 8, height = 6)


qsip_t30_time = qsip_t30.comparisons <- all_taxa_calcs(X.all = qsip_t30,
                                output_directory = "./processing/rScript/t30/",
                                comparisons = T30,
                                M.soil = sdat,
                                taxon_column = "taxon",
                                density_column = "Density",
                                copies_ul_column = "t.copies.ul",
                                tube_column = "unique.tube",
                                treatment_column = "unique.tmt",
                                soil_g_column = "g.dry.soil.tube",
                                growth.model = "linear",
                                prop.O.from.water = 0.60,
                                v.frac = 1,
                                copies.cell = 6,
                                pgC.cell = 0.1,
                                CI = 0.90,
                                draws = 1000, #Lowers run time
                                tailed.test = 1)

qsip_t30.comparisons %>% as_tibble() %>%
  select(taxonID, wad.diff.p.value, starts_with("ape")) %>%
  #filter(ape.boot.CI.L > 0) %>%
  ggplot(aes(x = reorder(taxonID, ape.obs),
             y = ape.boot.mean,
             color = wad.diff.p.value)) +
    plot_theme +
    geom_point() +
    geom_errorbar(aes(ymin = ape.boot.CI.L, ymax = ape.boot.CI.U)) +
    labs(title = "T30", x = "Taxa (reordered by ape.obs)") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
          legend.position = "right")
#Save the figure
ggsave("./processing/rScript/wad_diff_T30.svg", width = 8, height = 6)


df = bind_rows("T7" = qsip_t7.comparisons, 
               "T30" = qsip_t30.comparisons, 
               .id = "COMP") %>% 
  as_tibble()

df %>%
  #filter(ape.boot.CI.L > 0) %>%
  #filter(wad.diff.p.value < 0.1) %>%
  select(COMP, taxonID, ape.obs) %>%
  pivot_wider(names_from = COMP, values_from = ape.obs, values_fill = 0) %>%
  left_join(qsip_taxa, by = c("taxonID" = "taxon")) %>%
  separate(code, sep = ";", into = c("k", "p", "c", "o", "f", "g")) %>%
  ggplot(aes(x = T7, y = T30, fill = p)) +
    plot_theme +
    geom_point(pch = 21) +
    facet_wrap(~p) +
    scale_fill_manual(values = palette_jak$bay(28)) +
    theme(legend.position = "none") +
    labs(title = "ape.obs in T7 and T30, faceted by phylum")
#Save the figure
ggsave("./processing/rScript/ape_obs_T7_T30_phylum.svg",
  width = 10, height = 10, units = "in", dpi = 300)


#APE obs by phylum
qsip_t30.comparisons %>% as_tibble() %>%
  select(taxonID, wad.diff.p.value, starts_with("ape")) %>%
  left_join(qsip_taxa, by = c("taxonID" = "taxon")) %>%
  separate(code, sep = ";", into = c("k", "p", "c", "o", "f", "g")) %>%
  arrange(p, ape.obs) %>%
  mutate(taxonID = factor(taxonID, levels = taxonID)) %>%
  ggplot(aes(x = taxonID,
             y = ape.obs,
             color = p)) +
    plot_theme +
    geom_point() +
    geom_errorbar(aes(ymin = ape.boot.CI.L, ymax = ape.boot.CI.U)) +
    labs(title = "T30", x = "Taxa (reordered by ape.obs)") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = "right")
#Save the figure
ggsave("./processing/rScript/wad_diff_T30_phylum.svg")

#APE obs by phylum
qsip_t7.comparisons %>% as_tibble() %>%
  select(taxonID, wad.diff.p.value, starts_with("ape")) %>%
  left_join(qsip_taxa, by = c("taxonID" = "taxon")) %>%
  separate(code, sep = ";", into = c("k", "p", "c", "o", "f", "g")) %>%
  arrange(p, ape.obs) %>%
  mutate(taxonID = factor(taxonID, levels = taxonID)) %>%
  ggplot(aes(x = taxonID,
             y = ape.obs,
             color = p)) +
    plot_theme +
    geom_point() +
    geom_errorbar(aes(ymin = ape.boot.CI.L, ymax = ape.boot.CI.U)) +
    labs(title = "T7", x = "Taxa (reordered by ape.obs)") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = "right")
#Save the figure
ggsave("./processing/rScript/wad_diff_T7_phylum.svg", width = 14, height = 10)

#APE obs by phylum significant
qsip_t30.comparisons %>% as_tibble() %>%
  select(taxonID, wad.diff.p.value, starts_with("ape")) %>%
  left_join(qsip_taxa, by = c("taxonID" = "taxon")) %>%
  separate(code, sep = ";", into = c("k", "p", "c", "o", "f", "g")) %>%
  filter(ape.boot.CI.L > 0) %>%
  filter(wad.diff.p.value < 0.05) %>%
  mutate(p = gsub("p__", "", p)) %>%
  arrange(p, ape.obs) %>%
  mutate(taxonID = factor(taxonID, levels = taxonID)) %>%
  ggplot(aes(x = taxonID,
             y = ape.obs,
             color = p)) +
    facet_grid(. ~ p, scales = "free_x", space = "free_x", 
      labeller = label_wrap_gen(width = 10, multi_line = TRUE)) +
    plot_theme +
    geom_point() +
    geom_errorbar(aes(ymin = ape.boot.CI.L, ymax = ape.boot.CI.U)) +
    labs(title = "T30", x = "Taxa (reordered by ape.obs)") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = "right")
#Save the figure
ggsave("./processing/rScript/wad_diff_T30_phylum-005.svg", width = 14, height = 10)



#Numbers of growing taxa
qsip_t30.comparisons %>% as_tibble() %>%
          select(taxonID, wad.diff.p.value, starts_with("ape")) %>%
          left_join(qsip_taxa, by = c("taxonID" = "taxon")) %>%
          separate(code, sep = ";", into = c("k", "p", "c", "o", "f", "g")) %>%
          filter(ape.boot.CI.L > 0) %>%
          filter(wad.diff.p.value < 0.1)

qsip_t7.comparisons %>% as_tibble() %>%
          select(taxonID, wad.diff.p.value, starts_with("ape")) %>%
          left_join(qsip_taxa, by = c("taxonID" = "taxon")) %>%
          separate(code, sep = ";", into = c("k", "p", "c", "o", "f", "g")) %>%
          filter(ape.boot.CI.L > 0) %>%
          filter(wad.diff.p.value < 0.1)