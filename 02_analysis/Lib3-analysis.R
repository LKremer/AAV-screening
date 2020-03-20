library(tidyverse)
library(pheatmap)
library(viridis)
library(stringdist)
library(matrixStats)
library(magrittr)
library(ggbeeswarm)
require(svglite)

theme_set(theme_minimal())

#' Define all input files:
bc_count_file <- "../01_counting-barcodes/AAV-lib3/output/merged-barcode-counts.tsv.gz"
bc_info_file <- "Lib3_Overview_Sascha.csv"
sample_info_file <- "Lib3_cell_counts.tsv"

#' Read all files:
sample_info <- read_tsv( sample_info_file, col_types = "cccil")


raw_bc_cnt <- bc_count_file %>%
  read_tsv( col_names = c("barcode_count", "barcode_seq", "barcode_direction", "sample_ID", "read_mate",
                          "set", "celltype", "sample_n", "fastq_read_count"),
            col_types = "iccccccc--i") %>%
  mutate( set = case_when(
    set == "In vitro" & sample_n == "Sample1" ~ "In vitro 1",
    set == "In vitro" & sample_n == "Sample2" ~ "In vitro 2",
    TRUE ~ set)) %>%
  left_join( sample_info, by = c("set", "celltype", "sample_n")) %>%
  filter( !ignore) %>% dplyr::select( -n_facs_events, -ignore)  # remove some hand-selected samples, marked with ignore=TRUE
# A pipetting error occurred in these few tubes, hence we're excluding them from the analysis to be on the safe side

bc_info <- read_csv( bc_info_file,
                     col_names = c("AAV_ID", "barcode_seq", "AAV_lib_read_n"),
                     col_types = "-c-c-i-", skip = 1) %>%
  mutate( AAV_lib_scalefactor = mean(AAV_lib_read_n) / AAV_lib_read_n )




#' Since we got some sequencing errors:
#' Take all identified barcodes and find the most similar actual barcode (Levenshtein distance < 3)
nearest_bc <- tibble( unknown_barcode = unique(raw_bc_cnt$barcode_seq)) %>%
  mutate( matched_barcode = bc_info$barcode_seq[ amatch( unknown_barcode, bc_info$barcode_seq, maxDist = 2, method = "lv")] ) %>%
  filter( !is.na(matched_barcode))

#' Set 15-mers for which we can't find a known barcode within a Levenshtein distance of 2 to NA.
#' Change all barcodes where we do find a match to the matching barcode.
raw_bc_cnt_match <- raw_bc_cnt %>%
  left_join( nearest_bc, by = c("barcode_seq" = "unknown_barcode")) %>%
  dplyr::select( -barcode_seq, barcode_seq = matched_barcode)

#' merge all barcode counts within tubes (=technical replicates)
bc_cnt <- raw_bc_cnt_match %>%
  group_by( set, celltype, sample_n, barcode_seq) %>%
  summarise( barcode_count = sum(barcode_count)) %>% ungroup %>%
  filter( !is.na(barcode_seq)) %>%  #' remove the few barcodes we couldn't interpret cause of sequencing errors
  left_join( sample_info, by = c("set", "celltype", "sample_n")) %>%
  mutate( barcode_count = barcode_count * (n_facs_events/500) ) %>%  # downweight the tubes with < 500 cells cause they are noisier
  group_by( set, celltype, barcode_seq) %>%
  summarize( barcode_count = sum(barcode_count)) %>% ungroup

#' Sometimes an AAV has 0 barcode counts in a sample. In that case, it is not reported in the barcode count input file,
#' so we have to add the zero-counts manually (cause we still need the zeros to calculate e.g. the mean) 
not_observed <- bc_cnt %>% distinct(celltype, set) %>% unite(combo, sep="&") %>% pull(combo) %>%
  expand.grid(unique(bc_cnt$barcode_seq), stringsAsFactors=F) %>% as_tibble %>%
  separate(Var1, c("celltype", "set"), sep="&") %>% dplyr::rename(barcode_seq=Var2) %>%
  mutate( barcode_count = 0) %>% anti_join( bc_cnt, by = c("set", "celltype", "barcode_seq"))

#' Add the zero-count barcodes and find the AAV name that corresponds to each barcode:
bc_cnt <- bc_cnt %>%
  bind_rows( not_observed) %>%
  left_join( bc_info, by = "barcode_seq")

#' Calculate the raw proportions and the normalised (by input) proportions
bc_props_all <- bc_cnt %>% 
  mutate( barcode_count_scaled = barcode_count * AAV_lib_scalefactor) %>%
  group_by( set, celltype) %>%  # here, the tubes (~500 cells each) are merged
  mutate( total_barcode_count = sum(barcode_count),
          total_barcode_count_scaled = sum(barcode_count_scaled)) %>% ungroup %>%
  group_by( set, celltype, barcode_seq, AAV_ID) %>%
  mutate( proportion = barcode_count / total_barcode_count,
          norm_proportion = barcode_count_scaled / total_barcode_count_scaled) %>% ungroup

bc_props_invitro <- bc_props_all %>% filter( str_detect(set, "vitro"))
bc_props         <- bc_props_all %>% filter(!str_detect(set, "vitro"))

#' Sanity check: the barcode proportions (both normalised and raw) all add up to one!
bc_props %>% group_by( set, celltype) %>%
  summarise( total_prop = sum(proportion), total_norm_prop = sum(norm_proportion))

#' Check if our normalisation worked:
#' (Some AAVs were overrepresented in the AAV cocktail so they had an unfair advantage)
bc_props %>%  # no normalisation: overrepresented AAVs have an advantage
  ggplot( aes( x = proportion, y = AAV_lib_read_n)) +
  geom_point( size = .1) + geom_smooth( method = "lm", se = F) + scale_y_log10() + scale_x_log10()
bc_props %>%  # after normalisation: the effect is gone
  ggplot( aes( x = norm_proportion, y = AAV_lib_read_n)) +
  geom_point( size = .1) + geom_smooth( method = "lm", se = F) + scale_y_log10() + scale_x_log10()


#' Calculate the median normalised barcode proportion + the inter-quartile range (IQR)
summary_stats <- bc_props %>%
  group_by( AAV_ID) %>%
  summarize( mean = mean( norm_proportion),
             median = median( norm_proportion),
             IQR_lower = quantile( norm_proportion, 1/4),
             IQR_upper = quantile( norm_proportion, 3/4)) %>% ungroup

#' Get the top 10 AAVs according to median normalised barcode proportion:
topX <- summary_stats %>%
  arrange( -median) %>% pull( AAV_ID) %>% head(10)
topX_mean <- summary_stats %>%
  arrange( -mean) %>% pull( AAV_ID) %>% head(10)


#' # Plot the barplots
aav_colors <- rep("lightgray", length(bc_props$AAV_ID)) %>% set_names(bc_props$AAV_ID)
aav_colors_dark <- rep("black", length(bc_props$AAV_ID)) %>% set_names(bc_props$AAV_ID)
aav_colors[c("AAV1_P5", "AAV9_A2", "AAV9_wt", "AAV2_wt")] <- c("#800033ff", "#56b4e9ff", "#d55e00ff", "#009e73ff")
aav_colors_dark[c("AAV1_P5", "AAV9_A2", "AAV9_wt", "AAV2_wt")] <- c("#800033ff", "#56b4e9ff", "#d55e00ff", "#009e73ff")

#' Plot barplots of overall top 10 + data points (sorted by median)
bc_props %>% filter( AAV_ID %in% topX) %>%
  ggplot( aes( x = fct_relevel(AAV_ID, topX), y = norm_proportion, fill = AAV_ID)) +
  geom_boxplot( outlier.size = 0, outlier.alpha = 0) +
  geom_quasirandom( size = .75, color = "white") + 
  geom_quasirandom( size = .5) +
  theme( axis.text.x = element_text( angle = 30, hjust = 1)) +  # rotate x-axis labels
  labs( y = "normalized proportion", x = NULL) + guides( color=F, fill=F) + 
  scale_fill_manual( values = aav_colors) +
  ylim( 0, .5) +
  geom_point( data = bc_props %>% filter( AAV_ID %in% topX, norm_proportion > .5),
              mapping = aes( x = AAV_ID), y = .51, shape = 17, size = 1) +
  theme( panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
ggsave("plots/AAVlib3_overall-top10-boxplots.png", width=5, height=3, dpi=400)
ggsave("plots/AAVlib3_overall-top10-boxplots.svg", width=5, height=3, dpi=400)

#' Plot barplots of overall top 10 + data points (sorted by median)
bc_props %>% filter( AAV_ID %in% topX) %>%
  ggplot( aes( x = fct_relevel(AAV_ID, topX), y = norm_proportion, fill = AAV_ID)) +
  geom_col( data = filter( summary_stats, AAV_ID %in% topX),
            mapping = aes( x = fct_relevel(AAV_ID, topX), y = median, fill = AAV_ID)) +
  geom_quasirandom( size = .75, color = "white") + 
  geom_quasirandom( size = .5) +
  theme( axis.text.x = element_text( angle = 30, hjust = 1)) +  # rotate x-axis labels
  labs( y = "normalized proportion", x = NULL) + guides( color=F, fill=F) + 
  scale_fill_manual( values = aav_colors) +
  scale_y_continuous(limits=c(0, .5)) +
  geom_point( data = bc_props %>% filter( AAV_ID %in% topX, norm_proportion > .5),
              mapping = aes( x = AAV_ID), y = .51, shape = 17, size = 1) +
  theme( panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
ggsave("plots/AAVlib3_overall-top10-barplots-median.png", width=5, height=3, dpi=400)
ggsave("plots/AAVlib3_overall-top10-barplots-median.svg", width=5, height=3, dpi=400)

#' Plot boxplots of overall top 10 + data points (sorted by mean)
bc_props %>% filter( AAV_ID %in% topX_mean) %>%
  ggplot( aes( x = fct_relevel(AAV_ID, topX_mean), y = norm_proportion, fill = AAV_ID)) +
  geom_col( data = filter( summary_stats, AAV_ID %in% topX_mean),
            mapping = aes( x = fct_relevel(AAV_ID, topX_mean), y = mean, fill = AAV_ID)) +
  geom_quasirandom( size = .75, color = "white") + 
  geom_quasirandom( size = .5) +
  theme( axis.text.x = element_text( angle = 30, hjust = 1)) +  # rotate x-axis labels
  labs( y = "normalized proportion", x = NULL) + guides( color=F, fill=F) + 
  scale_fill_manual( values = aav_colors) +
  scale_y_continuous(limits=c(0, .5)) +
  geom_point( data = bc_props %>% filter( AAV_ID %in% topX_mean, norm_proportion > .5),
              mapping = aes( x = AAV_ID), y = .51, shape = 17, size = 1) +
  theme( panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
ggsave("plots/AAVlib3_overall-top10-barplots-mean.png", width=5, height=3, dpi=400)
ggsave("plots/AAVlib3_overall-top10-barplots-mean.svg", width=5, height=3, dpi=400)

#' Plotting by celltype:
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}

scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}  

topX_by_celltype <- bc_props %>% group_by( celltype, AAV_ID) %>% summarise( mean = mean(norm_proportion)) %>%
  arrange( celltype, -mean) %>% mutate( rank = row_number()) %>% ungroup %>% filter( rank <= 10) %>%
  unite( topX, celltype, AAV_ID) %>% pull(topX)

bc_props_celltype <- bc_props %>%
  unite( combination, celltype, AAV_ID, remove = F) %>% filter( combination %in% topX_by_celltype)

bc_props_celltype_mean <- bc_props_celltype %>%
  group_by( celltype, AAV_ID) %>% summarise( mean = mean(norm_proportion)) %>% ungroup

bc_props_celltype %>%
  ggplot( aes( x = reorder_within( AAV_ID, -norm_proportion, celltype), y = norm_proportion, color = AAV_ID)) +
  geom_col( data = bc_props_celltype_mean,
            mapping = aes( x = reorder_within( AAV_ID, -mean, celltype), y = mean, fill = AAV_ID)) +
  geom_quasirandom( size = 1.5, color = "white", width = .2) + 
  geom_quasirandom( size = 1, color = "black", width = .2) +
  geom_point( data = bc_props %>% filter( AAV_ID %in% topX_mean, norm_proportion > .57),
              y = .57, shape = 17, size = 1.5, color = "black") +
  scale_x_reordered() +
  facet_wrap( ~ celltype, scales = "free_x") +
  theme( axis.text.x = element_text( angle = 45, hjust = 1)) +
  theme( panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual( values = aav_colors) + scale_color_manual( values = aav_colors) +
  guides( fill=F, color=F) +
  scale_y_continuous(limits = c(0, .56)) +
  labs( y = "normalized proportion", x = NULL) + guides( color=F, fill=F)
ggsave("plots/AAVlib3_top10-perCelltype-barplots-mean.png", width=10, height=6, dpi=400)
ggsave("plots/AAVlib3_top10-perCelltype-barplots-mean.svg", width=10, height=6, dpi=400)


topX_invitro <- bc_props_invitro %>%
  group_by(AAV_ID) %>% summarise( mean = mean(norm_proportion)) %>% arrange(-mean) %>% head(10)

#' Plot in vitro sample:
bc_props_invitro %>%
  filter( AAV_ID %in% pull(topX_invitro, AAV_ID)) %>%
  ggplot( aes( x = fct_relevel(AAV_ID, pull(topX_invitro, AAV_ID)), y = norm_proportion, fill = AAV_ID)) +
  geom_col( data = topX_invitro,
            mapping = aes( x = fct_relevel(AAV_ID, pull(topX_invitro, AAV_ID)), y = mean, fill = AAV_ID)) +
  geom_point( size = .75, color = "white") + 
  geom_point( size = .5) +
  theme( axis.text.x = element_text( angle = 30, hjust = 1)) +  # rotate x-axis labels
  labs( y = "normalized proportion", x = NULL) + guides( color=F, fill=F) + 
  scale_fill_manual( values = aav_colors) +
  theme( panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
ggsave("plots/AAVlib3_inVitro-top10-barplots.png", width=5, height=3, dpi=400)
ggsave("plots/AAVlib3_inVitro-top10-barplots.svg", width=5, height=3, dpi=400)



#' Export all our info as CSV files:
bc_props %>% bind_rows( bc_props_invitro) %>%
  dplyr::rename(`AAV_nreads_in_transduction_library`=AAV_lib_read_n, AAV_normalisation_factor=AAV_lib_scalefactor,
                normalised_barcode_count=barcode_count_scaled, total_normalised_barcode_count=total_barcode_count_scaled,
                raw_proportion=proportion, normalised_proportion=norm_proportion) %>%
  arrange( celltype, set, -normalised_proportion) %>%
  write_csv("tables/AAVlib3_AAV-barcode-counts.csv")

bc_props %>%
  bind_rows( bc_props_invitro %>% mutate( celltype = paste0("in_vitro_", celltype)) ) %>%
  group_by( celltype, AAV_ID) %>%
  summarize( n = n(),
             mean = mean( norm_proportion, na.rm=T),
             sd = sd( norm_proportion, na.rm=T),
             median = median( norm_proportion, na.rm=T),
             IQR_lower = quantile( norm_proportion, 1/4, na.rm=T),
             IQR_upper = quantile( norm_proportion, 3/4, na.rm=T)) %>% ungroup %>%
  arrange( celltype, -mean) %>%
  write_csv("tables/AAVlib3_AAV-barcode-proportions-by-celltype.csv")

bc_props %>%
  group_by( AAV_ID) %>%
  summarize( n = n(),
             mean = mean( norm_proportion),
             sd = sd( norm_proportion),
             median = median( norm_proportion),
             IQR_lower = quantile( norm_proportion, 1/4),
             IQR_upper = quantile( norm_proportion, 3/4)) %>% ungroup %>%
  arrange( -mean) %>%
  write_csv("tables/AAVlib3_AAV-barcode-proportions_InVivo.csv")

bc_props_invitro %>%
  group_by( AAV_ID) %>%
  summarize( n = n(),
             mean = mean( norm_proportion),
             sd = sd( norm_proportion),
             median = median( norm_proportion),
             IQR_lower = quantile( norm_proportion, 1/4),
             IQR_upper = quantile( norm_proportion, 3/4)) %>% ungroup %>%
  arrange( -mean) %>%
  write_csv("tables/AAVlib3_AAV-barcode-proportions_InVitro.csv")

#' Write supplementary table
for (ct in unique(bc_props_all$celltype)) {
  if (ct == "NSC") {
    outf <- sprintf("tables/supplement/AAVlib3_%s_inVitro_serotype_proportions.tsv", ct)
  } else {
    outf <- sprintf("tables/supplement/AAVlib3_%s_serotype_proportions.tsv", ct)
  }
  bc_props_all %>% filter( celltype == ct) %>%
    dplyr::select( celltype, set, AAV_ID, norm_proportion) %>%
    group_by( celltype, AAV_ID) %>%
    mutate( average = mean( norm_proportion, na.rm=T)) %>% ungroup %>%
    spread( key = set, value = norm_proportion) %>% arrange(-average) %>%
    dplyr::rename( population=celltype, serotype=AAV_ID) %>%
    write_tsv( outf)
}

bc_props %>%
  unite( population, celltype, set, sep = " ") %>%
  dplyr::select( population, AAV_ID, norm_proportion) %>%
  group_by( AAV_ID) %>%
  mutate( average = mean( norm_proportion, na.rm=T)) %>% ungroup %>%
  spread( key = population, value = norm_proportion) %>% arrange(-average) %>%
  dplyr::select( serotype=AAV_ID, average, starts_with("aN"), starts_with("qN"), starts_with("TAP"),
                 starts_with("NB"), starts_with("AS"), starts_with("OL"), starts_with("EP")) %>%
  write_tsv( "tables/supplement/AAVlib3_overall_serotype_proportions.tsv")
