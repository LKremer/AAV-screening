library(tidyverse)
library(patchwork)

theme_set(theme_bw() + theme(text = element_text(family = 'Arial'),
                             panel.grid.minor = element_blank(),    # remove axis grid
                             panel.grid.major = element_blank(),
                             strip.background = element_blank(),    # remove silly gray bg of facet title
                             strip.text = element_text(hjust = 0),  # left-justify facets
                             plot.tag = element_text(size = 20, face = "bold")
))

# setwd("/home/lukas/code/AAV_screening/02_analysis/")

# get the names of the 68 shared AAV capsids
shared_aavs <- intersect(
  unique(pull(read_csv("tables/AAVlib1_AAV-barcode-counts.csv", col_types = "----c-------"))),
  unique(pull(read_csv("tables/AAVlib3_AAV-barcode-counts.csv", col_types = "----c-------") %>% mutate(AAV_ID = str_replace(AAV_ID, "wt", "_WT")))))

# get all shared AAVs of library 1 and calculate the proportion of AAV-barcodes
# in every sample. For every AAV, take the mean of this proportion.
lib1 <- read_csv("tables/AAVlib1_AAV-barcode-counts.csv", col_types = "cccdcdd-----") %>%
  filter(AAV_ID %in% shared_aavs) %>%
  group_by(set, celltype) %>%
  mutate(total_bc = sum(barcode_count)) %>% ungroup() %>%
  group_by(set, celltype, AAV_ID) %>%
  summarise(fraction = barcode_count / total_bc) %>% ungroup() %>%
  group_by(AAV_ID) %>%
  summarise(mean_fraction = mean(fraction)) %>% ungroup()

# same as above but with library 3
lib3 <- read_csv("tables/AAVlib3_AAV-barcode-counts.csv", col_types = "cccdcdd-----") %>%
  mutate(AAV_ID = str_replace(AAV_ID, "wt", "_WT")) %>%
  filter(AAV_ID %in% shared_aavs) %>%
  filter(!str_detect(set, "vitro")) %>%
  #mutate(barcode_count = barcode_count * AAV_normalisation_factor) %>% filter(celltype %in% c("NSC", "aNSC", "qNSC")) %>%
  group_by(set, celltype) %>%
  mutate(total_bc = sum(barcode_count)) %>% ungroup() %>%
  group_by(set, celltype, AAV_ID) %>%
  summarise(fraction = barcode_count / total_bc) %>% ungroup() %>%
  group_by(AAV_ID) %>%
  summarise(mean_fraction = mean(fraction)) %>% ungroup()

# merge the two tables
shared <- full_join(lib1, lib3, by = "AAV_ID", suffix = c("_lib1", "_lib3"))

# define a color scheme for our candidate AAVs
candidates <- shared %>%
  filter(AAV_ID %in% c("AAV1_P5", "AAV9_A2", "AAV9_WT", "AAV2_WT")) %>%
  mutate(AAV_ID = fct_relevel(str_replace(AAV_ID, "_WT", "_wt"), "AAV2_wt", "AAV9_wt", "AAV9_A2"))
aav_colors <- c("other" = "gray50")
aav_colors[c("AAV1_P5", "AAV9_A2", "AAV9_WT", "AAV2_WT", "AAV9_wt", "AAV2_wt")] <- c("#800033ff", "#56b4e9ff", "#d55e00ff", "#009e73ff", "#d55e00ff", "#009e73ff")

# plot the correlation of the two libraries
p <- shared %>%
  ggplot(aes(x = mean_fraction_lib1, y = mean_fraction_lib3, color = AAV_ID)) +
  geom_point(color = "black", shape = 1, alpha = .5) +
  geom_smooth(method = "lm", color = "gray30", se=F, fullrange=T, size = .51) +
  geom_point(data = candidates, size = 2) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), breaks = 0:10/20, limits = c(0, .3)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = 0:10/20, limits = c(0, .3)) +
  scale_color_manual(values = aav_colors) + 
  coord_fixed() +
  theme(panel.grid.minor = element_blank()) +
  ggtitle("barcode proportion\n(mean of all samples)")

p

# plot the correlation again, but this time on log-space
p_log <- shared %>%
  mutate(mean_fraction_lib1 = if_else(mean_fraction_lib1 < 1e-4, 1e-4, mean_fraction_lib1),
         mean_fraction_lib3 = if_else(mean_fraction_lib3 < 1e-4, 1e-4, mean_fraction_lib3)) %>%
  ggplot(aes(x = mean_fraction_lib1, y = mean_fraction_lib3, color = AAV_ID)) +
  geom_point(color = "black", shape = 1, alpha = .5) +
  geom_smooth(method = "lm", color = "gray30", se=F, fullrange=T, size = .5) +
  geom_point(data = candidates, size = 2) +
  scale_x_log10(labels = c("100%", "10%", "1%", "0.1%", "≤0.01%"), breaks = 10^-(0:4), limits = c(1e-04, 1)) +
  scale_y_log10(labels = c("100%", "10%", "1%", "0.1%", "≤0.01%"), breaks = 10^-(0:4), limits = c(1e-04, 1)) +
  scale_color_manual(values = aav_colors) + 
  coord_fixed() +
  theme(legend.position = "none") +
  ggtitle("barcode proportion\n(mean of all samples, log scale)")

p_log

# plot both plots side by side
p + p_log &
  labs(x = "library 1",
       y = "library 3",
       color = "")

# save this figure
ggsave("plots/Lib-agreement.png", dpi = 300,
       width = 20, height = 10, units = "cm")

# test spearman correlation
shared %>%
  {cor.test(.$mean_fraction_lib1, .$mean_fraction_lib3, alternative = "greater", method = "spearman", exact = F)}



# Adress a specific Reviewer comment:
# "For example, the fold increase of AAV9_A2 over AAV6_wt does not appear similar in the two libraries."

# highlight these two examples in our plot
p_log +
  ggrepel::geom_text_repel(
    data = shared %>% filter(AAV_ID %in% c("AAV6_WT", "AAV9_A2")),
    mapping = aes(x = mean_fraction_lib1, y = mean_fraction_lib3, label = AAV_ID),
    inherit.aes = F, nudge_x = -1
  ) +
  labs(x = "library #1",
       y = "library #3",
       color = "")
# they are pretty close to the diagonal, so the results of both libraries are roughly similar I'd say.

# check the proportions in percent:
shared %>% filter(AAV_ID %in% c("AAV6_WT", "AAV9_A2")) %>% mutate_at(2:3, function(x) {x*100})


# fold change of AAV9_A2 over AAV6_WT in library #1:
25.3 / 1.14   #  = 22.19
23.0 / 0.567  #  = 40.56
# these fold changes are off by a factor of two, but considering
# these are different mice and different AAV libraries I'd say this
# is tolerable.


sessionInfo()
"
R version 4.0.2 (2020-06-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.7 LTS

Matrix products: default
BLAS:   /usr/lib/openblas-base/libblas.so.3
LAPACK: /usr/lib/libopenblasp-r0.2.18.so

locale:
[1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_GB.UTF-8   
[6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] patchwork_1.0.1 forcats_0.5.0   stringr_1.4.0   dplyr_1.0.1     purrr_0.3.4     readr_1.3.1     tidyr_1.1.1     tibble_3.0.3   
[9] ggplot2_3.3.2   tidyverse_1.3.0

loaded via a namespace (and not attached):
[1] tidyselect_1.1.0 splines_4.0.2    haven_2.3.1      lattice_0.20-41  colorspace_1.4-1 vctrs_0.3.2      generics_0.0.2   yaml_2.2.1      
[9] mgcv_1.8-31      utf8_1.1.4       blob_1.2.1       rlang_0.4.7      pillar_1.4.6     glue_1.4.1       withr_2.2.0      DBI_1.1.0       
[17] dbplyr_1.4.4     modelr_0.1.8     readxl_1.3.1     lifecycle_0.2.0  munsell_0.5.0    gtable_0.3.0     cellranger_1.1.0 rvest_0.3.6     
[25] fansi_0.4.1      broom_0.7.0      Rcpp_1.0.5       scales_1.1.1     backports_1.1.8  jsonlite_1.7.0   farver_2.0.3     fs_1.5.0        
[33] hms_0.5.3        digest_0.6.25    stringi_1.4.6    ggrepel_0.8.2    grid_4.0.2       cli_2.0.2        tools_4.0.2      magrittr_1.5    
[41] crayon_1.3.4     pkgconfig_2.0.3  ellipsis_0.3.1   Matrix_1.2-18    xml2_1.3.2       reprex_0.3.0     lubridate_1.7.9  assertthat_0.2.1
[49] httr_1.4.2       rstudioapi_0.11  R6_2.4.1         nlme_3.1-148     compiler_4.0.2
"
