library(tidyverse)
library(patchwork)

theme_set(theme_minimal())

shared_aavs <- intersect(
  unique(pull(read_csv("tables/AAVlib1_AAV-barcode-counts.csv", col_types = "----c-------"))),
  unique(pull(read_csv("tables/AAVlib3_AAV-barcode-counts.csv", col_types = "----c-------") %>% mutate(AAV_ID = str_replace(AAV_ID, "wt", "_WT")))))

lib1 <- read_csv("tables/AAVlib1_AAV-barcode-counts.csv", col_types = "cccdcd------") %>%
  filter(AAV_ID %in% shared_aavs) %>%
  group_by(set, celltype) %>%
  mutate(total_bc = sum(barcode_count)) %>% ungroup() %>%
  group_by(set, celltype, AAV_ID) %>%
  summarise(fraction = barcode_count / total_bc) %>% ungroup() %>%
  group_by(AAV_ID) %>%
  summarise(mean_fraction = mean(fraction)) %>% ungroup()

lib3 <- read_csv("tables/AAVlib3_AAV-barcode-counts.csv", col_types = "cccdcd------") %>%
  mutate(AAV_ID = str_replace(AAV_ID, "wt", "_WT")) %>%
  filter(AAV_ID %in% shared_aavs) %>%
  filter(!str_detect(set, "vitro")) %>%
  group_by(set, celltype) %>%
  mutate(total_bc = sum(barcode_count)) %>% ungroup() %>%
  group_by(set, celltype, AAV_ID) %>%
  summarise(fraction = barcode_count / total_bc) %>% ungroup() %>%
  group_by(AAV_ID) %>%
  summarise(mean_fraction = mean(fraction)) %>% ungroup()

shared <- full_join(lib1, lib3, by = "AAV_ID", suffix = c("_lib1", "_lib3"))
candidates <- shared %>%
  filter(AAV_ID %in% c("AAV1_P5", "AAV9_A2", "AAV9_WT", "AAV2_WT"))

aav_colors <- c("other" = "gray50")
aav_colors[c("AAV1_P5", "AAV9_A2", "AAV9_WT", "AAV2_WT")] <- c("#800033ff", "#56b4e9ff", "#d55e00ff", "#009e73ff")

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
  ggtitle("barcode / library\n(mean of all samples)")

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
  ggtitle("barcode / library\n(mean of all samples, log scale)")

p + p_log &
  labs(x = "library 1",
       y = "library 3",
       color = "")

ggsave("plots/Lib-agreement.png", dpi = 300,
       width = 20, height = 10, units = "cm")

shared %>%
  #mutate_at(vars(contains("fraction")), log1p) %>%
  #filter(mean_fraction_lib1 > 0 & mean_fraction_lib3 > 0) %>%
  {cor.test(.$mean_fraction_lib1, .$mean_fraction_lib3, alternative = "greater", method = "spearman", exact = F)}

shared %>%
  {cor.test(.$mean_fraction_lib1, .$mean_fraction_lib3, alternative = "greater", method = "pearson")}
