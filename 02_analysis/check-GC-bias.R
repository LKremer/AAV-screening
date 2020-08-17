library(tidyverse)
library(patchwork)

# define the plotting style
theme_set(theme_bw() + theme(text = element_text(family = 'Arial'),
                             panel.grid.minor = element_blank(),    # remove axis grid
                             panel.grid.major = element_blank(),
                             strip.background = element_blank(),    # remove silly gray bg of facet title
                             strip.text = element_text(hjust = 0),  # left-justify facets
                             plot.tag = element_text(size = 20, face = "bold")
))
aav_colors <- c("other" = "gray50")
aav_colors[c("AAV1_P5", "AAV9_A2", "AAV9_WT", "AAV2_WT")] <- c("#800033ff", "#56b4e9ff", "#d55e00ff", "#009e73ff")

# Calculate GC-content and mean barcode proportion per AAV (library #1)
lib1 <- read_csv("tables/AAVlib1_AAV-barcode-counts.csv", col_types = "cccdcddddddd") %>%
  group_by(set, celltype) %>%
  mutate(total_bc = sum(barcode_count)) %>% ungroup() %>%
  group_by(set, celltype, AAV_ID, barcode_seq) %>%
  summarise(fraction = barcode_count / total_bc) %>% ungroup() %>%
  group_by(AAV_ID, barcode_seq) %>%
  summarise(mean_fraction = mean(fraction)) %>% ungroup() %>%
  mutate(gc_content = str_count(barcode_seq, "[GC]") / 15)

# Calculate GC-content and mean barcode proportion per AAV (library #3)
lib3 <- read_csv("tables/AAVlib3_AAV-barcode-counts.csv", col_types = "cccdcddddddd") %>%
  mutate(AAV_ID = str_replace(AAV_ID, "wt", "_WT")) %>%
  group_by(set, celltype) %>%
  mutate(total_bc = sum(barcode_count)) %>% ungroup() %>%
  group_by(set, celltype, AAV_ID, barcode_seq) %>%
  summarise(fraction = barcode_count / total_bc) %>% ungroup() %>%
  group_by(AAV_ID, barcode_seq) %>%
  summarise(mean_fraction = mean(fraction)) %>% ungroup() %>%
  mutate(gc_content = str_count(barcode_seq, "[GC]") / 15)

# plot GC-content against barcode proportion for both libraries, highlighting our candidates:
lib1_candidates <- lib1 %>% filter(AAV_ID %in% names(aav_colors))
p1 <- lib1 %>%
  ggplot(aes(x = gc_content, y = mean_fraction, color = AAV_ID, label = AAV_ID)) +
  geom_smooth(method = "lm", color = "gray70", se=F, fullrange=T, size = .5) +
  geom_point(color = "black", shape = 1, alpha = .5) +
  geom_point(data = lib1_candidates, size = 2) +
  ggrepel::geom_text_repel(data = lib1_candidates) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_log10(labels = c("100%", "10%", "1%", "0.1%", "0.01%", "1e-3%", "1e-4%", "1e-5%", "1e-6%"), breaks = 10^-(0:8), limits = c(NA, 1)) +
  scale_color_manual(values = aav_colors) +
  ggtitle("library #1")

lib3_candidates <- lib3 %>% filter(AAV_ID %in% names(aav_colors))
p3 <- lib3 %>%
  ggplot(aes(x = gc_content, y = mean_fraction, color = AAV_ID, label = AAV_ID)) +
  geom_smooth(method = "lm", color = "gray70", se=F, fullrange=T, size = .5) +
  geom_point(color = "black", shape = 1, alpha = .5) +
  geom_point(data = lib3_candidates, size = 2) +
  ggrepel::geom_text_repel(data = lib3_candidates) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_log10(labels = c("100%", "10%", "1%", "0.1%", "0.01%", "1e-3%", "1e-4%", "1e-5%", "1e-6%"), breaks = 10^-(0:8), limits = c(NA, 1)) +
  scale_color_manual(values = aav_colors) +
  ggtitle("library #3")

p1 + p3 &
  theme(legend.position = "none", panel.grid.minor = element_blank()) &
  labs(x = "barcode GC-content", y = "barcode proportion")

# save plots
ggsave("plots/GC-bias.png", dpi = 300,
       width = 18, height = 9, units = "cm")
ggsave("plots/GC-bias.svg",
       width = 18, height = 9, units = "cm")

# Calculate spearman correlation for both libraries:
lib1 %>%
  {cor.test(.$gc_content, .$mean_fraction, method = "spearman", exact = F)}

lib3 %>%
  {cor.test(.$gc_content, .$mean_fraction, method = "spearman", exact = F)}
