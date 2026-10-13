library(tidyverse)
library(magrittr)
library(KRIS)
library(ggtext)
library(cowplot)

selfed_id <-
  read_table("metadata/s132.fam.txt", col_names = c("FamID", "IndID"), show_col_types = FALSE) %>%
  filter(FamID == "NZ") %>%
  anti_join(read_table("metadata/s132o.fam.txt", col_names = c("FamID", "IndID"), show_col_types = FALSE), by = "IndID") %>%
  pull(IndID)

read_set <-
  function(species, impact){
    pfx <- paste0("plink/", species, "_", impact)
    plink_bed <- read.bed(paste0(pfx, ".bed"), paste0(pfx, ".bim"), paste0(pfx, ".fam"), only.snp = FALSE)
    keep <- !as.character(plink_bed$ind.info$IndID) %in% selfed_id
    plink_bed$snp <- plink_bed$snp[keep, , drop = FALSE]
    plink_bed$ind.info <- plink_bed$ind.info[keep, ]
    plink_bed
  }

ind_load <-
  function(plink_bed){
    X <- 2 - plink_bed$snp
    X[!X %in% 0:2] <- NA
    flip <- colMeans(X, na.rm = TRUE) > 1
    X[, flip] <- 2 - X[, flip]

    tibble(FamID = as.character(plink_bed$ind.info$FamID),
           IndID = as.character(plink_bed$ind.info$IndID),
           masked = rowSums(X == 1, na.rm = TRUE),
           realized = rowSums(X == 2, na.rm = TRUE))
  }

pop_tb <-
  tribble(
    ~FamID_raw,      ~Population,
    "CEu-Belgium",   "Central Europe\n(Belgium)",
    "NAm-MA",        "North America\n(MA)",
    "SAm-Bariloche", "South America\n(Bariloche)",
    "AU-WA",         "Australia\n(WA)",
    "AU-other",      "Australia\n(other)",
    "NZ",            "New Zealand",
    "Aphal-Drake",   "Point Reyes"
  )

load_tb <-
  expand_grid(species = c("suilu", "aphal"), impact = c("HIGH", "LOW")) %>%
  pmap_dfr(function(species, impact){
    read_set(species, impact) %>%
      ind_load() %>%
      mutate(species = species, impact = impact)
  }) %>%
  mutate(FamID_raw = FamID) %>%
  left_join(pop_tb, by = "FamID_raw") %>%
  select(-FamID_raw) %>%
  mutate(Population = factor(Population, levels = pop_tb$Population),
         species = factor(case_match(species, "suilu" ~ "_S. luteus_", "aphal" ~ "_A. phalloides_"),
                          levels = c("_S. luteus_", "_A. phalloides_")),
         impact = factor(impact, levels = c("HIGH", "LOW")))

pop_theme <-
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 18),
        plot.margin = margin(5, 10, 5, 30),
        panel.spacing = unit(1.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        strip.text = element_markdown(face = "bold", size = 16))

aphal_k <- 2.3
species_facet <-
  list(facet_grid(. ~ species, scales = "free_x"),
       ggh4x::force_panelsizes(cols = c(6, aphal_k)))

by_species <-
  function(layer_fun){
    list(layer_fun(\(d) filter(d, species == "_S. luteus_"), 1),
         layer_fun(\(d) filter(d, species == "_A. phalloides_"), 1 / aphal_k))
  }

# A: masked HIGH load
A <-
  load_tb %>%
  filter(impact == "HIGH") %>%
  ggplot(aes(x = Population, y = masked)) +
  by_species(\(sub, k) geom_boxplot(data = sub, outlier.shape = NA, width = 0.75 * k)) +
  by_species(\(sub, k) geom_jitter(data = sub, width = 0.15 * k, size = 0.8, shape = 1)) +
  ylim(0, NA) +
  species_facet +
  labs(x = NULL, y = "Masked HIGH load\n(heterozygous HIGH variants)") +
  pop_theme

# B: realized / masked, HIGH vs LOW
sig_tb <-
  load_tb %>%
  mutate(value = realized / masked) %>%
  select(species, Population, IndID, impact, value) %>%
  pivot_wider(names_from = impact, values_from = value) %>%
  group_by(species, Population) %>%
  summarise(p = binom.test(sum(HIGH < LOW), sum(HIGH != LOW), alternative = "greater")$p.value,
            y = max(HIGH, LOW) + 0.02, .groups = "drop") %>%
  group_by(species) %>%
  mutate(label = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"),
         x = rank(as.numeric(Population)),
         dx = if_else(species == "_A. phalloides_", 0.175 / aphal_k, 0.175)) %>%
  ungroup()

B <-
  load_tb %>%
  mutate(value = realized / masked) %>%
  ggplot(aes(x = Population, y = value, fill = impact)) +
  by_species(\(sub, k) geom_boxplot(data = sub, outlier.shape = NA, position = position_dodge(0.7 * k), width = 0.7 * k)) +
  by_species(\(sub, k) geom_point(data = sub, position = position_jitterdodge(jitter.width = 0.15 * k, dodge.width = 0.7 * k), size = 0.7, shape = 1)) +
  geom_segment(data = sig_tb, aes(x = x - dx, xend = x + dx, y = y, yend = y), inherit.aes = FALSE) +
  geom_segment(data = sig_tb, aes(x = x - dx, xend = x - dx, y = y, yend = y - 0.008), inherit.aes = FALSE) +
  geom_segment(data = sig_tb, aes(x = x + dx, xend = x + dx, y = y, yend = y - 0.008), inherit.aes = FALSE) +
  geom_text(data = sig_tb, aes(x = x, y = y + 0.004, label = label), inherit.aes = FALSE, vjust = 0) +
  scale_fill_manual(values = c(HIGH = "grey60", LOW = "white"), name = "Impact") +
  ylim(0, NA) +
  species_facet +
  labs(x = NULL, y = "Realized Load/Masked Load") +
  pop_theme

plot_grid(A, B, nrow = 1, labels = "AUTO", align = "h", rel_widths = c(1, 1.15))
ggsave("plots/load.pdf", width = 14*1.2, height = 5*1.2, units = "in")
ggsave("plots/load.png", width = 14*1.2, height = 5*1.2, units = "in", dpi = 200)
