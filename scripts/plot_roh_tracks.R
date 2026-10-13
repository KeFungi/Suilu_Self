library(tidyverse)
library(magrittr)
library(ggtext)

fai <-
  read_tsv("ref/Suilu4.fasta.fai",
           col_names = c("scaffold", "length", "off1", "off2", "off3"),
           show_col_types = FALSE) %>%
  arrange(desc(length)) %>%
  mutate(scaffold = factor(scaffold, levels = scaffold),
         cum_end   = cumsum(length),
         cum_start = lag(cum_end, default = 0),
         cum_mid   = (cum_start + cum_end) / 2,
         band      = as.integer(scaffold) %% 2)

top_scaffolds <- as.character(fai$scaffold)

roh <-
  list.files("bcfroh", pattern = "\\.roh_RG\\.txt$", full.names = TRUE) %>%
  map_dfr(~ read_tsv(.x,
             col_names = c("rg", "sample", "scaffold", "start", "end",
                           "length", "n_markers", "quality"),
             show_col_types = FALSE)) %>%
  filter(scaffold %in% top_scaffolds) %>%
  mutate(pop      = str_extract(sample, "^[^_]+"),
         scaffold = factor(scaffold, levels = top_scaffolds))

pop_tb <-
  tribble(~ "pop", ~ "FullName",
          "CEu",   "Central Europe (native)",
          "SAm",   "South America",
          "NAm",   "North America",
          "AU",    "Australia",
          "NZ",    "New Zealand",
          "AF",    "Africa (n=1)") %>%
  mutate(FullName = factor(FullName, levels = FullName))

roh %<>% left_join(pop_tb, by = "pop")

sample_order <-
  roh %>%
  group_by(FullName, sample) %>%
  summarise(total_roh = sum(length), .groups = "drop") %>%
  arrange(FullName, desc(total_roh)) %>%
  pull(sample)

roh %<>% mutate(sample = factor(sample, levels = sample_order))

fai_cum <- fai %>% filter(scaffold %in% top_scaffolds)
roh_cum <-
  roh %>%
  left_join(fai_cum %>% select(scaffold, cum_start), by = "scaffold") %>%
  mutate(gx_start = (start + cum_start) / 1e6,
         gx_end   = (end   + cum_start) / 1e6)

fai_cum %<>%
  mutate(cum_start_Mb = cum_start / 1e6,
         cum_end_Mb   = cum_end   / 1e6,
         cum_mid_Mb   = cum_mid   / 1e6)

roh_track <-
  ggplot() +
  geom_rect(data = fai_cum,
            aes(xmin = cum_start_Mb, xmax = cum_end_Mb,
                ymin = -Inf, ymax = Inf,
                fill = factor(band)),
            alpha = 0.35, inherit.aes = FALSE) +
  geom_segment(data = roh_cum,
               aes(x = gx_start, xend = gx_end,
                   y = sample, yend = sample,
                   color = FullName),
               linewidth = 0.35) +
  geom_vline(data = fai_cum,
             aes(xintercept = cum_end_Mb),
             color = "grey40", linewidth = 0.15) +
  scale_x_continuous(name = "Genome position (Mb, all 67 scaffolds concatenated)",
                     breaks = fai_cum$cum_mid_Mb,
                     labels = as.character(fai_cum$scaffold),
                     expand = c(0, 0)) +
  scale_fill_manual(values = c("0" = "grey95", "1" = "white"), guide = "none") +
  facet_grid(rows = vars(FullName), scales = "free_y", space = "free_y",
             switch = "y") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.spacing.y = unit(0.1, "lines"),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1),
        strip.background.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(5, 10, 5, 30))

ggsave("plots/roh_tracks.pdf", roh_track, width = 16, height = 12)
ggsave("plots/roh_tracks.png", roh_track, width = 16, height = 12, dpi = 200)

pop_summary <-
  roh %>%
  group_by(FullName) %>%
  summarise(n_samples       = n_distinct(sample),
            median_total_kb = round(median(tapply(length, sample, sum)) / 1e3, 1),
            median_n_roh    = round(median(tapply(length, sample, length)), 0),
            .groups = "drop") %>%
  arrange(FullName)

write_csv(pop_summary, "plots/roh_tracks_summary.csv")
