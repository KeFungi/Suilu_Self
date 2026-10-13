library(tidyverse)
library(pheatmap)
library(magrittr)
library(plotly)
library(cowplot)

slength <- 40027171 #non-repeat length of genome

s189_list <-
  read_table("plink/s189.fam", col_names = FALSE) %>%
  select(X1, X2)

s132_list <-
  read_table("metadata/s132.fam.txt", col_names = FALSE) %>%
  filter(X1!="NZ")

#introduced_tb <-
#tribble(~ "FID", ~"introduced", ~"FullName",
#        "CEu-Belgium", "Native", "Central Europe\n(Belgium)",
#        "AU-other", "Introduced", "Australia\n(other)",
#        "AU-WA", "Introduced", "Australia\n(WA)",
#        "NZ", "Introduced", "New Zealand",
#        "SAm-Bariloche", "Introduced", "South America\n(Bariloche)",
#        "SAm", "Introduced", "South America\n(all)",
#        "NAm-MA", "Introduced", "North America\n(MA)",
#        "NAm", "Introduced", "North America\n(all)",
#        "CEu", "Native", "Central Europe\n(all)",
#        "AU", "Introduced", "Australia\n(all)",
#        "AF", "Introduced", "Africa"
#)

introduced_tb <-
  tribble(~ "FID", ~"introduced", ~"FullName",
          "NZ", "Introduced", "New Zealand",
          "SAm", "Introduced", "South America",
          "NAm", "Introduced", "North America",
          "CEu", "Native", "Central Europe",
          "AU", "Introduced", "Australia",
          "AF", "Introduced", "Africa"
  ) %>%
  mutate(introduced=factor(introduced, levels=c("Native", "Introduced")))

meta_tb <-
  rbind(s189_list) %>%
  rename(FID=X1, IID=X2) %>%
  right_join(introduced_tb)

read_het <-
  function(hetfile){
    read_table(hetfile) %>%
      mutate(`H`= (`N(NM)`-`O(HOM)`)/`N(NM)`)
  }

sum_roh <-
  function(file, n) {
    read_lines(file) %>%
      grep("^RG",., value = TRUE) %>%
      paste0(collapse = "\n") %>%
      read_tsv(col_names = c("RG", "Sample", "Chromosome", "Start", "End", "Length", "Number of markers", "Quality")) %>%
      filter(Length > n) %>%
      pull(Length) %>%
      sum()
  }

unfilterd_het <-
  read_het("plink/s189.het")

het_plot <-
  unfilterd_het %>%
    select(-FID) %>%
    left_join(meta_tb) %>%
  ggplot() +
  geom_jitter(aes(x=FullName, y= H, color=FullName, group=IID), alpha=0.5, height = 0, width = 0.1) +
  geom_boxplot(aes(x=FullName, y=H), fill=NA, outlier.shape = NA) +
  facet_grid(cols=vars(introduced), scales = "free_x", space="free_x") +
  theme(axis.text.x=element_blank()) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("heterozygosity") +
  ylim(0, NA)

f_plot <-
  unfilterd_het %>%
    select(-FID) %>%
    left_join(meta_tb) %>%
  ggplot() +
  geom_jitter(aes(x=FullName, y= `F`, color=FullName, group=IID), alpha=0.5, height = 0, width = 0.1) +
  geom_boxplot(aes(x=FullName, y=`F`), fill=NA, outlier.shape = NA) +
  facet_grid(cols=vars(introduced), scales = "free_x", space="free_x") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none", ) +
  xlab("") +
  ylab("F_hom")

roh_tb <-
  list.files("bcfroh", pattern = paste0(".+\\.roh_RG\\.txt"), full.names=TRUE) %>%
  tibble(rohfile=.) %>%
  mutate(FID=str_match(rohfile, "^bcfroh/([^_]+)_.+\\.roh_RG\\.txt$")[,2],
         IID=str_match(rohfile, "^bcfroh/[^_]+_(.+)\\.roh_RG\\.txt$")[,2]) %>%
  select(-FID)

roh10000_tb <-
  roh_tb %>%
  left_join(meta_tb) %>%
  mutate(roh_len = map_dbl(rohfile, ~sum_roh(.x, 10000)))

roh_plot <-
  ggplot(roh10000_tb) +
  geom_jitter(aes(x=FullName, y=roh_len/slength, color=FullName, group=IID),alpha=0.5, height = 0, width = 0.1) +
  geom_boxplot(aes(x=FullName, y=roh_len/slength), fill=NA, outlier.shape = NA) +
  facet_grid(cols=vars(introduced), scales = "free_x", space="free_x") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none", ) +
  ylab("F_roh") +
  xlab("")

final_plot <-
  plot_grid(het_plot, f_plot, roh_plot, ncol=3, labels = "AUTO")

ggsave("plots/het_plot.pdf", height = 3, width=10)

ggplotly(het_plot)
ggplotly(f_plot)
ggplotly(roh_plot)
