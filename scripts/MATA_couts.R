library(tidyverse)
library(cowplot)
library(ggVennDiagram)
library(magrittr)
source("scripts/fun_MAT_allele_ML.R")

s189_list <-
  read_table("plink/s189.fam", col_names = FALSE) %>%
  select(FID=X1, IID=X2)

s132_list <-
  read_table("metadata/s132.fam.txt", col_names = FALSE) %>%
  filter(X1!="NZ")

introduce_tb <-
  tribble(~ "FID", ~"introduced", ~"FullName",
          "NZ", "Introduced", "New Zealand",
          "SAm", "Introduced", "South America",
          "NAm", "Introduced", "North America",
          "CEu", "Native", "Central Europe",
          "AU", "Introduced", "Australia",
          "AF", "Introduced", "Africa"
  ) %>%
mutate(introduced=factor(introduced, levels=c("Native", "Introduced")))

MATA_types <-
  read_csv("results/MATA_groups.csv") %>%
  mutate(sequence=str_replace(sequence, "_[^_]+$", "")) %>%
  select(IID=sequence, A_type=HD1_group) %>%
  right_join(s189_list) %>%
  left_join(introduce_tb)

di_inds <-
  MATA_types %>%
  group_by(IID) %>%
  summarise(n=n()) %>%
  filter(n==2) %>%
  pull(IID)

di_MATA_types <-
  MATA_types[MATA_types$IID %in% di_inds, ]

duplicat_MATA_cross_region <-
  di_MATA_types %>%
  group_by(IID, FID) %>%
  summarise(types=list(c(A_type))) %>%
  rowwise() %>%
  mutate(types=list(sort(types))) %>%
  group_by(types) %>%
  summarise(IIDs=list(c(IID)), regions=list(c(FID)), n=n()) %>%
  arrange(desc(n)) %>%
  rowwise() %>%
  mutate(regions=list(unique(regions))) %>%
  filter(length(regions) >=2) %>%
  mutate(types=paste(types, collapse = ", "),
         IIDs=paste(IIDs, collapse = ", "),
         regions=paste(regions, collapse = ", "))

duplicat_MATA_long_region <-
  di_MATA_types %>%
  group_by(IID, FID) %>%
  summarise(types=list(c(A_type))) %>%
  rowwise() %>%
  mutate(types=list(sort(types))) %>%
  group_by(types, FID) %>%
  summarise(IIDs=list(c(IID)), n=n()) %>%
  arrange(desc(n)) %>%
  rowwise() %>%
  mutate(types=paste(types, collapse = ", "),
         IIDs=paste(IIDs, collapse = ", "))

duplicate_tb <-
  duplicat_MATA_long_region %>%
  group_by(FID, n) %>%
  summarise(freq=n()) %>%
  pivot_wider(names_from = "n", values_from = "freq", values_fill = 0) %>%
  ungroup() %>%
  left_join(introduce_tb) %>%
  select(FullName, `Once`=`1`, `Twice`=`2`, `Thrice`=`3`, `Quarce`=`4`)

write_csv(duplicate_tb, "results/MATA_genotypeFreq.csv")

allele_counts <-
  di_MATA_types %>%
  group_by(A_type, FID) %>%
  summarise(n=n()) %>%
  pivot_wider(values_from = "n", names_from = "A_type", values_fill = 0)

allele_list <-
  di_MATA_types %>%
  group_by(A_type, FID) %>%
  summarise(n=n()) %>%
  group_by(FID) %>%
  summarise(allele_type = list(c(A_type))) %>%
  left_join(introduce_tb)

private_alleles <-
  allele_counts %>%
  column_to_rownames("FID") %>%
  mutate_all(~.x>0) %>%
  .[,colSums(.) == 1] %>%
  rownames_to_column("FID") %>%
  pivot_longer(-FID) %>%
  filter(value) %>%
  mutate(label="private") %>%
  arrange(FID)

MATA_stack_data <-
  di_MATA_types %>%
  group_by(A_type, FID) %>%
  summarise(n=n()) %>%
  mutate(A_type=as.character(A_type)) %>%
  left_join(private_alleles, by=c("A_type"="name", "FID"="FID")) %>%
  mutate(label=ifelse(is.na(label), A_type, label))

allele_dis_tb <-
  allele_counts %>%
  pmap_dfr(function(...){
    row <- tibble(...)
    allele_dis <-
      row[,-1] %>%
      reduce(cbind) %>%
      .[which(. != 0)]
    tibble(FID=row$FID, allele_dis=list(allele_dis))
  }
  ) %>%
  rowwise()

allele_k <-
  allele_dis_tb %>%
  mutate(max_dup=max(allele_dis)) %>%
  filter(max_dup>2) %>%
  mutate(k_dup=cal_k_dup(allele_dis),
         k_dup_number=cal_k_dup(allele_dis, FALSE),
         k_ML=cal_MLk_range(allele_dis, 1:100)
  ) %>%
  left_join(introduce_tb) %>%
  mutate(k=ifelse(k_ML==100,
                  paste0(round(k_dup, 1), "(", k_ML, ")*"),
                  paste0(round(k_dup, 1), "(", k_ML, ")")
                  )
  )

allele_k %>%
  mutate(n=sum(allele_dis)) %>%
  select(FullName, k_ML, k_dup_number, n) %>%
  unnest(k_dup_number) %>%
  select(Population=FullName,
         k=n,
         `incompatible crosses`=d,
         `total crosses`=pair,
         n_incom=k,
         n_MLE=k_ML) %>%
  mutate(n_incom=round(n_incom, 1)) %>%
  write_csv("results/MAT_n.csv")

allele_ML <-
  allele_dis_tb %>%
  mutate(max_dup=max(allele_dis)) %>%
  filter(max_dup>2) %>%
  mutate(k_ML=list(cal_Lk_range(allele_dis, 1:100))) %>%
  unnest(k_ML)

allele_L_plot <-
  allele_ML %>%
  left_join(introduce_tb) %>%
  ggplot() +
  geom_line(aes(x=k, y=L)) +
  theme_light() +
  facet_wrap(~FullName, scales = "free") +
  xlab("Total Allele Types in Population (n)") +
  ylab("Likelihood")

ggsave("plots/MAT_ML.pdf", allele_L_plot, width = 6, height =6)

MATA_stack_sumdata <-
  MATA_stack_data %>%
  group_by(FID) %>%
  summarise(total=n(),
            private=sum(label=="private"),
            shared=sum(label!="private")) %>%
  mutate(across(-1, ~ifelse(.x==0, NA, .x))) %>%
  left_join(introduce_tb)

#MATA_type_plot <-
#  MATA_stack_data %>%
#  left_join(introduce_tb) %>%
#  arrange(introduced) %>%
#  ggplot() +
#  geom_bar(aes(x=FullName, fill=label), stat="count", position = "stack") +
#  geom_text(data=MATA_stack_sumdata, aes(x=FullName, y= shared+1, label = shared), vjust = 1.5, colour = "black") +
#  geom_text(data=MATA_stack_sumdata, aes(x=FullName, y= 1, label = private), colour = "black") +
#  facet_grid(col=vars(introduced), scales = "free_x", space = "free_x") +
#  scale_fill_manual(values = c("private"="grey", "shared"="grey50")) +
#  theme_classic() +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1 ),
#        axis.title.x = element_blank(),
#        legend.title = element_blank()) +
#  ylab("Sampled Allele Types")

#legend1 <- get_legend(MATA_type_plot)

MATA_stack_sumdata <-
  MATA_stack_data %>%
  group_by(FID) %>%
  summarise(total=sum(n)) %>%
  left_join(introduce_tb)

MATA_count_plot <-
  MATA_stack_data %>%
  left_join(introduce_tb) %>%
  arrange(introduced) %>%
  group_by(FullName) %>%
  mutate(sum_n = sum(n)) %>%
  ungroup() %>%
  ggplot() +
  geom_bar(aes(x=FullName, fill=A_type, y=n), stat="identity", color="black", position = "stack") +
  #geom_text(data=MATA_stack_sumdata, aes(x=FullName, y=total, label=total), vjust=-0.4) +
  facet_grid(col=vars(introduced), scales = "free_x", space = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1 ),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Allele Counts") +
  ylim(0, 110)

MATA_n_data <-
  allele_k %>%
  left_join(introduce_tb) %>%
  arrange(introduced) %>%
  rename(`incompatible pairs`=k_dup,
         `multinomial distribution`=k_ML
  )

MATA_n_plot <-
  MATA_n_data %>%
  pivot_longer(c(`incompatible pairs`, `multinomial distribution`),
               names_to = "Method") %>%
  mutate(Method=case_match(Method, "incompatible pairs" ~ "incompatible\npairs\n", "multinomial distribution" ~ "multinomial\ndistribution\n")) %>%
  ggplot() +
  geom_col(aes(x=FullName, y=value, fill=Method), position = "dodge") +
  geom_text(data=MATA_n_data, aes(x=FullName, y=`multinomial distribution`+2, label = k), colour = "black", size=3) +
  facet_grid(col=vars(introduced), scales = "free_x", space = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1 ),
        axis.title.x = element_blank()) +
  ylab("Estimated total mating types")

legend2 <- get_legend(MATA_n_plot)

#final_plot <-
#  plot_grid(MATA_type_plot + theme(legend.position = "none"), legend1,
#            MATA_count_plot, ggplot()+theme_nothing(),
#            ncol = 2, nrow=2, rel_widths = c(1, 0.2),
#            labels=c("A", "", "B", "")
#  )

plot_grid(MATA_count_plot, ggplot()+theme_nothing(), MATA_n_plot + theme(legend.position = "none"), legend2,
          labels=c("A", "", "B", ""), ncol=2, rel_widths = c(1, 0.3))

ggsave("plots/MATA_count.pdf", width = 6, height = 8)
names(allele_list$allele_type) <- allele_list$FullName

# ---- Master allele-set list ---------------------------------
# A single named list `MATA_sets` containing the unique A_type values
# for every category found in BOTH s189 (broad) and s132 (sub-region)
# groupings. One library to look up anything by name.
#
#   names(MATA_sets)                  # see every available group
#   MATA_sets[["Australia"]]          # one set
#   MATA_sets[c("AU-WA", "NZ")]       # multi-set lookup for Venn
#
# Categories included:
#   s189 FID       : AU, NZ, SAm, NAm, CEu, AF
#   s189 FullName  : Australia, New Zealand, South America,
#                    North America, Central Europe, Africa
#   s132 sub-FID   : AU-WA, AU-other, NAm-MA, SAm-Bariloche,
#                    CEu-Belgium  (NZ re-attached from s189)

# Long form: one row per (IID, A_type, group_label) for every grouping.
MATA_membership <-
  bind_rows(
    # s189 by raw FID
    di_MATA_types %>%
      filter(!is.na(FID)) %>%
      transmute(IID, A_type, group = FID),
    # s189 by FullName
    di_MATA_types %>%
      filter(!is.na(FullName)) %>%
      transmute(IID, A_type, group = FullName),
    # s132 by sub-FID (X1 = sub-region, X2 = IID); NZ is filtered out of
    # s132_list at the top of the script, so add it back from s189
    di_MATA_types %>%
      inner_join(s132_list %>% select(IID = X2, group = X1), by = "IID") %>%
      select(IID, A_type, group),
    di_MATA_types %>%
      filter(FID == "NZ") %>%
      transmute(IID, A_type, group = "NZ")
  ) %>%
  distinct()

MATA_sets <-
  MATA_membership %>%
  group_by(group) %>%
  summarise(allele_type = list(sort(unique(A_type))), .groups = "drop") %>%
  { setNames(.$allele_type, .$group) }

# ---- Venn diagrams of MATA allele sharing -------------------
# All plots are simple lookups into MATA_sets.
venn_intro4 <-
  MATA_sets[c("South America", "North America", "Australia", "New Zealand")] %>%
  ggVennDiagram(label = "count", set_size = 4.5) +
  scale_fill_gradient(low = "white", high = "red", name = "count") +
  scale_color_manual(values = rep("black", 4)) +
  coord_equal(clip = "off") +
  theme(plot.margin = margin(10, 20, 10, 20))

ggsave("plots/MATA_venn.pdf", venn_intro4, height = 6, width = 7)

