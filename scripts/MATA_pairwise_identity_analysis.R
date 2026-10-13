library(Biostrings)
library(tidyverse)
library(cowplot)

fasta_file      <- "MATA/HD1_s189_aa alignment 200aa.fasta"
identity_cutoff <- 0.90

read_alignment <- function(path) {
  seqs <- readBStringSet(path)
  if (length(seqs) == 0L) stop("No sequences found in: ", path)

  # Sniff alphabet from a sample of the sequences
  sample_n  <- min(length(seqs), 20L)
  sample_chars <- unique(unlist(strsplit(
    toupper(as.character(seqs[seq_len(sample_n)])), ""
  )))
  nuc_alphabet <- c("A", "C", "G", "T", "U", "N", "-", ".")
  nuc_frac     <- mean(sample_chars %in% nuc_alphabet)
  seq_type     <- if (nuc_frac >= 0.9) "nucleotide" else "amino acid"

  attr(seqs, "seq_type") <- seq_type
  seqs
}

aln      <- read_alignment(fasta_file)
seq_type <- attr(aln, "seq_type")

if (length(unique(width(aln))) != 1L) {
  stop("Sequences are not the same length - is the file actually aligned?")
}

n   <- length(aln)
ids <- names(aln)

# n x L character matrix (uppercase so 'a' and 'A' are treated equally
# for both nucleotide and amino acid alignments)
char_mat <- toupper(as.matrix(aln))

# non-gap indicator and pairwise count of shared non-gap columns.
# Treat both '-' and '.' as gap characters (common in AA alignments).
gap_chars   <- c("-", ".")
nongap      <- (!(char_mat %in% gap_chars)) + 0   # numeric n x L
dim(nongap) <- dim(char_mat)
both_nongap <- tcrossprod(nongap)                 # n x n

# Sum of co-occurrence matrices over each non-gap residue
residues <- setdiff(unique(as.vector(char_mat)), gap_chars)
matches  <- matrix(0, n, n)
for (a in residues) {
  ind         <- (char_mat == a) + 0
  dim(ind)    <- dim(char_mat)
  matches <- matches + tcrossprod(ind)
}

pid_mat <- matches / both_nongap
pid_mat[!is.finite(pid_mat)] <- NA   # pairs with zero shared columns
diag(pid_mat) <- 1
dimnames(pid_mat) <- list(ids, ids)

# Long-format table of unique pairs (upper triangle)
pid_tbl <- which(upper.tri(pid_mat), arr.ind = TRUE) %>%
  as_tibble() %>%
  rename(i = row, j = col) %>%
  mutate(
    seq1     = ids[i],
    seq2     = ids[j],
    identity = round(pid_mat[cbind(i, j)], 3)
    )

# Save matrix and long table
as_tibble(round(pid_mat, 3), rownames = "sequence") %>%
  write_csv("results/HD1_pairwise_identity_matrix.csv")

pid_tbl
# ---- 3. Distribution plots ---------------------------------
pid_df <- pid_tbl %>% select(identity)

# Top: histogram
HD1_plot <-
  ggplot(pid_df, aes(x = identity)) +
  geom_histogram(breaks = seq(0.55, 1, 0.025)) +
  geom_jitter(aes(x=identity, y=-1000), alpha=0.02, width = 0, size=1, height = 500) +
  xlab("pairwise identity") +
  ylab("count") +
  theme_bw()

# ---- 4. Group sequences by an identity cutoff --------------
# Single-linkage clustering: two sequences are placed in the
# same group if a chain of pairs each with identity >= cutoff
# connects them (i.e., connected components of the graph whose
# edges are pairs with identity >= cutoff).
#
# We get this by feeding (1 - identity) as a distance to
# hclust(method = "single") and cutting at height (1 - cutoff).
# `identity_cutoff` is set near the top of this script (default 0.9)
# and may be overridden via the second command-line argument.

dist_mat <- 1 - pid_mat
dist_mat[is.na(dist_mat)] <- 1     # treat undefined pairs as fully diverged
diag(dist_mat) <- 0

hc      <- hclust(as.dist(dist_mat), method = "single")
cluster <- cutree(hc, h = 1 - identity_cutoff)

groups_tbl <-
  tibble(sequence = ids, group = unname(cluster)) %>%
  arrange(sequence) %>%
  mutate(group = factor(group, levels=unique(group)),
         group = as.numeric(group)
  )

# Per-group summary
# Two-column padded text file in the same style as results/MATA_types.txt
# (sequence name left-aligned, group number after)
name_width <- max(nchar(groups_tbl$sequence)) + 1
groups_tbl %>%
  arrange(group, sequence) %>%
  mutate(line = sprintf(paste0("%-", name_width, "s %d"), sequence, group)) %>%
  pull(line) %>%
  writeLines("results/HD1_groups.txt")

#####

library(Biostrings)
library(tidyverse)
library(patchwork)

fasta_file      <- "MATA/HD2_s189_aa alignment 200aa.fasta"
identity_cutoff <- 0.90

read_alignment <- function(path) {
  seqs <- readBStringSet(path)
  if (length(seqs) == 0L) stop("No sequences found in: ", path)

  # Sniff alphabet from a sample of the sequences
  sample_n  <- min(length(seqs), 20L)
  sample_chars <- unique(unlist(strsplit(
    toupper(as.character(seqs[seq_len(sample_n)])), ""
  )))
  nuc_alphabet <- c("A", "C", "G", "T", "U", "N", "-", ".")
  nuc_frac     <- mean(sample_chars %in% nuc_alphabet)
  seq_type     <- if (nuc_frac >= 0.9) "nucleotide" else "amino acid"

  attr(seqs, "seq_type") <- seq_type
  seqs
}

aln      <- read_alignment(fasta_file)
seq_type <- attr(aln, "seq_type")

if (length(unique(width(aln))) != 1L) {
  stop("Sequences are not the same length - is the file actually aligned?")
}

n   <- length(aln)
ids <- names(aln)

# n x L character matrix (uppercase so 'a' and 'A' are treated equally
# for both nucleotide and amino acid alignments)
char_mat <- toupper(as.matrix(aln))

# non-gap indicator and pairwise count of shared non-gap columns.
# Treat both '-' and '.' as gap characters (common in AA alignments).
gap_chars   <- c("-", ".")
nongap      <- (!(char_mat %in% gap_chars)) + 0   # numeric n x L
dim(nongap) <- dim(char_mat)
both_nongap <- tcrossprod(nongap)                 # n x n

# Sum of co-occurrence matrices over each non-gap residue
residues <- setdiff(unique(as.vector(char_mat)), gap_chars)
matches  <- matrix(0, n, n)
for (a in residues) {
  ind         <- (char_mat == a) + 0
  dim(ind)    <- dim(char_mat)
  matches <- matches + tcrossprod(ind)
}

pid_mat <- matches / both_nongap
pid_mat[!is.finite(pid_mat)] <- NA   # pairs with zero shared columns
diag(pid_mat) <- 1
dimnames(pid_mat) <- list(ids, ids)

# Long-format table of unique pairs (upper triangle)
pid_tbl <- which(upper.tri(pid_mat), arr.ind = TRUE) %>%
  as_tibble() %>%
  rename(i = row, j = col) %>%
  mutate(
    seq1     = ids[i],
    seq2     = ids[j],
    identity = round(pid_mat[cbind(i, j)], 3)
  )

# Save matrix and long table
as_tibble(round(pid_mat, 3), rownames = "sequence") %>%
  write_csv("results/HD2_pairwise_identity_matrix.csv")

pid_tbl
# ---- 3. Distribution plots ---------------------------------
pid_df <- pid_tbl %>% select(identity)

# Top: histogram
HD2_plot <-
ggplot(pid_df, aes(x = identity)) +
  geom_histogram(breaks = seq(0.55, 1, 0.025)) +
  geom_jitter(aes(x=identity, y=-1000), alpha=0.02, width = 0, size=1, height = 500) +
  xlab("pairwise identity") +
  ylab("count") +
  theme_bw()

# ---- 4. Group sequences by an identity cutoff --------------
# Single-linkage clustering: two sequences are placed in the
# same group if a chain of pairs each with identity >= cutoff
# connects them (i.e., connected components of the graph whose
# edges are pairs with identity >= cutoff).
#
# We get this by feeding (1 - identity) as a distance to
# hclust(method = "single") and cutting at height (1 - cutoff).
# `identity_cutoff` is set near the top of this script (default 0.9)
# and may be overridden via the second command-line argument.

dist_mat <- 1 - pid_mat
dist_mat[is.na(dist_mat)] <- 1     # treat undefined pairs as fully diverged
diag(dist_mat) <- 0

hc      <- hclust(as.dist(dist_mat), method = "single")
cluster <- cutree(hc, h = 1 - identity_cutoff)

groups_tbl <-
  tibble(sequence = ids, group = unname(cluster)) %>%
  arrange(sequence) %>%
  mutate(group = factor(group, levels=unique(group)),
         group = as.numeric(group)
  )

# Per-group summary
# Two-column padded text file in the same style as results/MATA_types.txt
# (sequence name left-aligned, group number after)
name_width <- max(nchar(groups_tbl$sequence)) + 1
groups_tbl %>%
  arrange(group, sequence) %>%
  mutate(line = sprintf(paste0("%-", name_width, "s %d"), sequence, group)) %>%
  pull(line) %>%
  writeLines("results/HD2_groups.txt")

#####

plot_grid(HD1_plot, HD2_plot, ncol=1, labels="AUTO")
ggsave("plots/MATA_identidy_dist.pdf", width = 6, height = 8)

full_join(read_table("results/HD1_groups.txt", col_names = c("sequence", "HD1_group")), read_table("results/HD2_groups.txt", col_names = c("sequence", "HD2_group"))) %>%
  arrange(sequence) %>%
  write_csv("results/MATA_groups.csv")
