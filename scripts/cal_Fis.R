library(tidyverse)
library(KRIS)

fis_pop <-
  function(plink_bed, pop1){
    pop1_pos <- which(plink_bed$ind.info$FamID==pop1)
    fis.nei(plink_bed$snp, pop1_pos)
  }

# Nei (1987) Fis, same as hierfstat basic.stats()$overall["Fis"]: 1 - mean(Ho)/mean(Hs)
fis.nei <-
  function (X, idx) {
    G = X[idx, , drop = FALSE]
    G[!G %in% 0:2] = NA
    n = colSums(!is.na(G))
    p1 = colSums(G, na.rm = T)/(2*n)
    Ho = colSums(G == 1, na.rm = T)/n
    Hs = n/(n-1) * (2*(p1-p1^2) - Ho/(2*n))
    ok = which(n >= 2)
    fis = 1 - mean(Ho[ok])/mean(Hs[ok])
    return(fis)
  }

fis_wrapper <-
  function(plink_bed){
    pop_min5 <-
      plink_bed$ind.info %>%
      group_by(FamID) %>%
      summarise(n=n()) %>%
      filter(n>=5) %>%
      pull(FamID)

    fis_vec <-
      map_dbl(pop_min5, ~fis_pop(plink_bed, .x))

    tibble(FamID=pop_min5) %>%
      bind_cols(tibble(Fis=fis_vec))
  }

s132ad_bed <-
  read.bed("plink/s132_admixture.bed",
           "plink/s132_admixture.bim",
           "plink/s132_admixture.fam",
           only.snp = FALSE)


s132oad_bed <-
  read.bed("plink/s132o_admixture.bed",
           "plink/s132o_admixture.bim",
           "plink/s132o_admixture.fam",
           only.snp = FALSE)

fis_s132 <-
  fis_wrapper(s132ad_bed)

fis_s132o <-
  fis_wrapper(s132oad_bed)

write_csv(rbind(fis_s132,fis_s132o) , "results/fis.csv")
