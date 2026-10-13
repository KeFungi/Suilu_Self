library(tidyverse)
library(adegenet)

genlight_dropInv <-
  function(genlight){
    genlight[,!(colSums(as.matrix(genlight)) %in% c(0, 2*dim(genlight)[1]))]
  }

s132ad_bed <-
  read.PLINK("plink/s132_admixture.raw")

s132oad_bed <-
  read.PLINK("plink/s132o_admixture.raw")

s132ad_pop <-
  seppop(s132ad_bed)

s132oad_pop <-
  seppop(s132oad_bed)

pop_list <-
  c("AU-other", "AU-WA", "CEu-Belgium", "NAm-MA", "NZ", "NZo", "SAm-Bariloche")

nsite <- 200
gl2RMEStext <-
  function(gl, size=1000){
    set.seed(12345)
    raw_matrix <-
      gl %>%
      genlight_dropInv() %>%
      .[,sample(1:dim(.)[2], size=size)] %>%
      as.matrix()

    het_matrix <-
     apply(raw_matrix, MARGIN = 2,
           function(x){
              case_when(x==0|x==2 ~ 0,
                        x==1 ~ 1,
                        is.na(x) ~ -99
              )
            }
     )

    text_block <-
      het_matrix %>%
      t() %>%
      as.data.frame() %>%
      map(~paste(.x, collapse = " ")) %>%
      reduce(~c(.x, .y))

    return(text_block)
  }

s_pops <-
  c(s132ad_pop, s132oad_pop) %>%
  .[pop_list]

blocks <-
  map(s_pops, ~gl2RMEStext(.x, nsite))

pop_size <-
  map(s_pops, ~dim(.x)[1])

title_line <-
  length(s_pops)

pop_lines <-
  mapply(c, pop_list, pop_size, nsite, SIMPLIFY = T) %>% as.vector()

all_lines <-
  c(title_line, pop_lines, unlist(blocks))

write_lines(all_lines, paste0("RMES/RMES_", nsite, ".txt"), sep="\r\n")
