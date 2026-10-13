require(tidyverse)
cal_Lk <-
  function(input_set, k) {
    data <-
      input_set[which(input_set != 0)]

    n <- sum(data)

    if (k<sum(data!=0)) {
      L <- 0
    } else {

      data_sort <-
        sort(data) %>%
        c(rep(0, k-length(.)), .)

      freq <-
        table(data_sort)

      m <- factorial(k)/reduce(factorial(freq), ~.x*.y)
      P <- dmultinom(data_sort, prob = rep(1/k, k))
      L <- m*P
    }
    return(L)
  }

cal_Lk_range <-
  function(input_set, krange){
    map_dfr(krange, ~tibble(k=.x, L=cal_Lk(input_set, .x)))
  }

cal_MLk_range <-
  function(input_set, krange, simplify=TRUE){
    if (simplify) {
    cal_Lk_range(input_set, krange) %>%
      .[which.max(.$L),] %>%
        pull(k)
    } else {
      cal_Lk_range(input_set, krange) %>%
        .[which.max(.$L),]
    }
  }

cal_k_dup <-
  function(input_set, simplify=TRUE){
    data <-
      input_set[which(input_set != 0)]

    n <- sum(data)

    denominator <-
      (n*n-n)/2

    numerator <- sum(data*(data-1)/2)

    k <- (numerator/denominator)^-1
    if (simplify) {
        return(k)
      } else {
        return(tibble(k=k, d=numerator, pair=denominator))
      }
  }
