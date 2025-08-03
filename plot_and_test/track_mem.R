library(pryr)
library(profmem)


output1 <- profmem({
  set.seed(10)
  FSM_ARG.decimal(20L, 10/5e6, 5e6L, bacteria = TRUE, delta = 1000, optimise_recomb = TRUE)
  })

outoutput2 <- profmem({
  set.seed(11)
  FSM_ARG.decimal(20L, 10/5e6, 5e6L, bacteria = TRUE, delta = 1000, optimise_recomb = TRUE)
})

output3 <- profmem({
  set.seed(10)
  FSM_ARG.decimal(20L, 20/5e6, 5e6L, bacteria = TRUE, delta = 1000, optimise_recomb = TRUE)
})

output4 <- profmem({
  set.seed(11)
  FSM_ARG.decimal(20L, 20/5e6, 5e6L, bacteria = TRUE, delta = 1000, optimise_recomb = TRUE)
})
