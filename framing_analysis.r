source('functions.r')

library(mediation)

data(framing)
W <- framing %>% dplyr:::select(age, educ, gender, income)
W <- data.frame(model.matrix( ~ 0 + ., W))
colnames(W) <- paste0('W', 1:ncol(W))
A <- pull(framing, tone) + pull(framing, eth)
M <- pull(framing, emo) %>% as.numeric()#cut(breaks = c(0, 5, 8, 13)) %>% as.numeric()
M <- M - 2
Z <- pull(framing, p_harm) %>% as.numeric()#cut(breaks = c(0, 4, 6, 8)) %>% as.numeric()
Z <- Z - 1
Y <- pull(framing, immigr)

set.seed(2983645)
estimates <- estimate(Y, M, Z, A, W, 'xgbTree', 100)

library(xtable)
xtable(estimates, digits = c(0,0,3,3,3))

save.image('framing.rda')
