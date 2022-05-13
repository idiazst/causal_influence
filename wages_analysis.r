source('functions.r')

library(causalweight)

data(wexpect)
W <- wexpect %>% dplyr:::select(age, swiss, hassiblings, motherhighedu, fatherhighedu, motherworkedfull, motherworkedpart,
                                matwellbeing, homeowner, treatmentinformation)
W <- data.frame(model.matrix( ~ 0 + ., W))
colnames(W) <- paste0('W', 1:ncol(W))
A <- pull(wexpect, male)
## Msectr <- wexpect %>% select(starts_with('sector')) %>% max.col() %>% as.factor()
M <- wexpect %>% select(starts_with('plans')) %>% max.col()
## M <- model.matrix(~.^2, data = data.frame(Msectr, Mplans)) %>% max.col()
Z <- wexpect %>% select(business, econ, communi, businform) %>% max.col()
Y <- pull(wexpect, wexpect2)

set.seed(2983645)
estimates <- estimate(Y, M, Z, A, W, 'xgbTree', 100)

library(xtable)
xtable(estimates, digits = c(0,0,3,3,3))

save.image('wages.rda')
