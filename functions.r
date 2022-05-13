library(tidyverse)
library(origami)
library(caret)
library(data.table)
library(caretEnsemble)
library(mvtnorm)

estimate <- function(Y, M, Z, A, W, model, p) {


    cf_fun <- function(fold, data, outcome, predictors, type, method = 'xgbTree', tuneLength){

        if(outcome == 'A') pred_vals <- NULL
        if(outcome == 'Z') pred_vals <- expand.grid(A = vA)
        if(outcome == 'M') pred_vals <- expand.grid(A = vA, Z = vZ)
        if(outcome == 'Y') pred_vals <- expand.grid(A = vA, Z = vZ, M = vM)

        predictors <- c(predictors, names(data)[substr(names(data), 1, 1) == 'W'])

        train_data <- training(data)
        valid_data <- validation(data)

        train_X <- train_data[, predictors]
        valid_X <- valid_data[, predictors]
        if(type == 'reg'){
            train_Y <- train_data[, outcome]
            type_pred <- 'raw'
            classProbs <- FALSE
        } else {
            train_Y <- factor(train_data[, outcome])
            ## levels(train_Y) <- paste0(outcome, levels(train_Y))
            type_pred <- 'prob'
            classProbs <- TRUE
        }

        fit <- train(train_X, train_Y, method = method, tuneLength = tuneLength,
                     trControl =  caret::trainControl(method = "cv", number = 5, search = 'random', verboseIter = TRUE))

        if(outcome != 'A') {
            pred_datasets <- lapply(1:nrow(pred_vals), function(i) {
                if(outcome == 'Z')pred_data <- valid_data %>% mutate(A = pred_vals[i, 'A'])
                if(outcome == 'M')pred_data <- valid_data %>% mutate(Z = pred_vals[i, 'Z'], A = pred_vals[i, 'A'])
                if(outcome == 'Y')pred_data <- valid_data %>% mutate(M = pred_vals[i, 'M'],
                                                                     Z = pred_vals[i, 'Z'], A = pred_vals[i, 'A'])
                pred_data <- pred_data[, predictors]
                pred <- predict(fit, pred_data, type = type_pred)
                return(pred)
            })
        } else {
            pred_datasets <- predict(fit, valid_data[, predictors], type = type_pred)
        }

        out <- list(pred = pred_datasets, id = valid_data[, 'id'], pred_vals = pred_vals)

        return(out)

    }

    regress <- function(data, outcome, predictors, type, method, tuneLength, folds){

        cv_results <- cross_validate(
            cv_fun = cf_fun, folds = folds, data = data,
            outcome = outcome, predictors = predictors, type = type, method = method, tuneLength = tuneLength
        )

    }

    post_processM <- function(pred){
        obj <- function(a, z, m) {
            vals <- unique(pred$pred_vals)
            i <- which(vals$A == a & vals$Z == z)
            out <- lapply(pred$pred, function(x)x[[i]])
            out <- rbindlist(out)
            out <- data.frame(id = pred$id, out)
            return(out[order(out$id), m + 1])
        }
    }

    post_processY <- function(pred){
        obj <- function(a, z, m) {
            vals <- unique(pred$pred_vals)
            i <- which(vals$A == a & vals$Z == z & vals$M == m)
            out <- lapply(pred$pred, function(x)x[[i]])
            out <- unlist(out)
            out <- data.frame(id = pred$id, out)
            return(out[order(out$id), 2])
        }
    }

    post_processZ <- function(pred){
        obj <- function(a, z) {
            vals <- unique(pred$pred_vals)
            i <- which(vals$A == a)
            out <- lapply(pred$pred, function(x)x[[i]])
            out <- rbindlist(out)
            out <- data.frame(id = pred$id, out)
            return(out[order(out$id), z + 1])
        }
    }

    post_processA <- function(pred){
        obj <- function(a) {
            out <- data.frame(id = pred$id, pred$pred)
            return(out[order(out$id), a + 2])
        }
    }

    mY <- function(A, Z, M) {
        out <- rowSums(
            sapply(1:nrow(vals), function(i){
                a <- vals[i, 'A']
                z <- vals[i, 'Z']
                m <- vals[i, 'M']
                funY(a, z, m) * (A == a) * (Z == z) * (M == m)
            })
        )
        return(out)
    }

    pM <- function(A, Z, M) {
        out <- rowSums(sapply(1:nrow(vals), function(i){
            a <- vals[i, 'A']
            z <- vals[i, 'Z']
            m <- vals[i, 'M']
            funM(a, z, m) * (A == a) * (Z == z) * (M == m)
        }))
        return(out)
    }

    pZ <- function(A, Z) {
        out <- rowSums(
            sapply(1:nrow(valsAZ), function(i){
                z <- valsAZ[i, 'Z']
                a <- valsAZ[i, 'A']
                funZ(a, z) * (A == a) * (Z == z)
            })
        )
        return(out)
    }

    pA <- function(A) {
        out <- rowSums(sapply(vA, function(a){
            funA(a) * (A == a)
        }))
        return(out)
    }

    pM.Z <- function(A, M) rowSums(sapply(vA, function(z) pM(A, z, M) * pZ(A, z)))

    ## Compute EIFs

    data <- data.frame(id = 1:length(Y), Y, M, Z, A, W)
    folds <- make_folds(data)

    vA <- unique(A)
    vZ <- unique(Z)
    vM <- unique(M)

    vals <- expand.grid(A = vA, Z = vZ, M = vM)
    valsAZ <- expand.grid(A = vA, Z = vZ)
    valsAM <- expand.grid(A = vA, M = vZ)

    ## Fit models
    library(doParallel)
    cl <- makePSOCKcluster(5)
    registerDoParallel(cl)

    predA <- regress(data, 'A', NULL, 'class', model, p, folds)
    predZ <- regress(data, 'Z', 'A', 'class', model, p, folds)
    predM <- regress(data, 'M', c('A', 'Z'), 'class', model, p, folds)
    predY <- regress(data, 'Y', c('A', 'Z', 'M'), 'reg', model, p, folds)
    funA <- post_processA(predA)
    funZ <- post_processZ(predZ)
    funM <- post_processM(predM)
    funY <- post_processY(predY)

    stopCluster(cl)

    h10A <- rowSums(sapply(vA, function(a) a * pM(a, Z, M) * pZ(a, Z) * pA(a)))
    h10  <- rowSums(sapply(vA, function(a) pM(a, Z, M) * pZ(a, Z) * pA(a)))

    h11A <- rowSums(sapply(vA, function(a) a * pM.Z(a, M) * pZ(a, Z) * pA(a)))
    h11  <- rowSums(sapply(vA, function(a) pM.Z(a, M) * pZ(a, Z) * pA(a)))

    h21A <- rowSums(sapply(vA, function(a) a * pM.Z(a, M) * pA(a)))
    h21  <- rowSums(sapply(vA, function(a) pM.Z(a, M) * pA(a)))

    h32A <- rowSums(sapply(vA, function(a) rowSums(sapply(vZ, function(z)a * pM(a, z, M) * pZ(a, z) * pA(a)))))
    h32 <- rowSums(sapply(vA, function(a) rowSums(sapply(vZ, function(z)pM(a, z, M) * pZ(a, z) * pA(a)))))

    h30A <- rowSums(sapply(vA, function(a) a * pM(a, Z, M) * pA(a)))
    h30  <- rowSums(sapply(vA, function(a) pM(a, Z, M) * pA(a)))

    h40A <- rowSums(sapply(vA, function(a) a * pA(a)))
    h40  <- rowSums(sapply(vA, function(a) pA(a)))


    eif10A <- h10A / (pM(A, Z, M) * pZ(A, Z)) * (Y - mY(A, Z, M)) +
        rowSums(sapply(1:nrow(vals), function(i){
            vals[i, 'A'] * mY(A, vals[i, 'Z'], vals[i, 'M']) * pM(vals[i, 'A'], vals[i, 'Z'], vals[i, 'M']) *
                pZ(vals[i, 'A'], vals[i, 'A']) * pA(vals[i, 'A'])
        })) -
        rowSums(sapply(vA, function(aprime){
            rowSums(sapply(1:nrow(vals), function(i){
                vals[i, 'A'] * mY(aprime, vals[i, 'Z'], vals[i, 'M']) * pM(vals[i, 'A'], vals[i, 'Z'], vals[i, 'M']) *
                    pZ(vals[i, 'A'], vals[i, 'A']) * pA(vals[i, 'A'])
            })) * pA(aprime)
        })) +
        A * rowSums(sapply(vA, function(aprime)mY(aprime, Z, M) * pA(aprime)))

    eif10 <- h10 / (pM(A, Z, M) * pZ(A, Z)) * (Y - mY(A, Z, M)) +
        rowSums(sapply(1:nrow(vals), function(i){
            mY(A, vals[i, 'Z'], vals[i, 'M']) * pM(vals[i, 'A'], vals[i, 'Z'], vals[i, 'M']) *
                pZ(vals[i, 'A'], vals[i, 'A']) * pA(vals[i, 'A'])
        })) -
        rowSums(sapply(vA, function(aprime){
            rowSums(sapply(1:nrow(vals), function(i){
                mY(aprime, vals[i, 'Z'], vals[i, 'M']) * pM(vals[i, 'A'], vals[i, 'Z'], vals[i, 'M']) *
                    pZ(vals[i, 'A'], vals[i, 'A']) * pA(vals[i, 'A'])
            })) * pA(aprime)
        })) +
        rowSums(sapply(vA, function(aprime)mY(aprime, Z, M) * pA(aprime)))

    eif11A <- h11A / (pM(A, Z, M) * pZ(A, Z)) * (Y - mY(A, Z, M)) +
        rowSums(sapply(1:nrow(vals), function(i){
            vals[i, 'A'] * mY(A, vals[i, 'Z'], vals[i, 'M']) * pM.Z(vals[i, 'A'], vals[i, 'M']) *
                pZ(vals[i, 'A'], vals[i, 'A']) * pA(vals[i, 'A'])
        })) -
        rowSums(sapply(vA, function(aprime){
            rowSums(sapply(1:nrow(vals), function(i){
                vals[i, 'A'] * mY(aprime, vals[i, 'Z'], vals[i, 'M']) * pM.Z(vals[i, 'A'], vals[i, 'M']) *
                    pZ(vals[i, 'A'], vals[i, 'A']) * pA(vals[i, 'A'])
            })) * pA(aprime)
        })) +
        A * rowSums(sapply(1:nrow(valsAM), function(i)mY(valsAM[i, 'A'], Z, valsAM[i, 'M']) *
                                                      pM.Z(A, valsAM[i, 'M']) * pA(valsAM[i, 'A']))) -
        A * rowSums(sapply(1:nrow(vals), function(i)mY(vals[i, 'A'], vals[i, 'Z'], vals[i, 'M']) *
                                                    pM.Z(A, vals[i, 'M']) * pZ(A, vals[i, 'Z']) * pA(vals[i, 'A']))) +
        mean(A * rowSums(sapply(1:nrow(valsAZ), function(i)mY(valsAZ[i, 'A'], valsAZ[i, 'Z'], M) *
                                                           pZ(A, valsAZ[i, 'Z']) * pA(valsAZ[i, 'A']))))

    eif11 <- h11 / (pM(A, Z, M) * pZ(A, Z)) * (Y - mY(A, Z, M)) +
        rowSums(sapply(1:nrow(vals), function(i){
            mY(A, vals[i, 'Z'], vals[i, 'M']) * pM.Z(vals[i, 'A'], vals[i, 'M']) *
                pZ(vals[i, 'A'], vals[i, 'A']) * pA(vals[i, 'A'])
        })) -
        rowSums(sapply(vA, function(aprime){
            rowSums(sapply(1:nrow(vals), function(i){
                mY(aprime, vals[i, 'Z'], vals[i, 'M']) * pM.Z(vals[i, 'A'], vals[i, 'M']) *
                    pZ(vals[i, 'A'], vals[i, 'A']) * pA(vals[i, 'A'])
            })) * pA(aprime)
        })) +
        rowSums(sapply(1:nrow(valsAM), function(i)mY(valsAM[i, 'A'], Z, valsAM[i, 'M']) *
                                                  pM.Z(A, valsAM[i, 'M']) * pA(valsAM[i, 'A']))) -
        rowSums(sapply(1:nrow(vals), function(i)mY(vals[i, 'A'], vals[i, 'Z'], vals[i, 'M']) *
                                                pM.Z(A, vals[i, 'M']) * pZ(A, vals[i, 'Z']) * pA(vals[i, 'A']))) +
        rowSums(sapply(1:nrow(valsAZ), function(i)mY(valsAZ[i, 'A'], valsAZ[i, 'Z'], M) *
                                                  pZ(A, valsAZ[i, 'Z']) * pA(valsAZ[i, 'A'])))

    eif21A <- h21A / pM(A, Z, M) * (Y - mY(A, Z, M)) +
        rowSums(sapply(1:nrow(valsAM), function(i){
            valsAM[i, 'A'] * mY(A, Z, valsAM[i, 'M']) * pM.Z(valsAM[i, 'A'], valsAM[i, 'M']) * pA(valsAM[i, 'A'])
        })) -
        rowSums(sapply(1:nrow(valsAZ), function(j){
            rowSums(sapply(1:nrow(valsAM), function(i){
                valsAM[i, 'A'] * mY(valsAZ[j, 'A'], valsAZ[j, 'Z'], valsAM[i, 'M']) *
                    pM.Z(valsAM[i, 'A'], valsAM[i, 'M']) * pA(valsAM[i, 'A']) *
                    pZ(valsAZ[j, 'A'], valsAZ[j, 'Z']) * pA(valsAZ[j, 'A'])
            }))
        })) +
        A * rowSums(sapply(1:nrow(valsAZ), function(j){
            mY(valsAZ[j, 'A'], valsAZ[j, 'Z'], M) * pZ(valsAZ[j, 'A'], valsAZ[j, 'Z']) * pA(valsAZ[j, 'A'])
        }))

    eif21 <- h21 / pM(A, Z, M) * (Y - mY(A, Z, M)) +
        rowSums(sapply(1:nrow(valsAM), function(i){
            mY(A, Z, valsAM[i, 'M']) * pM.Z(valsAM[i, 'A'], valsAM[i, 'M']) * pA(valsAM[i, 'A'])
        })) -
        rowSums(sapply(1:nrow(valsAZ), function(j){
            rowSums(sapply(1:nrow(valsAM), function(i){
                mY(valsAZ[j, 'A'], valsAZ[j, 'Z'], valsAM[i, 'M']) *
                    pM.Z(valsAM[i, 'A'], valsAM[i, 'M']) * pA(valsAM[i, 'A']) *
                    pZ(valsAZ[j, 'A'], valsAZ[j, 'Z']) * pA(valsAZ[j, 'A'])
            }))
        })) +
        rowSums(sapply(1:nrow(valsAZ), function(j){
            mY(valsAZ[j, 'A'], valsAZ[j, 'Z'], M) * pZ(valsAZ[j, 'A'], valsAZ[j, 'Z']) * pA(valsAZ[j, 'A'])
        }))

    eif32A <- h32A / pM(A, Z, M) * (Y - mY(A, Z, M)) +
        rowSums(sapply(1:nrow(vals), function(i){
            vals[i, 'A'] * mY(A, Z, vals[i, 'M']) * pM(vals[i, 'A'], vals[i, 'Z'], vals[i, 'M']) *
                pZ(A, vals[i, 'Z']) * pA(vals[i, 'A'])
        })) -
        rowSums(sapply(1:nrow(valsAZ), function(j){
            rowSums(sapply(1:nrow(vals), function(i){
                vals[i, 'A'] * mY(valsAZ[j, 'A'], valsAZ[j, 'Z'], vals[i, 'M']) * pM(vals[i, 'A'], vals[i, 'Z'], vals[i, 'M']) *
                    pZ(valsAZ[j, 'A'], vals[i, 'Z']) * pA(vals[i, 'A']) *
                    pZ(valsAZ[j, 'A'], valsAZ[j, 'Z']) * pA(valsAZ[j, 'A'])
            }))
        })) +
        A / pZ(A, Z) * rowSums(sapply(1:nrow(valsAZ), function(j){
            mY(valsAZ[j, 'A'], valsAZ[j, 'Z'], M) * pZ(valsAZ[j, 'A'], valsAZ[j, 'Z']) *
                pZ(valsAZ[j, 'A'], Z) * pA(valsAZ[j, 'A'])
        })) -
        A / pZ(A, Z) * rowSums(sapply(1:nrow(vals), function(j){
            mY(vals[j, 'A'], vals[j, 'Z'], vals[j, 'M']) * pZ(vals[j, 'A'], vals[j, 'Z']) *
                pM(A, Z, vals[j, 'M']) * pZ(vals[j, 'A'], Z) * pA(vals[j, 'A'])
        })) +
        rowSums(sapply(1:nrow(vals), function(j){
            vals[j, 'A'] * mY(A, vals[j, 'Z'], vals[j, 'M']) * pZ(A, vals[j, 'Z']) *
                pM(vals[j, 'A'], Z, vals[j, 'M']) * pA(vals[j, 'A'])
        })) -
        rowSums(sapply(vZ, function(zprime) {
            rowSums(sapply(1:nrow(vals), function(j){
                vals[j, 'A'] * mY(A, vals[j, 'Z'], vals[j, 'M']) * pZ(A, vals[j, 'Z']) *
                    pM(vals[j, 'A'], zprime, vals[j, 'M']) * pA(vals[j, 'A']) * pZ(A, zprime)
            }))
        })) +
        A * rowSums(sapply(vZ, function(zprime) {
            rowSums(sapply(1:nrow(vals), function(j){
                mY(vals[j, 'A'], vals[j, 'Z'], vals[j, 'M']) * pZ(vals[j, 'A'], vals[j, 'Z']) *
                    pM(A, zprime, vals[j, 'M']) * pA(vals[j, 'A']) * pZ(vals[j, 'A'], zprime)
            }))
        }))

    eif32 <- h32 / pM(A, Z, M) * (Y - mY(A, Z, M)) +
        rowSums(sapply(1:nrow(vals), function(i){
            mY(A, Z, vals[i, 'M']) * pM(vals[i, 'A'], vals[i, 'Z'], vals[i, 'M']) *
                pZ(A, vals[i, 'Z']) * pA(vals[i, 'A'])
        })) -
        rowSums(sapply(1:nrow(valsAZ), function(j){
            rowSums(sapply(1:nrow(vals), function(i){
                mY(valsAZ[j, 'A'], valsAZ[j, 'Z'], vals[i, 'M']) * pM(vals[i, 'A'], vals[i, 'Z'], vals[i, 'M']) *
                    pZ(valsAZ[j, 'A'], vals[i, 'Z']) * pA(vals[i, 'A']) *
                    pZ(valsAZ[j, 'A'], valsAZ[j, 'Z']) * pA(valsAZ[j, 'A'])
            }))
        })) +
        1 / pZ(A, Z) * rowSums(sapply(1:nrow(valsAZ), function(j){
            mY(valsAZ[j, 'A'], valsAZ[j, 'Z'], M) * pZ(valsAZ[j, 'A'], valsAZ[j, 'Z']) *
                pZ(valsAZ[j, 'A'], Z) * pA(valsAZ[j, 'A'])
        })) -
        1 / pZ(A, Z) * rowSums(sapply(1:nrow(vals), function(j){
            mY(vals[j, 'A'], vals[j, 'Z'], vals[j, 'M']) * pZ(vals[j, 'A'], vals[j, 'Z']) *
                pM(A, Z, vals[j, 'M']) * pZ(vals[j, 'A'], Z) * pA(vals[j, 'A'])
        })) +
        rowSums(sapply(1:nrow(vals), function(j){
            mY(A, vals[j, 'Z'], vals[j, 'M']) * pZ(A, vals[j, 'Z']) *
                pM(vals[j, 'A'], Z, vals[j, 'M']) * pA(vals[j, 'A'])
        })) -
        rowSums(sapply(vZ, function(zprime) {
            rowSums(sapply(1:nrow(vals), function(j){
                mY(A, vals[j, 'Z'], vals[j, 'M']) * pZ(A, vals[j, 'Z']) *
                    pM(vals[j, 'A'], zprime, vals[j, 'M']) * pA(vals[j, 'A']) * pZ(A, zprime)
            }))
        })) +
        rowSums(sapply(vZ, function(zprime) {
            rowSums(sapply(1:nrow(vals), function(j){
                mY(vals[j, 'A'], vals[j, 'Z'], vals[j, 'M']) * pZ(vals[j, 'A'], vals[j, 'Z']) *
                    pM(A, zprime, vals[j, 'M']) * pA(vals[j, 'A']) * pZ(vals[j, 'A'], zprime)
            }))
        }))

    eif30A <- h30A / pM(A, Z, M) * (Y - mY(A, Z, M)) +
        A / pZ(A, Z) * rowSums(sapply(vA, function(a){
            mY(a, Z, M) * pZ(a, Z) * pA(a)
        })) -
        A / pZ(A, Z) * rowSums(sapply(1:nrow(valsAM), function(j){
            mY(valsAM[j, 'A'], Z, valsAM[j, 'M']) * pM(A, Z, valsAM[j, 'M']) *
                pZ(valsAM[j, 'A'], Z) * pA(valsAM[j, 'A'])
        })) +
        rowSums(sapply(1:nrow(valsAM), function(j){
            valsAM[j, 'A'] * mY(A, Z, valsAM[j, 'M']) * pM(valsAM[j, 'A'], Z, valsAM[j, 'M']) *
                pA(valsAM[j, 'A'])
        })) -
        rowSums(sapply(1:nrow(valsAZ), function(i){
            rowSums(sapply(1:nrow(valsAM), function(j){
                valsAM[j, 'A'] * mY(valsAZ[i, 'A'], valsAZ[i, 'Z'], valsAM[j, 'M']) *
                    pM(valsAM[j, 'A'], valsAZ[i, 'Z'], valsAM[j, 'M']) * pA(valsAM[j, 'A'])  *
                    pZ(valsAZ[i, 'A'], valsAZ[i, 'Z']) * pA(valsAZ[i, 'A'])
            }))
        })) +
        A * rowSums(sapply(1:nrow(vals), function(j){
            mY(vals[j, 'A'], vals[j, 'Z'], vals[j, 'M']) * pM(A, vals[j, 'Z'], vals[j, 'M']) *
                pZ(vals[j, 'A'], vals[j, 'Z']) * pA(vals[j, 'A'])
        }))


    eif30 <- h30 / pM(A, Z, M) * (Y - mY(A, Z, M)) +
        1 / pZ(A, Z) * rowSums(sapply(vA, function(a){
            mY(a, Z, M) * pZ(a, Z) * pA(a)
        })) -
        1 / pZ(A, Z) * rowSums(sapply(1:nrow(valsAM), function(j){
            mY(valsAM[j, 'A'], Z, valsAM[j, 'M']) * pM(A, Z, valsAM[j, 'M']) *
                pZ(valsAM[j, 'A'], Z) * pA(valsAM[j, 'A'])
        })) +
        rowSums(sapply(1:nrow(valsAM), function(j){
            mY(A, Z, valsAM[j, 'M']) * pM(valsAM[j, 'A'], Z, valsAM[j, 'M']) *
                pA(valsAM[j, 'A'])
        })) -
        rowSums(sapply(1:nrow(valsAZ), function(i){
            rowSums(sapply(1:nrow(valsAM), function(j){
                mY(valsAZ[i, 'A'], valsAZ[i, 'Z'], valsAM[j, 'M']) *
                    pM(valsAM[j, 'A'], valsAZ[i, 'Z'], valsAM[j, 'M']) * pA(valsAM[j, 'A'])  *
                    pZ(valsAZ[i, 'A'], valsAZ[i, 'Z']) * pA(valsAZ[i, 'A'])
            }))
        })) +
        rowSums(sapply(1:nrow(vals), function(j){
            mY(vals[j, 'A'], vals[j, 'Z'], vals[j, 'M']) * pM(A, vals[j, 'Z'], vals[j, 'M']) *
                pZ(vals[j, 'A'], vals[j, 'Z']) * pA(vals[j, 'A'])
        }))

    mYW <- rowSums(sapply(1:nrow(vals), function(j) {
        mY(vals[j, 'A'], vals[j, 'Z'], vals[j, 'M']) * pM(vals[j, 'A'], vals[j, 'Z'], vals[j, 'M']) *
            pZ(vals[j, 'A'], vals[j, 'Z']) * pA(vals[j, 'A'])

    }))

    eif40A <- h40A * (Y - mYW) + A * mYW
    eif40  <- Y

    eif_fun <- function(eifA, eif) {
        eifA - mean(eifA) - (mean(A) * (eif - mean(eif)) + mean(eif) * (A - mean(A)))
    }

    param_fun <- function(AY, Y) {
        mean(AY) - mean(A) * mean(Y)
    }

    theta1 <- param_fun(A * Y, Y) - param_fun(eif10A, eif10)
    eif1 <- eif_fun(A * Y, Y) - eif_fun(eif10A, eif10)

    theta2 <- param_fun(eif11A, eif11) - param_fun(eif21A, eif21)
    eif2 <- eif_fun(eif11A, eif11) - eif_fun(eif21A, eif21)

    theta3 <- param_fun(eif21A, eif21) - param_fun(eif32A, eif32)
    eif3 <- eif_fun(eif21A, eif21) - eif_fun(eif32A, eif32)

    theta4 <- param_fun(eif30A, eif30) - param_fun(eif40A, eif40)
    eif4 <- eif_fun(eif30A, eif30) - eif_fun(eif40A, eif40)

    theta <- param_fun(A * Y, Y) - param_fun(eif40A, eif40)
    eif <- eif_fun(A * Y, Y) - eif_fun(eif40A, eif40)

    theta23 <- theta - (theta1 + theta2 + theta3 + theta4)
    eif23 <- eif - (eif1 + eif2 + eif3 + eif4)

    eifs <- data.frame(eif, eif1, eif2, eif3, eif4, eif23)
    ses <- sqrt(diag(var(eifs)) / length(Y))

    Smat <- cor(eifs)
    nsamples <- rmvnorm(1e6, mean = rep(0, 6), sigma = Smat)
    maxstat <- apply(nsamples, 1, function(x)max(abs(x)))
    ## qq <- quantile(maxstat, 1 - 0.05 / 2)
    qq <- qnorm(0.975)

    path  <-  c('Total', 'A->Y', 'A->Z->Y', 'A->Z->M->Y', 'A->M->Y', 'Int-Conf')
    effect  <-  c(theta, theta1, theta2, theta3, theta4, theta23)
    lci  <- effect - qq * ses
    uci  <- effect + qq * ses

    out <- data.frame(path, effect, lci, uci)

    return(out)

}
