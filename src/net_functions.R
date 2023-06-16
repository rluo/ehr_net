
require(IsingSampler)
require(tidyverse)
require(infotheo)
require(glmnet)
require(qgraph)


generate_nextstatus <- function(current_status, true_beta, thresh) {
    # p = exp(current_status%*%true_beta)/(1+exp(current_status%*%true_beta))
    p <- exp(current_status %*% true_beta + thresh) / (1 + exp(current_status %*% true_beta + thresh))
    next_dis <- rbinom(1, 1, p)
    return(next_dis)
}

fn_simulate <- function(init_state, n_obs, true_beta_m, thresh) {
    current_status <- init_state
    sim <- init_state
    for (i in 1:(n_obs - 1)) {
        thresh_dis <- thresh[i]
        next_status <- apply(true_beta_m, MARGIN = 2, FUN = generate_nextstatus, current_status = current_status, thresh = thresh_dis)
        sim <- rbind(sim, next_status)
        current_status <- next_status
    }

    sim <- cbind(Time = 1:n_obs, sim)
    rownames(sim) <- NULL
    return(sim)
}

IsingCogitator <- function(nSample, Graph, Thresh_store) {
    # Gen's method
    dat.sim <- NULL
    # dat.sim = rbind(dat.sim, rep(0, dim(Graph)[1]), rep(1, dim(Graph)[1]))

    # generate the initial state, make sure at least one disease
    for (i in 1:nSample) {
        repeat{
            repeat{
                init_state <- rbinom(dim(Graph)[1], 1, 0.2)
                if (sum(init_state) != 0) break
            }
            n_obs <- sample(3:6, 1, prob = c(0.3, 0.3, 0.2, 0.2))
            pts_encounter1 <- fn_simulate(init_state, n_obs, Graph, Thresh_store)
            pts_encounter2 <- tail(pts_encounter1, 2)
            pts_encounter3 <- pts_encounter2[, -1]
            pts_encounter4 <- pts_encounter3[1, ] + pts_encounter3[2, ]
            pts_encounter5 <- as.integer(pts_encounter4 > 0)

            if (sum(pts_encounter5) != 0 & sum(pts_encounter5) != dim(Graph)[1]) break
        }
        dat.sim <- rbind(dat.sim, pts_encounter5)
    }
    rownames(dat.sim) <- NULL
    return(dat.sim)
}





# Elastic Net eBIC
Elastic_fit <- function(x, family = "binomial", AND = TRUE, gamma = 0.25,
                        plot = TRUE, progressbar = TRUE, lowerbound.lambda = NA, alpha = 0.9,
                        ...) {
    t0 <- Sys.time()
    xx <- x
    if (family != "binomial") {
          stop("This procedure is currently only supported for binary (family='binomial') data")
      }
    checklognet <- function(y) {
        res <- c()
        y <- as.factor(y)
        ntab <- table(y)
        minclass <- min(ntab)
        if (minclass <= 1) {
              res <- 0
          } else {
            res <- 1
        }
        return(res)
    }
    NodesToAnalyze <- apply(x, 2, checklognet) != 0
    names(NodesToAnalyze) <- colnames(x)
    if (!any(NodesToAnalyze)) {
          stop("No variance in dataset")
      }
    if (any(!NodesToAnalyze)) {
        warning(paste(
            "Nodes with too little variance (not allowed):",
            paste(colnames(x)[!NodesToAnalyze], collapse = ", ")
        ))
    }
    x <- as.matrix(x)
    allthemeans <- colMeans(x)
    x <- x[, NodesToAnalyze, drop = FALSE]
    nvar <- ncol(x)
    p <- nvar - 1
    intercepts <- betas <- lambdas <- list(vector, nvar)
    nlambdas <- rep(0, nvar)

    for (i in 1:nvar) {
        a <- glmnet(x[, -i], x[, i], family = family, alpha = alpha)
        intercepts[[i]] <- a$a0
        betas[[i]] <- a$beta
        lambdas[[i]] <- a$lambda
        nlambdas[i] <- length(lambdas[[i]])
    }

    if (progressbar == TRUE) {
          pb <- txtProgressBar(max = nvar, style = 3)
      }
    P <- logl <- sumlogl <- J <- matrix(0, max(nlambdas), nvar)
    for (i in 1:nvar) {
        J[1:ncol(betas[[i]]), i] <- colSums(betas[[i]] != 0)
    }
    logl_M <- P_M <- array(0, dim = c(
        nrow(x), max(nlambdas),
        nvar
    ))
    N <- nrow(x)
    for (i in 1:nvar) {
        betas.ii <- as.matrix(betas[[i]])
        int.ii <- intercepts[[i]]
        y <- matrix(0, nrow = N, ncol = ncol(betas.ii))
        xi <- x[, -i]
        NB <- nrow(betas.ii)
        for (bb in 1:NB) {
            y <- y + betas.ii[rep(bb, N), ] * xi[, bb]
        }
        y <- matrix(int.ii, nrow = N, ncol = ncol(y), byrow = TRUE) +
            y
        n_NA <- max(nlambdas) - ncol(y)
        if (n_NA > 0) {
            for (vv in 1:n_NA) {
                y <- cbind(y, NA)
            }
        }
        P_M[, , i] <- exp(y * x[, i]) / (1 + exp(y))
        logl_M[, , i] <- log(P_M[, , i])
        if (progressbar == TRUE) {
              setTxtProgressBar(pb, i)
          }
    }
    logl_Msum <- colSums(logl_M, 1, na.rm = FALSE)
    if (progressbar == TRUE) {
          close(pb)
      }
    sumlogl <- logl_Msum
    sumlogl[sumlogl == 0] <- NA
    penalty <- J * log(nrow(x)) + 2 * gamma * J * log(p)
    EBIC <- -2 * sumlogl + penalty
    lambda.mat <- matrix(NA, nrow(EBIC), ncol(EBIC))
    for (i in 1:nvar) {
        lambda.mat[, i] <- c(lambdas[[i]], rep(NA, nrow(EBIC) -
            length(lambdas[[i]])))
    }
    if (!is.na(lowerbound.lambda)) {
        EBIC <- EBIC / (lambda.mat >= lowerbound.lambda) * 1
    }
    lambda.opt <- apply(EBIC, 2, which.min)
    lambda.val <- rep(NA, nvar)
    thresholds <- 0
    for (i in 1:length(lambda.opt)) {
        lambda.val[i] <- lambda.mat[lambda.opt[i], i]
        thresholds[i] <- intercepts[[i]][lambda.opt[i]]
    }

    weights.opt <- matrix(, nvar, nvar)
    for (i in 1:nvar) {
        weights.opt[i, -i] <- betas[[i]][, lambda.opt[i]]
    }
    asymm.weights <- weights.opt
    diag(asymm.weights) <- 0

    # absolute min
    AB_MIN_graph <- asymm.weights
    AB_MIN_graph_ele <- ifelse(abs(AB_MIN_graph[lower.tri(AB_MIN_graph)]) < abs(t(AB_MIN_graph)[lower.tri(t(AB_MIN_graph))]),
        AB_MIN_graph[lower.tri(AB_MIN_graph)], t(AB_MIN_graph)[lower.tri(t(AB_MIN_graph))]
    )
    AB_MIN_graph[lower.tri(AB_MIN_graph)] <- AB_MIN_graph_ele
    AB_MIN_graph[upper.tri(AB_MIN_graph)] <- t(AB_MIN_graph)[upper.tri(AB_MIN_graph)]
    AB_MIN_graph_final <- AB_MIN_graph


    if (AND == TRUE) {
        adj <- weights.opt
        adj <- (adj != 0) * 1
        EN.weights <- adj * t(adj)
        EN.weights <- EN.weights * weights.opt
        meanweights.opt <- (EN.weights + t(EN.weights)) / 2
        diag(meanweights.opt) <- 0
    }
    else {
        meanweights.opt <- (weights.opt + t(weights.opt)) / 2
        diag(meanweights.opt) <- 0
    }
    graphNew <- matrix(0, length(NodesToAnalyze), length(NodesToAnalyze))
    graphNew[NodesToAnalyze, NodesToAnalyze] <- meanweights.opt
    colnames(graphNew) <- rownames(graphNew) <- colnames(xx)
    threshNew <- ifelse(allthemeans > 0.5, -Inf, Inf)
    threshNew[NodesToAnalyze] <- thresholds
    if (plot == TRUE) {
          notplot <- FALSE
      } else {
        notplot <- TRUE
    }
    q <- qgraph(graphNew,
        layout = "spring", labels = names(NodesToAnalyze),
        DoNotPlot = notplot, ...
    )
    Res <- list(
        weiadj = graphNew, thresholds = threshNew, q = q,
        gamma = gamma, AND = AND, time = Sys.time() - t0, asymm.weights = asymm.weights,
        lambda.values = lambda.val, AB_MIN_graph_final = AB_MIN_graph_final
    )
    class(Res) <- "IsingFit"
    return(Res)
}





raw_weight <- function(Ising_Data) {
    N <- dim(Ising_Data)[2]
    raw_w <- matrix(0, nrow = N, ncol = N)
    for (i in 1:nrow(Ising_Data)) {
        encounter_1 <- Ising_Data[i, ]
        encounter_2 <- Ising_Data[i, ]
        raw_w[encounter_1 > 0, encounter_2 > 0] <- raw_w[encounter_1 > 0, encounter_2 > 0] + 1
    }
    return(raw_w)
}



convention_methods <- function(Ising_Data) {
    raw_weight_matrix <- raw_weight(Ising_Data)
    N <- dim(Ising_Data)[2]


    OR_2by2_matrix_raw <- matrix(NA, N, N)
    OR_P_matrix <- matrix(NA, N, N)

    Phi_matrix_raw <- matrix(NA, N, N)
    Phi_P_matrix <- matrix(NA, N, N)

    OER_matrix_raw <- matrix(NA, N, N)
    OER_P_matrix <- matrix(NA, N, N)

    for (i in 1:N) {
        for (j in 1:N) {
            if (i == j) {
                OR_2by2_matrix_raw[i, j] <- 0
                Phi_matrix_raw[i, j] <- 0
                OER_matrix_raw[i, j] <- 0
            }
            else {
                a <- raw_weight_matrix[i, j]
                b <- raw_weight_matrix[i, i] - raw_weight_matrix[i, j]
                c <- raw_weight_matrix[j, j] - raw_weight_matrix[i, j]
                d <- dim(Ising_Data)[1] - a - b - c
                if (a == 0) {
                    a <- a + 1
                }
                if (b == 0) {
                    b <- b + 1
                }
                if (c == 0) {
                    c <- c + 1
                }
                if (d == 0) {
                    d <- d + 1
                }

                con_table <- matrix(c(a, b, c, d), 2, 2, byrow = TRUE)
                OR_2by2_matrix_raw[i, j] <- log((a * d) / (b * c))
                OR_P_matrix[i, j] <- fisher.test(con_table)$p.value


                Phi_matrix_raw[i, j] <- (a * d - b * c) / sqrt((a + b) * (a + c) * (b + d) * (c + d))
                Phi_n <- max((a + b), (a + c))
                Phi_t <- Phi_matrix_raw[i, j] * sqrt(Phi_n - 2) / sqrt(1 - Phi_matrix_raw[i, j]**2)
                Phi_P_matrix[i, j] <- 2 * pt(-abs(Phi_t), df = Phi_n - 1)



                OER_matrix_raw[i, j] <- log((a * (a + b + c + d)) / ((a + b) * (a + c)))
                OER_log_std <- 1 / a + 1 / ((a + b) * (a + c)) - 1 / dim(Ising_Data)[1] - 1 / (dim(Ising_Data)[1]**2)
                OER_Z <- OER_matrix_raw[i, j] / OER_log_std
                OER_P_matrix[i, j] <- 2 * pnorm(-abs(OER_Z))
            }
        }
    }

    OR_P_vector_adj <- p.adjust(as.vector(OR_P_matrix), method = "bonferroni", n = sum(!is.na(OR_P_matrix))) # method = "fdr"
    OR_P_matrix_adj <- matrix(OR_P_vector_adj, nrow = N, ncol = N)
    OR_P_matrix_adj[is.na(OR_P_matrix_adj)] <- 1
    OR_P_Pvalue_mask <- ifelse(OR_P_matrix_adj < 0.05, 1, 0)
    OR_2by2_matrix <- OR_2by2_matrix_raw * OR_P_Pvalue_mask

    OR_out_adj_M <- OR_2by2_matrix
    OR_out_adj_M_ele <- ifelse(abs(OR_out_adj_M[lower.tri(OR_out_adj_M)]) < abs(t(OR_out_adj_M)[lower.tri(t(OR_out_adj_M))]),
        OR_out_adj_M[lower.tri(OR_out_adj_M)], t(OR_out_adj_M)[lower.tri(t(OR_out_adj_M))]
    )
    OR_out_adj_M[lower.tri(OR_out_adj_M)] <- OR_out_adj_M_ele
    OR_out_adj_M[upper.tri(OR_out_adj_M)] <- t(OR_out_adj_M)[upper.tri(OR_out_adj_M)]
    OR_2by2_matrix_final <- OR_out_adj_M



    Phi_P_vector_adj <- p.adjust(as.vector(Phi_P_matrix), method = "bonferroni", n = sum(!is.na(Phi_P_matrix))) # method = "fdr"
    Phi_P_matrix_adj <- matrix(Phi_P_vector_adj, nrow = N, ncol = N)
    Phi_P_matrix_adj[is.na(Phi_P_matrix_adj)] <- 1
    Phi_P_Pvalue_mask <- ifelse(Phi_P_matrix_adj < 0.05, 1, 0)
    Phi_matrix <- Phi_matrix_raw * Phi_P_Pvalue_mask

    Phi_out_adj_M <- Phi_matrix
    Phi_out_adj_M_ele <- ifelse(abs(Phi_out_adj_M[lower.tri(Phi_out_adj_M)]) < abs(t(Phi_out_adj_M)[lower.tri(t(Phi_out_adj_M))]),
        Phi_out_adj_M[lower.tri(Phi_out_adj_M)], t(Phi_out_adj_M)[lower.tri(t(Phi_out_adj_M))]
    )
    Phi_out_adj_M[lower.tri(Phi_out_adj_M)] <- Phi_out_adj_M_ele
    Phi_out_adj_M[upper.tri(Phi_out_adj_M)] <- t(Phi_out_adj_M)[upper.tri(Phi_out_adj_M)]
    Phi_matrix_final <- Phi_out_adj_M



    OER_P_vector_adj <- p.adjust(as.vector(OER_P_matrix), method = "bonferroni", n = sum(!is.na(OER_P_matrix))) # method = "fdr"
    OER_P_matrix_adj <- matrix(OER_P_vector_adj, nrow = N, ncol = N)
    OER_P_matrix_adj[is.na(OER_P_matrix_adj)] <- 1
    OER_P_Pvalue_mask <- ifelse(OER_P_matrix_adj < 0.05, 1, 0)
    OER_matrix <- OER_matrix_raw * OER_P_Pvalue_mask

    OER_out_adj_M <- OER_matrix
    OER_out_adj_M_ele <- ifelse(abs(OER_out_adj_M[lower.tri(OER_out_adj_M)]) < abs(t(OER_out_adj_M)[lower.tri(t(OER_out_adj_M))]),
        OER_out_adj_M[lower.tri(OER_out_adj_M)], t(OER_out_adj_M)[lower.tri(t(OER_out_adj_M))]
    )
    OER_out_adj_M[lower.tri(OER_out_adj_M)] <- OER_out_adj_M_ele
    OER_out_adj_M[upper.tri(OER_out_adj_M)] <- t(OER_out_adj_M)[upper.tri(OER_out_adj_M)]
    OER_matrix_final <- OER_out_adj_M

    return(list(OR_2by2_matrix_final, Phi_matrix_final, OER_matrix_final))
}

mi_pair <- function(v1, v2, n_perm = 5000) {
    val <- mutinformation(v1, v2)
    perm_vals <- sapply(1:n_perm, function(j) {
        set.seed(j + 100)
        mutinformation(v1, sample(v2))
    })
    return(c(val, mean(perm_vals > val)))
}

mi_conn <- function(Ising_Data) {
    ## nsubs vs p codes
    p <- ncol(Ising_Data)
    ret <- array(0, c(2, p, p))
    for (i in 1:(p - 1)) {
        for (j in (i + 1):p) {
            tmp <- mi_pair(Ising_Data[, i], Ising_Data[, j])
            ret[, i, j] <- tmp
            ret[, j, i] <- tmp
        }
    }
    return(ret[1, , ] * (ret[2, , ] < 0.05))
}


simu_data <- function(iter, N = 20, nSample = 5000, sparse_rate = 0.2) {
    sample_multiplier <- 1.5
    random_seed <- 50 + iter
    set.seed(random_seed)
    # Ising parameters:
    Graph1 <- matrix(runif(N^2, -2, 2), N, N)
    Graph2 <- pmax(Graph1, t(Graph1)) ## symmetrize matrix

    sparse_filter <- matrix(0, N, N)
    sparse_filter[lower.tri(sparse_filter)] <- sample(0:1, (N * (N - 1) / 2), TRUE, prob = c((1 - sparse_rate), sparse_rate))
    sparse_filter[upper.tri(sparse_filter)] <- t(sparse_filter)[upper.tri(sparse_filter)]

    Graph <- sparse_filter * Graph2
    Thresh_store <- -rowSums(Graph) / 3

    # Simulate:, method ='CFTP'
    Ising_mtx <- IsingSampler((nSample * sample_multiplier), Graph, Thresh_store)

    Ising_Data <- as.data.frame(Ising_mtx)
    Ising_Data["Total"] <- rowSums(Ising_Data)
    Ising_Data <- Ising_Data %>%
        filter(Total < 100 & Total > 0) %>%
        select(-Total) %>%
        slice_head(n = nSample)
    Ising_Data <- as.matrix(Ising_Data)

    return(list(graph = Graph, data = Ising_Data))
}