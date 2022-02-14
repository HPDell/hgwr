library(lmerTest)
library(MASS)
library(dplyr)

generate.data <- function(n, each, order = c(0,1,2), random = c(0,0,0,0)) {
    k.g <- 1
    k.h <- 2
    s <- ifelse(length(each) == 1, each, length(each))
    l <- ifelse(length(each) == 1, n*n*each, sum(each))
    if (length(each) == 1) {
        group <- rep(1:(n*n), each = each)
    } else {
        group <- rep(1:(n*n), each)
    }
    X <- mvrnorm(l, rep(0, k.g + k.h), diag(1, k.g + k.h, k.g + k.h))
    for (i in 1:(n*n)) {
        group.i <- which(group == i)
        X[group.i,1] <- mean(X[group.i,1])
    }
    X <- cbind(1, X)
    positions.data <- matrix(0, l, 2)
    positions.beta <- matrix(0, n*n, 2)
    for (i in 1:n) {
        for (j in 1:n) {
            p1 = (i - 1)*n + j
            positions.beta[p1,] = c(i,j)
            for (k in 1:each[p1]) {
                p2 = ifelse(p1 > 1, sum(each[1:(p1 - 1)]) + k, k)
                positions.data[p2,] = c(i,j)
            }
        }
    }
    beta.g <- matrix(0, n*n, k.g)
    for (i in 1:n) {
        for (j in 1:n) {
            p = (i - 1)*n + j
            beta = 0
            order.g = order[2]
            for (o1 in 0:order.g) {
                for (o2 in 0:order.g) {
                    set.seed(200 + i + j)
                    beta = beta + rnorm(1) * i^o1 * j^o2
                }
            }
            beta.g[p, 1] = 10 * beta / (n^(2*order.g))
        }
    }
    beta.g[, 1] = scale(beta.g[, 1] + random[2] * rnorm(n*n, sd = 0.1))
    beta.h <- matrix(0, n*n, k.h)
    for (i in 1:n) {
        for (j in 1:n) {
            p = (i - 1)*n + j
            beta = 0
            order.h = order[3]
            for (o1 in 0:order.h) {
                for (o2 in 0:order.h) {
                    set.seed(300 + i + j)
                    beta = beta + rnorm(1) * i^o1 * j^o2
                }
            }
            beta.h[p, 1] = beta / (n^(2*order.h))
        }
    }
    set.seed(12)
    beta.h[, 1] = scale(beta.h[, 1] + random[3] * rnorm(n*n, sd = 0.1))
    set.seed(13)
    beta.h[, 2] = 1 + random[4] * rnorm(n*n, sd = 1)
    beta.intercept <- numeric(n*n)
    for (i in 1:n) {
        for (j in 1:n) {
            p = (i - 1)*n + j
            beta = 0
            order.intercept = order[1]
            for (o1 in 0:order.intercept) {
                for (o2 in 0:order.intercept) {
                    set.seed(100 + i + j)
                    beta = beta + rnorm(1) * i^o1 * j^o2
                }
            }
            beta.intercept[p] = beta / (n^(2*order.intercept))
        }
    }
    set.seed(10)
    beta.intercept = scale(beta.intercept + random[1] * rnorm(n*n, sd = 0.1))
    beta <- cbind(beta.intercept, beta.g, beta.h)
    y <- numeric(l)
    for (i in 1:(n*n)) {
        beta.i <- beta[i,]
        group.i <- which(group == i)
        X.i <- X[group.i,]
        y[group.i] <- X.i %*% beta.i
    }
    y = y + rnorm(l)
    data <- as.data.frame(cbind(group, positions.data, y, X[,-1]))
    colnames(data) <- c("group", "lon", "lat", "y", paste0(rep("g", k.g), 1:k.g), paste0(rep("h", k.h), 1:k.h))
    beta <- as.data.frame(cbind(positions.beta, beta))
    colnames(beta) <- c("lon", "lat", "Intercept", paste0(rep("g", k.g), 1:k.g), paste0(rep("h", k.h), 1:k.h))
    list(
        data = data,
        beta = beta
    )
}

oi = 3
og = 2
oh = 2
ri = 0
rg = 0
rh1 = 1
rh2 = 0

set.seed(100)
gen.each <- floor(runif(15*15, min = 10, max = 50))
gen <- generate.data(15, gen.each, order = c(oi, og, oh), random = c(ri, rg, rh1, rh2))

data.group <- gen$data %>% group_by(group) %>% summarise(g1 = mean(g1), lon = mean(lon), lat = mean(lat)) %>% as.data.frame()

write.table(gen$data["y"], "./data/hlmgwr_y.csv", sep = ",", row.names = F, col.names = F)
write.table(cbind(1, gen$data[c("h1", "h2")]), "./data/hlmgwr_x.csv", sep = ",", row.names = F, col.names = F)
write.table(gen$data["group"] - 1, "./data/hlmgwr_group.csv", sep = ",", row.names = F, col.names = F)
write.table(data.group[c("lon", "lat")], "./data/hlmgwr_u.csv", sep = ",", row.names = F, col.names = F)
write.table(cbind(1, data.group["g1"]), "./data/hlmgwr_g.csv", sep = ",", row.names = F, col.names = F)

write.table(gen$beta, "./data/hlmgwr_beta.csv", sep = ",", row.names = F, col.names = F)

ngroup = 15*15
ndata = nrow(gen$data)
p = 2
q = 1

mask_np = matrix(c(1, 0), nrow = p, ncol = ngroup, byrow = F)
mask = rbind(mask_np)
write.table(mask, "./data/hlmgwr_mask.csv", sep = ",", row.names = F, col.names = F)



hlm.f <- y ~ 1 + g1 + h1 + h2 + (h1 | group)
hlm.model <- lmer(hlm.f, gen$data)

y <- gen$data[["y"]]
y.tss <- sum((y - mean(y))^2)

hlm.rss <- sum(residuals(hlm.model)^2)
hlm.pr2 <- 1 - hlm.rss / y.tss
cat("HLM R2:", hlm.pr2, "\n")

library(GWmodel)
gen.data <- gen$data
coordinates(gen.data) <- ~ lon + lat
gwr.formula <- y ~ g1 + h1 + h2
gwr.bw <- bw.gwr(gwr.formula, gen.data, kernel = "gaussian", adaptive = F, parallel.method = "omp")
gwr.model <- gwr.basic(gwr.formula, gen.data, bw = gwr.bw, kernel = "gaussian", adaptive = F, parallel.method = "omp")
gwr.R2 <- gwr.model$GW.diagnostic$gw.R2
cat("GWR R2:", gwr.R2, "\n")

