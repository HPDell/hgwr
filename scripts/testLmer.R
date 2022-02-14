library(lmerTest)

generate.data <- function(size.group) {
  size.x <- sum(size.group)
  group.num <- length(size.group)
  set.seed(21)
  x1 <- rnorm(n = size.x, mean = 0, sd = 1)
  set.seed(22)
  x2 <- rnorm(n = size.x, mean = 0, sd = 1)
  set.seed(23)
  x3 <- rnorm(n = size.x, mean = 0, sd = 1)
  set.seed(31)
  g1 <- rnorm(n = group.num, mean = 0, sd = 1)
  set.seed(32)
  g2 <- rnorm(n = group.num, mean = 0, sd = 1)
  gamma00 <- 1
  gamma01 <- 2
  gamma02 <- 3
  gamma10 <- 4
  gamma20 <- 5
  gamma30 <- 6
  set.seed(40)
  nita0 <- rnorm(n = group.num, mean = 0, sd = 1)
  set.seed(41)
  nita1 <- numeric(group.num)
  set.seed(42)
  nita2 <- numeric(group.num)
  set.seed(42)
  nita3 <- rnorm(n = group.num, mean = 0, sd = 1)
  beta_list <- list()
  for (j in 1:group.num) {
    beta_0j <- gamma00 + gamma01 * g1[j] + gamma02 * g2[j] + nita0[j]
    beta_1j <- gamma10 + nita1[j]
    beta_2j <- gamma20 + nita2[j]
    beta_3j <- gamma30 + nita3[j]
    beta_j <- c(beta_0j, beta_1j, beta_2j, beta_3j)
    beta_list[[j]] <- matrix(beta_j, nrow = size.group[j], ncol = 4, byrow = T)
  }
  beta <- Reduce(rbind,  beta_list)
  set.seed(0)
  epsilon <- rnorm(n = size.x, mean = 0, sd = 1)
  y <- rowSums(beta * as.matrix(cbind(1, x1, x2, x3))) + epsilon
  list(
    y = y,
    x = cbind(1, x1, x2, x3),
    g = cbind(1, g1, g2),
    beta = beta,
    group = Reduce(c, Map(rep, 1:group.num, size.group))
  )
}

size.group <- c(10, 20, 30, 40)
simulate.data <- generate.data(size.group)

write.table(simulate.data$x, "./data/hlm_x.csv", row.names = F, col.names = F, sep = ",")
write.table(simulate.data$g, "./data/hlm_g.csv", row.names = F, col.names = F, sep = ",")
write.table(simulate.data$y, "./data/hlm_y.csv", row.names = F, col.names = F, sep = ",")
write.table(simulate.data$beta, "./data/hlm_beta.csv", row.names = F, col.names = F, sep = ",")
write.table(simulate.data$group - 1, "./data/hlm_group.csv", row.names = F, col.names = F, sep = ",")

hlm.data.g <- Reduce(rbind, apply(cbind(size.group, simulate.data$g), MARGIN = 1, FUN = function(x) {
  gn = x[1]
  g = x[-1]
  matrix(g, nrow = gn, ncol = length(g), byrow = T)
}))

hlm.data <- as.data.frame(cbind(as.numeric(simulate.data$y), hlm.data.g[,-1], simulate.data$x[,-1]))
colnames(hlm.data) <- c("y", "g1", "g2", "x1", "x2", "x3")
hlm.data$group = paste("group", simulate.data$group, sep = "")
head(hlm.data)

hlm.f <- formula(y ~ g1 + g2 + x1 + x2 + x3 + (x3 | group))
hlm <- lmer(hlm.f, data = hlm.data)
coefficients(hlm)$group

data("sleepstudy", package="lme4")
sleepstudy
