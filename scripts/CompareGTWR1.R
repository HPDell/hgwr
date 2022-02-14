# ## Experiment 1
experiment.current <- "compare-gtwr-1"
generate.data <- function (size) {
  set.seed(11)
  U <- rnorm(n = size, mean = 3000, sd = 100)
  set.seed(12)
  V <- rnorm(n = size, mean = 3000, sd = 100)
  set.seed(21)
  x1 <- rnorm(n = size, mean = 0, sd = 1)
  set.seed(22)
  x2 <- rnorm(n = size, mean = 0, sd = 1)
  set.seed(23)
  x3 <- rnorm(n = size, mean = 0, sd = 1)
  b0 <- ((U - 3000)/100) + ((V - 3000)/100)^2 - 3
  b1 <- ((U - 3000)/100) + ((V - 3000)/100)^2
  b2 <- ((U - 3000)/100) + 5 * ((V - 3000)/100)^2 + 2
  b3 <- -((U - 3000)/100) + ((V - 3000)/100)^2
  set.seed(1)
  y <- b0 + b1 * x1 + b2 * x2 + b3 *x3 + rnorm(n = size, mean = 0, sd = 1)
  list(
    data = data.frame(y = y, x1 = x1, x2 = x2, x3 = x3),
    coord = cbind(U = U, V = V),
    beta = data.frame(Intercept = b0, x1 = b1, x2 = b2, x3 = b3)
  )
}
simulate.data <- generate.data(100)
simulate.formula <- y ~ x1 + x2 + x3
simulate <- cbind(simulate.data$data, simulate.data$coord, simulate.data$beta)
colnames(simulate) <- c("y", paste("x", c(1:3), sep = ""), "U", "V", "Intercept", paste("b", c(1:3), sep = ""))
head(simulate)

simulate.r2 <- 1 - 
  sum((simulate.data$data$y - rowSums(as.matrix(cbind(1, simulate.data$data[,2:4])) * as.matrix(simulate.data$beta)))^2) / 
  sum((simulate.data$data$y - mean(simulate.data$data$y))^2)
simulate.r2

simulate.coord.minmax <- apply(simulate.data$coord, 2, max) - apply(simulate.data$coord, 2, min)
kernel.list <- list(
  gwdr.make.kernel(simulate.coord.minmax[1] * 0.618, adaptive = FALSE),
  gwdr.make.kernel(simulate.coord.minmax[2] * 0.618, adaptive = FALSE)
)

lambda <- 0.9999


write.table(cbind(1.0,simulate.data$data[,2:4]), "gwr_x.csv", sep = " ", row.names = F, col.names = F)
write.table(simulate.data$data[,1], "gwr_y.csv", sep = " ", row.names = F, col.names = F)
write.table(simulate.data$coord, "gwr_u.csv", sep = " ", row.names = F, col.names = F)
write.table(simulate.data$beta, "gwr_beta.csv", sep = " ", row.names = F, col.names = F)


library(GWmodel)
data.gwr <- simulate[,1:6]
coordinates(data.gwr) <- ~ U + V

gwr.result <- gwr.basic(simulate.formula, data.gwr, bw = 270.0, kernel = "gaussian")
gwr.result$GW.diagnostic
head(gwr.result$SDF@data)
