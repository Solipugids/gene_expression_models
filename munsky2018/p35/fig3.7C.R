library(phaseR)

# initialized ODE
derivative <- function(t, y, parameters){
  k_1 <- parameters[1]
  k_2 <- parameters[2]
  k_3 <- parameters[3]
  k_4 <- parameters[4]
  K_4 <- parameters[5]
  k_5 <- parameters[6]
  K_5 <- parameters[7]
  S <- parameters[8]
  Y <- parameters[9]
  E <- parameters[10]
  dx <- k_1 * S + k_2 * y[2] - k_3 * y[1]
  dy_p <- (k_4 * y[1] * (Y - y[2]) / (K_4 + Y - y[2])) - (k_5 * E * y[2] / (K_5 + y[2]))
  list(c(dx, dy_p))
}

S_con <- seq(0, 0.4, 0.001)
y0_m <- matrix(c(0, 0, 0.5, 0.5, 1 ,1), 3, 2, byrow = TRUE)
df <- data.frame(S = 0, x = 0)
for (i in 1:length(S_con)){
  for (j in 1:3) {
    fp <- findEquilibrium(derivative, 
                          y0 = y0_m[j, ],
                          parameters = c(1, 0.8, 1.2, 1, 0.05, 1, 0.05, S_con[i], 1, 0.5),
                          system = "two.dim",
                          summary = FALSE,
                          plot.it = FALSE
    )
    if (is.list(fp)){
      if (fp$classification=="Stable node"){
        df <- rbind(df, c(S_con[i], fp$ystar[1, 1]))
      }
    }
  }
}

plot(df, pch = 20, cex = 0.5, ylab = "Steady state response concentration x", main = "Hysteresis in a bistable system", xaxs="i", yaxs="i")