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

# Figure 3.5 (B)
v_field <- flowField(derivative, 
                     parameters = c(1, 0.8, 1.2, 1, 0.05, 1, 0.05, 0.2, 1, 0.5),
                     xlim = c(0, 1), 
                     ylim = c(0, 1),
                     system = "two.dim",
                     add = FALSE
                     )
traj <- trajectory(derivative,
                   y0 = c(1, 0.1),
                   parameters = c(1, 0.8, 1.2, 1, 0.05, 1, 0.05, 0.2, 1, 0.5),
                   tlim = c(0, 800),
                   system = "two.dim"
)

# Figure 3.5 (C)
nullc <- nullclines(derivative, 
                    parameters = c(1, 0.8, 1.2, 1, 0.05, 1, 0.05, 0.2, 1, 0.5),
                    xlim = c(0, 1), 
                    ylim = c(0, 1),
                    system = "two.dim",
                    col = c("black", "black"),
                    lwd = 2,
                    add.legend = FALSE,
                    add = FALSE
                    )
mf <- drawManifolds(derivative,
                    parameters = c(1, 0.8, 1.2, 1, 0.05, 1, 0.05, 0.2, 1, 0.5),
                    y0 = c(0.1, 0.1),
                    tstep = 0.01,
                    tend = 200,
                    add.legend = FALSE,
                    col = c("black", rgb(red = 0, green = 0, blue = 0, alpha = 0)),
                    lty = 2
)
y0 <- matrix(c(0, 0.9, 0, 0.7, 0, 0.4, 1, 0.1, 1, 0.4, 1, 0.8), 6, 2, byrow = TRUE)
traj <- trajectory(derivative,
                   y0 = y0,
                   parameters = c(1, 0.8, 1.2, 1, 0.05, 1, 0.05, 0.2, 1, 0.5),
                   tlim = c(0, 800),
                   system = "two.dim",
                   col = "gray"
                   )
