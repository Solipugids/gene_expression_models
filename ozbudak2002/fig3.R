## monte carlo simulation
library(GillespieSSA2)
initial_state <- c(DNA = 1, mRNA = 0, Protein = 0)
params <- c(k_R = 0.01, k_P = 1, delta_R = 0.1, delta_P = 0.002)

reactions <- list(
  reaction(~k_R * DNA, c(mRNA = +1)),
  reaction(~k_P * mRNA, c(Protein = +1)),
  reaction(~delta_R * mRNA, c(mRNA = -1)),
  reaction(~delta_P * Protein, c(Protein = -1))
)
out <- 
  ssa(
    initial_state = initial_state,
    reactions = reactions,
    params = params,
    method = ssa_exact(),
    final_time = 3 * 60 * 60,
    census_interval = .01,
    verbose = FALSE
  )
plot_ssa(out)

## deterministic model
tlast <- 3 * 60 * 60
dt <- 1
gamma_R <- 0.1
gamma_P <- 0.002
k_R <- 0.01
k_P <- 1

timepoints <- seq(1, tlast)
R_total <- rep(0, length(timepoints))
P_total <- rep(0, length(timepoints))

R <- R_total[1]
P <- P_total[1]
for (i in seq(2, tlast)){
  dr_dt <- k_R - (gamma_R * R)
  dr <- dr_dt * dt
  R <- dr + R
  dp_dt <- k_P * R - gamma_P * P
  dp <- dp_dt * dt
  P <- dp + P
  R_total[i] <- R
  P_total[i] <- P
}
plot(timepoints, P_total)
