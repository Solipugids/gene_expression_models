library(GillespieSSA2)
initial_state <- c(R = 1, A = 0, M = 0, P = 0)
params <- c(k_on = 10, k_off = 10, s_A = 50, s_R = 5, s_P = 0.2, delta_M = 0.1, delta_P = 0.05)

reactions <- list(
  reaction(~k_on * R, c(R = -1, A = +1)),
  reaction(~k_off * A, c(R = +1, A = -1)),
  reaction(~s_A * A, c(M = +1)),
  reaction(~s_R * R, c(M = +1)),
  reaction(~s_P * M, c(P = +1)),
  reaction(~delta_M * M, c(M = -1)),
  reaction(~delta_P * P, c(P = -1))
)
out <- 
  ssa(
    initial_state = initial_state,
    reactions = reactions,
    params = params,
    method = ssa_exact(),
    final_time = 1200,
    census_interval = .001,
    verbose = FALSE
  )
