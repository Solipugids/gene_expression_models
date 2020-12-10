k_1 <- 1
k_2 <- 0.8
k_3 <- 1.2
k_4 <- 1
K_4 <- 0.05
k_5 <- 1
K_5 <- 0.05
S <- 0.2
Y <- 1
E <- 0.5
dx <- 0.01
dy <- 0.01
dt <- 0.01

# vector field
df <- data.frame(expand.grid(rep(list(seq(0.1, 1, 0.1)), 2)))
colnames(df) <- c("x", "y")
plot.new()
plot(10, xlab="x", ylab="y", xlim=c(0, 1.1), ylim=c(0, 1.1), xaxs="i", yaxs="i")
for (i in 1:nrow(df)) {
  x_i <- df[i, "x"]
  y_i <- df[i, "y"]
  dx <- k_1 * S + k_2 * y_i - k_3 * x_i 
  dy <- (k_4 * x_i * (Y - y_i) / (K_4 + Y - y_i)) - (k_5 * E * y_i) / (K_5 + y_i)
  x_j <- x_i + (dx * dt)
  y_j <- y_i + (dy * dt)
  if (dx != 0 & dy != 0){
    arrows(x_i, y_i, x_j, y_j, code = 2, angle = 15, length = 0.1)
  }
}

t_0 <- 0
tlast <- 10000
timepoint <- seq(t_0, (t_0 + dt * tlast), dt)
x_i <- 1
y_i <- 0
df <- data.frame(x = x_i, y = y_i)
for (m in 1:length(timepoint)){
  dx <- k_1 * S + k_2 * y_i - k_3 * x_i 
  dy <- (k_4 * x_i * (Y - y_i) / (K_4 + Y - y_i)) - (k_5 * E * y_i) / (K_5 + y_i)
  x_i <- x_i + dx * dt
  y_i <- y_i + dy * dt
  df <- rbind(df, c(x_i, y_i))
}
points(df, type = "l", lwd = 2)

# nullclines
x <- seq(0, 1, 0.1)
y_p <- (k_3 * x - k_1 * S) / k_2
plot(x, y_p, xlab = "x", ylab = "y", type = "l", xaxs="i", yaxs="i", xlim = c(0, 1), ylim = c(0, 1), lty = 2, lwd = 2)

y_p <- seq(0, 1, 0.01)
x <- ((k_5* E * y_p) * (K_4 + Y - y_p))/((K_5 + y_p) * (Y - y_p) * k_4)
points(x, y_p, type = "l", lwd = 2)

# different initial conditions of x and y
t_0 <- 0
tlast <- 800
timepoint <- seq(t_0, (t_0 + dt * tlast), dt)
init_cond <- data.frame(x = c(0, 0, 0, 1, 1, 1), y = c(0.4, 0.5, 0.9, 0.1, 0.5, 0.6))
for (m in 1:nrow(init_cond)) {
  x_i <- init_cond[m, "x"]
  y_i <- init_cond[m, "y"]
  df <- data.frame(x = x_i, y = y_i)
  for (n in 1:length(timepoint)){
    dx <- k_1 * S + k_2 * y_i - k_3 * x_i 
    dy <- (k_4 * x_i * (Y - y_i) / (K_4 + Y - y_i)) - (k_5 * E * y_i) / (K_5 + y_i)
    x_i <- x_i + dx * dt
    y_i <- y_i + dy * dt
    df <- rbind(df, c(x_i, y_i))
  }
  points(df, type = "l", lwd = 2, col = "gray")
  arrows(df[50, "x"], df[50, "y"], df[51, "x"], df[51, "y"], angle = 15, code = 2, length = 0.1, col = "gray")
}
