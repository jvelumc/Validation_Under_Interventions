library(data.table)

parameters <- list(
  n = 1000,
  k = 5,
  # paramaters
  U.mu = 0,
  U.sd = 0.1,
  L0.mu = 0,
  L0.sd = 1,
  
  gamma.0 = -1,
  gamma.L = 0.5,
  
  #model for hazard
  alpha.0 = -2,
  alpha.A = -0.5,
  alpha.L = 0.5,
  alpha.U = 0.5
)
list2env(parameters, .GlobalEnv)

generate_events <- function(A, L, U) {
  
  T.obs <- rep(NA_real_, n)
  for (i in 1:k) {
    haz <- exp(alpha.0+alpha.A*A[,i]+alpha.L*L[,i]+alpha.U*U)
    potential_event_time <- rexp(n, rate = haz)
    T.obs <- fifelse(
      test = is.na(T.obs) & potential_event_time < 1, 
      yes = i - 1 + potential_event_time,
      no = T.obs
    )
  }
  D.obs <- fifelse(is.na(T.obs), 0, 1)
  T.obs <- fifelse(is.na(T.obs), 5, T.obs)
  
  return(list("D" = D.obs, "time" = T.obs))
}

build_data <- function(seed = 123) {
  
  set.seed(seed)
  
  # simulate ----------------------------------------------------------------
  
  U <- rnorm(n, U.mu, U.sd)
  A <- matrix(nrow = n, ncol = k)
  L <- matrix(nrow = n, ncol = k)
  
  L[, 1] <- rnorm(n, L0.mu + U, L0.sd)
  A[, 1] <- rbinom(n, 1, plogis(gamma.0 + gamma.L * L[, 1]))
  
  for (i in 2:k) {
    L[, i] <- rnorm(n, 0.8 * L[, i-1] - A[, i-1] + 0.1*i + U, 1)
    A[, i] <- fifelse(
      test = A[, i-1] == 1,
      yes = 1,
      no = rbinom(n, 1, plogis(gamma.0 + gamma.L * L[,i]))
    )
  }
  
  # simulate events
  events <- generate_events(A, L, U)
  D.obs <- events$D
  T.obs <- events$time
  
  # build dataset
  colnames(A) <- paste0("A", 0:4)
  colnames(L) <- paste0("L", 0:4)
  data <- data.table(id = 1:n, T.obs, D.obs, A, L, U)
  data[, `:=`(
    A1 = fifelse(T.obs < 1, 0, A1),
    A2 = fifelse(T.obs < 2, 0, A2),
    A3 = fifelse(T.obs < 3, 0, A3),
    A4 = fifelse(T.obs < 4, 0, A4)
  )] # no treatment after censoring/event
  return(data[])
}

make_long <- function(data) {
  Anames <- paste0("A", 0:4)
  Lnames <- paste0("L", 0:4)
  data_long <- melt(
    data, 
    id.vars = c("id", "T.obs", "D.obs", "U"), 
    measure.vars = list(A = Anames, L = Lnames),
    variable.name = "visit"
  )
  data_long[, visit := as.numeric(visit)]
  data_long[, time := visit - 1]
  data_long <- data_long[order(id, time)]
  
  data_long <- data_long[time < T.obs]
  data_long[, time.stop := pmin(time + 1, T.obs)]
  data_long[, event := fifelse(time.stop == T.obs & D.obs == 1, 1, 0)]
  
  data_long[, `:=`(
    Alag1 = shift(A, n = 1, fill = 0),
    Alag2 = shift(A, n = 2, fill = 0),
    Alag3 = shift(A, n = 3, fill = 0),
    Alag4 = shift(A, n = 4, fill = 0),
    Llag1 = shift(L, n = 1, fill = 0),
    Llag2 = shift(L, n = 2, fill = 0),
    Llag3 = shift(L, n = 3, fill = 0),
    Llag4 = shift(L, n = 4, fill = 0)
  ), by = id]
  
  data_long[, X := first(L), by = id]
  return(data_long[])
}

build_cf_data <- function(df_dev, seed = 345) {

  set.seed(seed)

  A0 <- matrix(0, nrow = n, ncol = k)
  A1 <- matrix(1, nrow = n, ncol = k)
  
  L0 <- matrix(nrow = n, ncol = k)
  L1 <- matrix(nrow = n, ncol = k)
  
  L0[,1] <- df_dev$L0
  L1[,1] <- df_dev$L0
  
  for (i in 2:k) {
    L0[, i] <- rnorm(n, 0.8 * L0[, i-1] - A0[, i-1] + 0.1*i + df_dev$U, 1)
    A0[, i] <- 0
    L1[, i] <- rnorm(n, 0.8 * L1[, i-1] - A1[, i-1] + 0.1*i + df_dev$U, 1)
    A1[, i] <- 1
  }
  
  events_0 <- generate_events(A0, L0, df_dev$U)  
  events_1 <- generate_events(A1, L1, df_dev$U)  
  
  cf_table <- data.table(
    id = 1:n,
    X = df_dev$L0,
    T0 = events_0$time,
    D0 = events_0$D,
    T1 = events_1$time,
    D1 = events_1$D
  )
  
  return(cf_table)
}
