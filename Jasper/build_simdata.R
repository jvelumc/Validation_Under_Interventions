library(data.table)

n <- 1000
k <- 5
# paramaters
U.mu <- 0
U.sd <- 2
L0.mu <- 10 
L0.sd <- 4

gamma.0 <- -2
gamma.L <- 0.1

#model for hazard
alpha.0 <- 0.2
alpha.A <- -0.04
alpha.L <- 0.01
alpha.U <- 0.01
tfac <- 0.2

censoring <- TRUE

set.seed(123)

# simulate ----------------------------------------------------------------

U <- rnorm(n, U.mu, U.sd)
A <- matrix(nrow = n, ncol = k)
L <- matrix(nrow = n, ncol = k)


L[, 1] <- rnorm(n, L0.mu + U, L0.sd)
A[, 1] <- rbinom(n, 1, plogis(gamma.0 + gamma.L * L[, 1]))

for (i in 2:k) {
  L[, i] <- rnorm(n, 0.8 * L[, i-1] - A[, i-1] + 0.1*i + U, 4)
  A[, i] <- fifelse(
    test = A[, i-1] == 1,
    yes = 1,
    no = rbinom(n, 1, plogis(gamma.0 + gamma.L * L[,i]))
  )
}

# simulate events
T.obs <- rep(NA_real_, n)

for (i in 1:k) {
  haz <- alpha.0+alpha.A*A[,i]+alpha.L*(1-(i-1)*tfac)*L[,i]+alpha.U*U
  potential_event_time <- rexp(n, rate = haz)
  T.obs <- fifelse(
    test = is.na(T.obs) & potential_event_time < 1, 
    yes = i - 1 + potential_event_time,
    no = T.obs
  )
}
D.obs <- fifelse(is.na(T.obs), 0, 1)
T.obs <- fifelse(is.na(T.obs), 5, T.obs)

# simulate noninformative censoring
if (censoring == TRUE) {
  C <- rexp(n, rate = 0.1)
  T.obs <- fifelse(T.obs > C , C, T.obs)
  D.obs <- fifelse(T.obs < C & T.obs < 5, 1, 0)
}

# build dataset
colnames(A) <- paste0("A", 0:4)
colnames(L) <- paste0("L", 0:4)
data <- data.table(id = 1:n, T.obs, D.obs, A, L)
data[, `:=`(
  A1 = fifelse(T.obs < 1, 0, A1),
  A2 = fifelse(T.obs < 2, 0, A2),
  A3 = fifelse(T.obs < 3, 0, A3),
  A4 = fifelse(T.obs < 4, 0, A4)
)] # no treatment after censoring/event

# make long
data_long <- melt(
  data, 
  id.vars = c("id", "T.obs", "D.obs"), 
  measure.vars = list(A = colnames(A), L = colnames(L)),
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
