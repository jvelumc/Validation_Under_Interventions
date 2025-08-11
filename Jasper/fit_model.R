library(survival)

source("Jasper/build_simdata.R")

run <- function(seed) {
  df_dev <- build_data(seed)
  df_dev_long <- make_long(df_dev)
  
  # weird weighting model due to way treatment is assigned
  weight_model <- glm(A ~ L, data = df_dev_long[Alag1 == 0, ], family = "binomial")
  pred <- predict.glm(weight_model, newdata = df_dev_long, type = "response")
  df_dev_long[, w := fifelse(A == 1, pred, 1 - pred)]
  df_dev_long[, w := fifelse(Alag1 == 1, 1, w)]
  df_dev_long[, ipw := 1/cumprod(w), by = id]
  
  cox.msm <- coxph(
    formula = Surv(time, time.stop, event) ~ A + Alag1 + Alag2 + Alag3 + Alag4 + X,
    data = df_dev_long,
    weights = df_dev_long$w.cum
  )
  
  pred_int_risk <- function(msm, A, X) {
    
    baselinehazard <- survfit(
      msm, 
      newdata = data.frame(A = 0, Alag1 = 0, Alag2 = 0, Alag3 = 0, Alag4 = 0, X = 0),
    )
    cumbasehaz <- stepfun(baselinehazard$time, c(0, baselinehazard$cumhaz))
    
    if (A == 0) {
      cumhazard <- cumbasehaz(5) * exp(msm$coefficients[["X"]]*X)
    }
    else if (A == 1) {
      cumhazard <- 
        cumbasehaz(1) * exp(
          cox.msm$coefficients[["X"]]*X + 
            sum(cox.msm$coefficients[1])
        ) + 
        ( cumbasehaz(2) - cumbasehaz(1) ) * exp(
          cox.msm$coefficients[["X"]]*X + 
            sum(cox.msm$coefficients[1:2]) 
        ) + 
        ( cumbasehaz(3) - cumbasehaz(2) ) * exp(
          cox.msm$coefficients[["X"]]*X + 
            sum(cox.msm$coefficients[1:3])
        ) + 
        ( cumbasehaz(4) - cumbasehaz(3) ) * exp(
          cox.msm$coefficients[["X"]]*X + 
            sum(cox.msm$coefficients[1:4])
        ) + 
        ( cumbasehaz(5) - cumbasehaz(4) ) * exp(
          cox.msm$coefficients[["X"]]*X + 
            sum(cox.msm$coefficients[1:5])
        )
    }
    else stop("unkown A")
    
    surv <- exp(-cumhazard)
    return(1 - surv)
  }
  
  pred_int_risk(cox.msm, 0, 0)
  pred_int_risk(cox.msm, 1, 0)
  
  # should be different script from this point onwards
  
  df_val <- build_data(seed = seed + 1000)
  df_val_long <- make_long(df_val)
  
  # 'observed' risk in 'untreated' pop
  
  # first make weighting model: weighting the individuals * in such a way that it
  # represents the pop when all patients has been untreated
  # * in a population where patients get censored as soon as they get treated.
  
  # Fit it in Alag1 == 0, because when Alag = 1 then patient remains treated 
  # (prob = 1), due to simulation set up. 
  # 
  
  # df_val_long[, artificial_censor_time := min(.SD[A == 1, visit], Inf), by = id]
  
  
  mod.wt <- glm(A ~ L, 
                family = "binomial",
                data = df_val_long[Alag1 == 0, ])
  
  # probability of remaining untreated (this doesnt make sense for patients who 
  # are already treated, but these are artificially censored in next step anyway)
  
  # df_val_long[, wt0 := 1 - predict(mod.wt, type = "response", newdata = df_val_long)]
  # df_val_long[, ipw0 := 1/cumprod(wt0), by = id]
  # df_val_long[, in.dat.0 := A == 0]
  
  
  pred <- predict.glm(mod.wt, newdata = df_val_long, type = "response")
  df_val_long[, w.0 := fifelse(A == 1, 0, 1 - pred)]
  df_val_long[, w.0 := fifelse(Alag1 == 1, 1, w.0)]
  df_val_long[, ipw0 := 1/cumprod(w.0), by = id]
  
  df_val_long[, w.1 := fifelse(A == 1, pred, 0)]
  df_val_long[, w.1 := fifelse(Alag1 == 1, 1, w.1)]
  df_val_long[, ipw1 := 1/cumprod(w.1), by = id]
  
  
  df_val_long[, in_val_0 := A == 0]
  df_val_long[, in_val_1 := .SD[, first(A)] == 1, by = id]
  
  # estimate observed risk using KM
  
  km.0 <- survfit(
    Surv(time, time.stop, event) ~ 1,
    data = df_val_long[in_val_0 == TRUE, ],
    weights = ipw0
  )
  
  km.1 <- survfit(
    Surv(time, time.stop, event) ~ 1,
    data = df_val_long[in_val_1 == TRUE, ],
    weights = ipw1
  )
  
  observed_risk0 <- 1 - summary(km.0, time = 5)$surv
  observed_risk1 <- 1 - summary(km.1, time = 5)$surv
  df_cfdata <- build_cf_data(df_val, seed = seed + 2000)
  cf_real_risk0 <- 1-(survfit(
    Surv(T0, D0) ~ 1,
    data = df_cfdata,
  ) |> summary(time = 5))$surv
  
  cf_real_risk1 <- 1-(survfit(
    Surv(T1, D1) ~ 1,
    data = df_cfdata,
  ) |> summary(time = 5))$surv
  
  return(data.table(observed_risk0, cf_real_risk0, observed_risk1, cf_real_risk1))
}

l_results <- lapply(1:100, function(x) {print(x); run(x)})
df_results <- rbindlist(l_results)

df_results[, lapply(.SD, mean)]
