new_model <- odin::odin({
  
  # define approach-dependent parameters
  h_0    <- if (AON > 0) adjusted_failure else 1
  ff_app <- if (AON > 0) ff               else ff * adjusted_failure
  
  # define time-dependent parameters
  lambda_t  <- if (t < tP) 0 else lambda
  ff_t      <- if (t < tP) 0 else ff_app
  
  
  # ODE system
  deriv(H) <- - (lambda_t + ff_t + gamma) * H
  deriv(S) <-                      gamma  * H - lambda_t * S
  deriv(I) <-   (lambda_t + ff_t        ) * H + lambda_t * S
  
  
  # initial conditions
  initial(H) <- h_0
  initial(S) <- 1 - h_0
  initial(I) <- 0
  
  
  # parameter values
  tP               <- user()
  AON              <- user()
  adjusted_failure <- user()
  ff               <- user()
  lambda           <- user()
  gamma            <- user()
  
})


# H + S
one_minus_bigF <- function (delta_t, h_tP, AON, gamma, lambda, ff, adjusted_failure) {
  
  if (delta_t < 0) {
    
    1
    
  } else {
    
    # define approach-dependent parameters
    if (AON > 0) {
      h_tP <- h_tP * adjusted_failure
    } else {
      ff   <- ff   * adjusted_failure
    }
    
    
    # result
    ((ff * (1 - h_tP) + gamma) * exp(-lambda * delta_t) + ff * h_tP * exp(-(gamma + lambda + ff) * delta_t)) / (ff + gamma)
    
  }
  
}

AON_one_minus_bigF <- function (delta_t, h_tP, gamma, lambda, ff, adjusted_failure) {
  
  # define approach-dependent parameters
  h_tP <- h_tP * adjusted_failure
  
  # result
  ((ff * (1 - h_tP) + gamma) * exp(-lambda * delta_t) + ff * h_tP * exp(-(gamma + lambda + ff) * delta_t)) / (ff + gamma)
  
}

leaky_one_minus_bigF <- function (delta_t, h_tP, gamma, lambda, ff, adjusted_failure) {
  
  # define approach-dependent parameters
  ff <- ff * adjusted_failure
  
  # result
  ((ff * (1 - h_tP) + gamma) * exp(-lambda * delta_t) + ff * h_tP * exp(-(gamma + lambda + ff) * delta_t)) / (ff + gamma)
  
}


# I = 1 - (H + S)
bigF <- function (delta_t, h_tP, AON, gamma, lambda, ff, adjusted_failure) {
  
  if (delta_t < 1) {
    
    0
    
  } else {
    
    # define approach-dependent parameters
    if (AON > 0) {
      h_tP <- h_tP * adjusted_failure
    } else {
      ff   <- ff   * adjusted_failure
    }
    
    
    # result
    1 - ((ff * (1 - h_tP) + gamma) * exp(-lambda * delta_t) + ff * h_tP * exp(-(gamma + lambda + ff) * delta_t)) / (ff + gamma)
    
  }
  
}

AON_bigF <- function (delta_t, h_tP, gamma, lambda, ff, adjusted_failure) {
  
  # define approach-dependent parameters
  h_tP <- h_tP * adjusted_failure
  
  # result
  1 - ((ff * (1 - h_tP) + gamma) * exp(-lambda * delta_t) + ff * h_tP * exp(-(gamma + lambda + ff) * delta_t)) / (ff + gamma)
  
}

leaky_bigF <- function (delta_t, h_tP, gamma, lambda, ff, adjusted_failure) {
  
  # define approach-dependent parameters
  ff <- ff * adjusted_failure
  
  # result
  1 - ((ff * (1 - h_tP) + gamma) * exp(-lambda * delta_t) + ff * h_tP * exp(-(gamma + lambda + ff) * delta_t)) / (ff + gamma)
  
}


# F(t) - F(t-1)
diff_bigF <- function (delta_t, h_tP, AON, gamma, lambda, ff, adjusted_failure) {
  
  if (delta_t < 1) {
    
    0
    
  } else {
    
    # define approach-dependent parameters
    if (AON > 0) {
      h_tP <- h_tP * adjusted_failure
    } else {
      ff   <- ff   * adjusted_failure
    }
    
    
    # result
    ((ff * (1 - h_tP) + gamma) * exp(-lambda * delta_t) * (exp(lambda) - 1) + ff * h_tP * exp(-(gamma + lambda + ff) * delta_t) * (exp(gamma + lambda + ff) - 1)) / (ff + gamma)
    
  }
  
}

AON_diff_bigF <- function (delta_t, h_tP, gamma, lambda, ff, adjusted_failure) {
  
  # define approach-dependent parameters
  h_tP <- h_tP * adjusted_failure
  
  # result
  ((ff * (1 - h_tP) + gamma) * exp(-lambda * delta_t) * (exp(lambda) - 1) + ff * h_tP * exp(-(gamma + lambda + ff) * delta_t) * (exp(gamma + lambda + ff) - 1)) / (ff + gamma)
  
}

leaky_diff_bigF <- function (delta_t, h_tP, gamma, lambda, ff, adjusted_failure) {
  
  # define approach-dependent parameters
  ff <- ff * adjusted_failure
  
  # result
  ((ff * (1 - h_tP) + gamma) * exp(-lambda * delta_t) * (exp(lambda) - 1) + ff * h_tP * exp(-(gamma + lambda + ff) * delta_t) * (exp(gamma + lambda + ff) - 1)) / (ff + gamma)
  
}



# - (H + S)' / (H + S)
hazard_rate <- function (delta_t, h_tP, AON, gamma, lambda, ff, adjusted_failure) {
  
  if (delta_t < 0) {
    
    0
    
  } else {
    
    # define approach-dependent parameters
    if (AON > 0) {
      h_tP <- h_tP * adjusted_failure
    } else {
      ff   <- ff   * adjusted_failure
    }
    
    
    # result
    lambda + ( (ff + gamma) * ff * h_tP ) / ( (ff * (1 - h_tP) + gamma) * exp((gamma + ff) * delta_t) + ff * h_tP )
    
  }
  
}

