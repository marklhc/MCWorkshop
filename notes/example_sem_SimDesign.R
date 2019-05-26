library(tidyverse)
library(SimDesign)
library(mnormt)
library(lavaan)

# Design factors:
DESIGNFACTOR <- expand.grid(
  N = c(50, 100, 200), 
  phi22 = c(0.1, 0.5), 
  alpha2 = c(0, 0.5)
)
# Add condition number:
DESIGNFACTOR <- rowid_to_column(DESIGNFACTOR, "cond")
DESIGNFACTOR

# Function for generating data:
gen_lgm_data <- function(condition, fixed_objects = NULL) {
  N <- condition$N
  phi22 <- condition$phi22
  alpha2 <- condition$alpha2
  alpha <- c(1, alpha2)
  Phi <- matrix(c(1, phi22 / 2,
                  phi22 / 2, phi22), nrow = 2)
  Lambda <- cbind(c(1, 1, 1, 1), 
                  c(0, 1, 2, 3))
  Theta <- diag(0.5, nrow = 4)
  # Generate latent factor scores
  eta <- rmnorm(N, mean = alpha, varcov = Phi)
  # Generate residuals:
  e <- rmnorm(N, varcov = Theta)
  # Compute outcome scores
  y <- tcrossprod(eta, Lambda) + e
  colnames(y) <- paste0("y", 1:4)
  # Make it a data frame
  as.data.frame(y)
}
# Test: generate data
# test_df <- gen_lgm_data(DESIGNFACTOR[1, ])

# lavaan syntax (to be passed through `fixed_objects`)
# True model
M1 <- 'i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4
       s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4
       i ~~ s'
# Model without random slopes
M2 <- 'i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4
       s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4
       s ~~ 0*i + 0*s'

# Analysis function
run_lgm <- function(condition, dat, fixed_objects = NULL) {
  # Run model 1
  m1_fit <- growth(fixed_objects$M1, data = dat)
  # Run model 2
  m2_fit <- growth(fixed_objects$M2, data = dat)
  # Extract parameter estimates and standard errors
  ret <- c(coef(m1_fit)["s~1"], 
           sqrt(vcov(m1_fit)["s~1", "s~1"]),
           coef(m2_fit)["s~1"], 
           sqrt(vcov(m2_fit)["s~1", "s~1"]))
  names(ret) <- c("m1_est", "m1_se", "m2_est", "m2_se")
  ret
}
# Test: run analyses for simulated data
# run_lgm(dat = test_df)

# Helper function for computing relative SE bias
rse_bias <- function(est_se, est) {
  est_se <- as.matrix(est_se)
  est <- as.matrix(est)
  est_se <- colMeans(est_se)
  emp_sd <- apply(est, 2L, sd)
  est_se / emp_sd - 1
}

# Evaluative function
evaluate_lgm <- function(condition, results, fixed_objects = NULL) {
  alpha2 <- condition$alpha2
  c(bias = bias(results[ , c(1, 3)], parameter = alpha2), 
    std_bias = bias(results[ , c(1, 3)], parameter = alpha2, 
                    type = "standardized"), 
    rmse = RMSE(results[ , c(1, 3)], parameter = alpha2), 
    rse_bias = rse_bias(results[ , c(2, 4)], results[ , c(1, 3)])
  )
}

# Put all together
sim_results <- runSimulation(DESIGNFACTOR, 
                             replications = 500, 
                             generate = gen_lgm_data, 
                             analyse = run_lgm, 
                             summarise = evaluate_lgm, 
                             fixed_objects = 
                               list(M1 = M1, M2 = M2))

# With parallel processing (4 cores in this example; uncomment the following to
# run)
# sim_results <- runSimulation(DESIGNFACTOR, 
#                              replications = 500, 
#                              generate = gen_lgm_data, 
#                              analyse = run_lgm, 
#                              summarise = evaluate_lgm, 
#                              packages = c("mnormt", "lavaan"), 
#                              parallel = TRUE, 
#                              ncores = 4L, 
#                              fixed_objects = 
#                                list(M1 = M1, M2 = M2))