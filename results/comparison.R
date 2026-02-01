
# Case: Using TMLE if it is necessary

sample_sizes <- c(50, 100, 200, 500, 1000, 1500, 2000)

var_tmle_misspec <- numeric(length(sample_sizes))
var_mle_misspec <- numeric(length(sample_sizes))

bias_tmle_misspec <- numeric(length(sample_sizes))
bias_mle_misspec <- numeric(length(sample_sizes))

ci_tmle_misspec <- numeric(length(sample_sizes))

index <- 1:length(sample_sizes)

for (i in index) {
  
  # TMLE with correct models
  res_tmle_misspec <- tmle_simulation(
    n_sim = 1000,          
    n = sample_sizes[i],  
    outcome_model = Y ~ A,   # misspecified outcome
    exposure_model = A ~ w1  # correctly specified exposure
  )
  var_tmle_misspec[i] <- res_tmle_misspec$Var_eff
  bias_tmle_misspec[i] <- res_tmle_misspec$mean_bias
  ci_tmle_misspec[i] <- res_tmle_misspec$mean_coverage
  
  # MLE (outcome regression, no targeting step)
  res_mle_misspec <- tmle_simulation(
    n_sim = 1000,          
    n = sample_sizes[i],  
    outcome_model = Y ~ A,   # misspecified outcome
    exposure_model = A ~ 1 # MLE 
  )
  var_mle_misspec[i] <- res_mle_misspec$Var_eff
  bias_mle_misspec[i] <- res_mle_misspec$mean_bias
  
  
}
