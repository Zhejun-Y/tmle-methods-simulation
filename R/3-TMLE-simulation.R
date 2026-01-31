
tmle_simulation <- function(n_sim, n, outcome_model, exposure_model) {
  
  # Store results
  ATE_estimates <- numeric(n_sim)
  Std_Errors <- numeric(n_sim)
  Biases <- numeric(n_sim)
  coverage <- numeric(n_sim)
  
  # Compute True ATE 
  true_ATE <- true_ATE_compute(1e6)
  
  # TMLE Algorithm
  for (i in 1:n_sim) {
    set.seed(6001 + i)
    ObsData <- generateData(n)
    
    # Step 1: Fit outcome model
    Q_model <- glm(outcome_model, family = binomial, data = ObsData)
    QAw <- predict(Q_model, type = "response")
    Q1w <- predict(Q_model, newdata = transform(ObsData, A = 1), type = "response")
    Q0w <- predict(Q_model, newdata = transform(ObsData, A = 0), type = "response")
    
    # Step 2: Fit propensity score model
    g_model <- glm(exposure_model, family = binomial, data = ObsData)
    gw <- predict(g_model, type = "response")
    
    # Step 3: Updating step
    H1w <- ObsData$A / gw
    H0w <- (1 - ObsData$A) / (1 - gw)
    epsilon <- coef(glm(ObsData$Y ~ -1 + H0w + H1w + 
                          offset(qlogis(pmin(pmax(QAw, 1e-8), 1 - 1e-8))), family = binomial))
    
    Q1w_star <- plogis(qlogis(pmin(pmax(Q1w, 1e-8), 1 - 1e-8)) + epsilon["H1w"] / gw)
    Q0w_star <- plogis(qlogis(pmin(pmax(Q0w, 1e-8), 1 - 1e-8)) + epsilon["H0w"] / (1 - gw))
    
    # TMLE ATE
    EY1_tmle <- mean(Q1w_star)
    EY0_tmle <- mean(Q0w_star)
    ATE_tmle <- EY1_tmle - EY0_tmle
    
    # Influence curve and SE
    D1 <- ObsData$A / gw * (ObsData$Y - Q1w_star) + Q1w_star - EY1_tmle
    D0 <- (1 - ObsData$A) / (1 - gw) * (ObsData$Y - Q0w_star) + Q0w_star - EY0_tmle
    ATE_IC <- D1 - D0
    ATE_var <- var(ATE_IC) / n
    Std_Err <- sqrt(ATE_var)
    
    # Save results
    ATE_estimates[i] <- ATE_tmle
    Std_Errors[i] <- Std_Err
    Biases[i] <- ATE_tmle - true_ATE
    
    CI_low <- ATE_tmle - 1.96 * Std_Err
    CI_high <- ATE_tmle + 1.96 * Std_Err
    coverage[i] <- (true_ATE >= CI_low) & (true_ATE <= CI_high)
    
  }
  
  # Summary
  
  
  result <- list(
    true_ATE = true_ATE,
    mean_tmle = mean(ATE_estimates),
    mean_bias = mean(Biases),
    mean_est_SE = mean(Std_Errors),
    empirical_SE = sd(ATE_estimates),
    mean_coverage = mean(coverage),
    
    ATE_Estimates = ATE_estimates,
    Var_eff = (sd(ATE_estimates))^2
    
  )
  
  return(result)
}
