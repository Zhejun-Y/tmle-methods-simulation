
# Function to compute true ATE
true_ATE_compute <- function(n){
  
  set.seed(123)
  true_data <- generateData(n)
  
  return(mean(true_data$Y.1-true_data$Y.0))
}