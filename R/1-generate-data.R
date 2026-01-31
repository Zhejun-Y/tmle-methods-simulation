
# Function to generate data
generateData <- function(n) {
  
  error <- rnorm(n,mean=0,sd=1)
  
  # two covariates w1 and w2
  w1 <- rnorm(n,mean=0,sd=1)
  w2 <- w1 + error
  
  # binary treatment A generated based on w1 
  A<-rbinom(n,size=1,prob=plogis(w1))
  
  # potential outcome Y(1) and Y(0) generated based on w2 and A
  Y.1<-rbinom(n,size=1,prob=plogis(w2+1))
  Y.0<-rbinom(n,size=1,prob=plogis(w2+0))
  
  # ensure Consistency
  Y<-Y.1*A + Y.0*(1-A)
  
  # return data.frame
  data.frame(w1,w2,A,Y,Y.1,Y.0)
}

