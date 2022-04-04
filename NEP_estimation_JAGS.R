model{
  # Observation level priors
  # gamma is the conjugate for the precision of the normal distribution
  # is defined for positive real numbers
  precision.delta ~ dgamma(0.001, 0.001)

  # Group level priors
  mu.alpha ~ dunif(prior.min, prior.max)
  precision.alpha ~ dgamma(0.001, 0.001)

  for (i in 1:n){
      precision.y[i]<-1/(sigma[i]^2)
      
      # Observations drawn from distribution of study-specific mean and variance
      y[i]~dnorm(delta[i], precision.y[i]) 
      
      # Study specific mean drawn from distribution of system-specific mean 
      delta[i]~dnorm(mu[i], precision.delta) 
      mu[i] <- alpha[group[i]] 
  }
  
  for(j in 1:n.group){
      # System specific mean drawn from a normal distribution with vague priors
      alpha[j]~dnorm(mu.alpha, precision.alpha)
  }
}
