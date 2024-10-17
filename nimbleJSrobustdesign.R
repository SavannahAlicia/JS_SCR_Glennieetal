
#nonspatial
phipenttp_popan <- function() {
  # Priors and constraints
  for (i in 1:M){
    for (t in 1:(n.occasions-1)){
      phi[i,t] <- mean.phi
    } #t
    for (t in 1:n.occasions){
      p[i,t] <- mean.p
    } #t
  } #i
  mean.phi ~ dunif(0, 1) # Prior for mean survival
  mean.p ~ dunif(0, 1) # Prior for mean capture
  psi ~ dunif(0, 1) # Prior for inclusion probability
  # Dirichlet prior for entry probabilities
  for (t in 1:n.occasions){
    beta[t] ~ dgamma(1, 1)
    b[t] <- beta[t] / sum(beta[1:n.occasions])
  }
  # Convert entry probs to conditional entry probs
  nu[1] <- b[1]
  for (t in 2:n.occasions){
    nu[t] <- b[t] / (1 - sum(b[1:(t-1)]))
  } #t
  # Likelihood
  for (i in 1:M){
    # First occasion
    # State process
    w[i] ~ dbern(psi) # Draw latent inclusion
    z[i,1] ~ dbern(nu[1])
    # Observation process
    mu1[i] <- z[i,1] * p[i,1] * w[i]
    y[i,1] ~ dbern(mu1[i])
    # Subsequent occasions
    for (t in 2:n.occasions){
      # State process
      q[i,t-1] <- 1 - z[i,t-1]
      mu2[i,t] <- phi[i,t-1] * z[i,t-1] + nu[t] * prod(q[i,1:(t-1)])
      z[i,t] ~ dbern(mu2[i,t])
      # Observation process
      mu3[i,t] <- z[i,t] * p[i,t] * w[i]
      y[i,t] ~ dbern(mu3[i,t])
    } #t
  } #i
  # Calculate derived population parameters
  for (i in 1:M){
    for (t in 1:n.occasions){
      u[i,t] <- z[i,t] * w[i] # Deflated latent state (u)
    }
  }
  for (i in 1:M){
    recruit[i,1] <- u[i,1]
    for (t in 2:n.occasions){
      recruit[i,t] <- (1 - u[i,t-1]) * u[i,t]
    } #t
  } #i
  for (t in 1:n.occasions){
    N[t] <- sum(u[1:M,t]) # Actual population size
    B[t] <- sum(recruit[1:M,t]) # Number of entries
  } #t
  for (i in 1:M){
    Nind[i] <- sum(u[i,1:n.occasions])
    Nalive[i] <- 1 - equals(Nind[i], 0)
  } #i
  Nsuper <- sum(Nalive[]) # Superpopulation size
}

chs_noneuclidean_extraprims2000 <- readRDS("~/Documents/UniStAndrews/Dolphins/Charleston/data/onison/all_occasions/model_objs/chs_noneuclidean_extraprims2000.Rds")
ch <- chs_noneuclidean_extraprims2000$capthist()
traps <- chs_noneuclidean_extraprims2000$traps()
distmat <- chs_noneuclidean_extraprims2000$userdistmat()

#spatial
phipenttp_popan <- function() {
  # Priors and constraints
  for (i in 1:M){
    for (tp in 1:(n.prim.occasions-1)){
      phi[i,tp] <- mean.phi
    } #t
  } #i
  mean.phi ~ dunif(0, 1) # Prior for mean survival
  mean.lambda0 ~ dunif(0, 1) # Prior for mean capture
  mean.sigma ~ dunif(1, 10000)
  psi ~ dunif(0, 1) # Prior for inclusion probability
  # Dirichlet prior for entry probabilities
  for (t in 1:n.occasions){
    beta[t] ~ dgamma(1, 1)
    b[t] <- beta[t] / sum(beta[1:n.occasions])
  }
  # Convert entry probs to conditional entry probs
  nu[1] <- b[1]
  for (t in 2:n.occasions){
    nu[t] <- b[t] / (1 - sum(b[1:(t-1)]))
  } #t
  ##Density and home range centers???
  
  # Likelihood
  for (i in 1:M){
    # First occasion
    # State process
    w[i] ~ dbern(psi) # how to think of transition per season???
    z[i,1] ~ dbern(nu[1])
    # Observation process
    #Detection at multi-catch trap
    hus <- array(0, dim = J)
    exphus <- array(0, dim = J)
    for(j in 1:J){
      ##need to figure out how to calculate noneuclidean distances
      dij <- distmat[i, j] #?? doesn't exist yet
      hij <- mean.lambda0 * exp(-(dij^2)/(2*mean.sigma^2))
      ujt <- usage(traps)[j, 1]
      hus[j] <- hij * ujt
      exphus[j] <- exp(-hij * ujt)
    }
    Jprobs <- array(0, dim = J+1)
    for(j in 1:J){
      Jprobs[j] <- hus[j]/sum(hus)*(1 - prod(exphus))
    }
    Jprobs[J+1] <- prod(exphus) #J+1 indexed outcome is no detection
    p[i,1] <- dmultinom(Jprobs) #how to make sure the probability corresponds to the right index???
    mu1[i] <- z[i,1] * p[i,1] * w[i]
    y[i,1] ~ #???
    
    # Subsequent occasions
    for (t in 2:n.occasions){
      # State process
      q[i,t-1] <- 1 - z[i,t-1]
      mu2[i,t] <- phi[i,t-1] * z[i,t-1] + nu[t] * prod(q[i,1:(t-1)])
      z[i,t] ~ dbern(mu2[i,t])
      # Observation process
      
      
      mu3[i,t] <- z[i,t] * p[i,t] * w[i]
      y[i,t] ~  #multinomial, either detected at one trap j or not detected
    } #t
  } #i
  # Calculate derived population parameters
  for (i in 1:M){
    for (t in 1:n.occasions){
      u[i,t] <- z[i,t] * w[i] # Deflated latent state (u)
    }
  }
  for (i in 1:M){
    recruit[i,1] <- u[i,1]
    for (t in 2:n.occasions){
      recruit[i,t] <- (1 - u[i,t-1]) * u[i,t]
    } #t
  } #i
  for (t in 1:n.occasions){
    N[t] <- sum(u[1:M,t]) # Actual population size
    B[t] <- sum(recruit[1:M,t]) # Number of entries
  } #t
  for (i in 1:M){
    Nind[i] <- sum(u[i,1:n.occasions])
    Nalive[i] <- 1 - equals(Nind[i], 0)
  } #i
  Nsuper <- sum(Nalive[]) # Superpopulation size
}


