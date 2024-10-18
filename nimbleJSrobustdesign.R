
chs_noneuclidean_extraprims2000 <- readRDS("~/Documents/UniStAndrews/Dolphins/Charleston/data/onison/all_occasions/model_objs/chs_noneuclidean_extraprims2000.Rds")
ch <- chs_noneuclidean_extraprims2000$capthist()
traps <- chs_noneuclidean_extraprims2000$traps()
distmat <- chs_noneuclidean_extraprims2000$userdistmat()

##testing 
rho_ua_Intercept <- apply(as.array(1:77), 1, function(x){m$state()$trm(k = x)[2,1]})
rho_au_Intercept <- apply(as.array(1:77), 1, function(x){m$state()$trm(k = x)[1,2]})
beta <- m$get_par("beta", m = 1, s = 1)
phi <- m$get_par("phi", m = 1, s = 1)
delta <- m$state()$delta()[1]
n.prim.occasions <- length(beta)
dt <- diff(m$data()$time())


#spatial
phipenttp_popan <- function() {
  
  # Priors and constraints
  phi ~ dunif(0, 1) # Prior for survival
  lambda0 ~ dunif(0, 1) # 
  sigma ~ dunif(1, 10000)
  # 2 alive states, available and unavailable
  delta ~ dunif(0,1) # determines relative recruitment into alive state 
  rho_au_Intercept ~ dunif(0, 20) # Prior for transition from available to unavailable Intercept (season2summerfall) (so I can add other season later)
  rho_ua_Intercept ~ dunif(0, 20) #Prior for transition rate from unavailable to available Intercept (season2summerfall)
  
    #survival and transition rates
    for (tp in 1:(n.prim.occasions-1)){
      phi[tp] <- phi
      rho_au_Intercept[tp] <- rho_au_Intercept
      rho_ua_Intercept[tp] <- rho_ua_Intercept
    } #t

  # Dirichlet prior for entry probabilities
  for (t in 1:n.prim.occasions){
    beta1[t] ~ dgamma(1, 1)
    beta[t] <- beta1[t] / sum(beta1[1:n.occasions])
  }
  # Convert entry probs to conditional entry probs
  b[1] <- beta[1]
  for (t in 2:n.prim.occasions){
    prod1mb <- 1
    for(tt in 1:(t-1)){
      prod1mb <- prod1mb * (1 - b[tt])
    }
   b[t] <- beta[t] / prod1mb 
  #b[t] <- beta[t]/(1 - sum(beta[1:(t-1)])) #is this the same?? rounding error?
  } #t
  
  # For prob available given alive, I need state transitions, including available and unavalaible states
    # state membership probability
  stateprob <- array(data = 0, dim = c(n.prim.occasions,4))
    #first primary is entry probability and delta
  stateprob[1,] <- array(c(1 - b[1], b[1]*delta, b[1]*(1-delta), 0))
  Q <- array(data = 0, dim = c(4,4,(n.prim.occasions - 1))) #transition rates (per year, since dt is time difference in years between prims)
  G <- array(data = 0, dim = c(4,4,(n.prim.occasions - 1))) #transition probabilities
  for (tp in 1:(n.prim.occasions - 1)){
    Q[,,tp] <- t(array(c( #note recruitment index starts at 1, but transitions happen between occasions 
    log(1 - b[(tp+1)])/dt[tp],   -log(1 - b[(tp+1)]*delta)/dt[tp],       -log(1 - b[(tp+1)]*(1 - delta))/dt[tp],  0,
    0,                           -rho_au_Intercept[tp] + log(phi[tp]),   rho_au_Intercept[tp],                    -log(phi[tp]),
    0,                           rho_ua_Intercept[tp],                   -rho_ua_Intercept[tp] + log(phi[tp]),    -log(phi[tp]),
    0, 0, 0, 0
    ), dim = c(4, 4)))
  G[1, 1:4, tp] <- c(1-b[(tp+1)], b[(tp+1)]*delta, b[(tp+1)]*(1-delta), 0) #from unborn
  G[2:4,2:4,tp] <- as.array(expm::expm((Q[2:4,2:4,tp] * dt[tp]))) #matrix exponential for competing states
  stateprob[(tp+1),] <- stateprob[(tp),] %*% G[,,tp] #membership probability 
  }
  # Psi
  for(tp in 1:(n.prim.occasions - 1)){
    psi[tp] <- stateprob[tp,2]/sum(stateprob[tp,2:3]) #probability of being available given alive
  }
  
  ##Density and home range centers (to do)
  
  # Likelihood
  for (i in 1:M){
    # First occasion
    # State process
    w[i,1] ~ dbern(psi[1]) #available given alive
    z[i,1] ~ dbern(b[1]) #alive in first occasion
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
    ##how to make sure the probability corresponds to the right index???
    p[i,1] <- dmultinom(Jprobs) #probability of capture history for i,k given i alive and available
    mu1[i] <- p[i,1] * z[i,1] * w[i] #capthist and alive and available
    y[i,1] ~ ()#???
    
    # Subsequent occasions
    for (t in 2:n.occasions){
      # State process

      # Observation process
      
      
      mu2[i,t] <- z[i,t] * p[i,t] * w[i,t]
      y[i,t] ~ () #multinomial, either detected at one trap j or not detected
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


