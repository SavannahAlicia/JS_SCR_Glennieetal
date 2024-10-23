library(nimble)
chs_noneuclidean_extraprims2000 <- readRDS("~/Documents/UniStAndrews/Dolphins/Charleston/data/onison/all_occasions/model_objs/chs_noneuclidean_extraprims2000.Rds")
ch <- chs_noneuclidean_extraprims2000$capthist()
traps <- chs_noneuclidean_extraprims2000$traps()
distmat <- chs_noneuclidean_extraprims2000$userdistmat()
mesh <- chs_noneuclidean_extraprims2000$mesh()
primary <- m$data()$primary()
usage.traps <- usage(traps)

#create capture history that has a row for no detections per occasion
ch0 <- array(data = NA, dim = c(dim(ch)[1:2], dim(ch)[3]+1))
ch0[1:dim(ch)[1], 1:dim(ch)[2], 1:dim(ch)[3]] <- ch
ch0[1:dim(ch)[1], 1:dim(ch)[2], (dim(ch)[3]+1)] <- ifelse(apply(ch, c(1,2), sum)>0, 0, 1)
dimnames(ch0) <- list(dimnames(ch)[[1]], dimnames(ch)[[2]], c(dimnames(ch)[[3]], "no_detection"))

##testing 
rho_ua_Intercept <- apply(as.array(1:77), 1, function(x){m$state()$trm(k = x)[2,1]})
rho_au_Intercept <- apply(as.array(1:77), 1, function(x){m$state()$trm(k = x)[1,2]})
beta <- m$get_par("beta", m = 1, s = 1)
phi <- m$get_par("phi", m = 1, s = 1)
delta <- m$state()$delta()[1]
n.prim.occasions <- length(beta)
dt <- diff(m$data()$time())


#need functions sum, prod, which, expm::expm, usage,
expmexpm <- function(x){
  expm::expm(x)
}
Rexpm <- nimbleRcall(function(x = double(2)){}, Rfun = 'expmexpm',
                     returnType = double(2))
Rt <- nimbleRcall(function(x = double(2)){}, Rfun = 't',
                  returnType = double(2))

#spatial
JS_SCR <- nimbleCode({
  
  # Priors and constraints
  lambda0 ~ dunif(0, 1) # tilda takes left thing and puts it in x spot
  sigma ~ dunif(0, 10000)
  sigma2 <- sigma * sigma
  # 2 alive states, available and unavailable
  delta ~ dunif(0,1) # determines relative recruitment into alive state 
  
  #survival and transition rates
    for (tp in 1:(n.prim.occasions-1)){
      phi[tp] ~ dunif(0, 1)
      rho_au_Intercept[tp] ~ dunif(0, 20) # Prior for transition from available to unavailable Intercept (season2summerfall) (so I can add other season later)
      rho_ua_Intercept[tp] ~ dunif(0, 20)#Prior for transition rate from unavailable to available Intercept (season2summerfall)
    } #tp

  # Dirichlet prior for entry probabilities
  ##Do I need a prior for alpha?
  alpha[1:n.prim.occasions] <- rep(1, n.prim.occasions) 
  beta[1:n.prim.occasions] ~ ddirch(alpha[1:n.prim.occasions])
  
  # Convert entry probs to conditional entry probs
  b[1] <- beta[1]
  for (t in 2:n.prim.occasions){
    b[t] <- beta[t]/prod(1 - b[1:(t-1)]) #is this vectorizing correctly?
  } #t

  # For prob available given alive, I need state transitions, including available and unavalaible states
    # state membership probability
 # stateprob <- array(data = 0, dim = c(n.prim.occasions,4))
    #first primary is entry probability and delta
  stateprob[1,1:4] <- c(1 - b[1], b[1]*delta, b[1]*(1-delta), 0)
 # Q <- array(data = 0, dim = c(4,4,(n.prim.occasions - 1))) #transition rates (per year, since dt is time difference in years between prims)
#  G <- array(data = 0, dim = c(4,4,(n.prim.occasions - 1))) #transition probabilities
  for (tp in 1:(n.prim.occasions - 1)){
    #get rid of array and t 
    Q[1:4,1:4,tp] <- c( #note recruitment index starts at 1, but transitions happen between occasions 
      log(1 - b[(tp+1)])/dt[tp], 0,0,0,
      -log(1 - b[(tp+1)]*delta)/dt[tp], -rho_au_Intercept[tp] + log(phi[tp]), rho_ua_Intercept[tp],          0,
      -log(1 - b[(tp+1)]*(1 - delta))/dt[tp],  rho_au_Intercept[tp],      -rho_ua_Intercept[tp] + log(phi[tp]),  0,
      0, -log(phi[tp]), -log(phi[tp]), 0
    )
      #have to fill matrix by column, not row
    #log(1 - b[(tp+1)])/dt[tp],   -log(1 - b[(tp+1)]*delta)/dt[tp],       -log(1 - b[(tp+1)]*(1 - delta))/dt[tp],  0,
    #0,                           -rho_au_Intercept[tp] + log(phi[tp]),   rho_au_Intercept[tp],                    -log(phi[tp]),
    #0,                           rho_ua_Intercept[tp],                   -rho_ua_Intercept[tp] + log(phi[tp]),    -log(phi[tp]),
    #0, 0, 0, 0
    #)
  G[1, 1:4, tp] <- c(1-b[(tp+1)], b[(tp+1)]*delta, b[(tp+1)]*(1-delta), 0) #from unborn
  G[2:4, 1, tp] <- c(0,0,0) #cannot transition back to unborn
  #set expm as R function??
  Qd[1:3,1:3,tp] <- (Q[2:4,2:4,tp] * dt[tp])

  G[2:4,2:4,tp] <- Rexpm(Qd[1:3,1:3,tp]) #matrix exponential for competing states
#  stateprob[(tp+1),1:4] <- stateprob[(tp),1:4] %*% G[1:4,1:4,tp] #membership probability 
  }
  # Psi
#  for(tp in 1:(n.prim.occasions - 1)){
#    psi[tp] <- stateprob[tp,2]/sum(stateprob[tp,2:3]) #probability of being available given alive
#    alive[tp] <- stateprob[tp,2] + stateprob[tp,3]
#  }
}) 
  # Likelihood
  for (i in 1:M){
    Smesh[i] ~ dunif(1, nrow(mesh))
    for (t in 1:n.prim.occasions){
      # State process
      w[i,t] ~ dbern(psi[t]) #available given alive
      z[i,t] ~ dbern(alive[t]) #alive in primary
      # Observation process
      sec_in_prim <- which(primary == t)
      #Detection at multi-catch trap
      for (s in sec_in_prim){
        hus <- array(0, dim = J)
        exphus <- array(0, dim = J)
        for(j in 1:J){
          dij <- distmat[j, Smesh[i]]
          hij <- lambda0 * exp(-(dij * dij)/(2*sigma2))
          ujt <- usage.traps[j, t]
          hus[j] <- hij * ujt
          exphus[j] <- exp(-hij * ujt)
          Jprobs[i,s,j] <- hus[j]/sum(hus)*(1 - prod(exphus))
        } #j
        Jprobs[i,s,J+1] <- prod(exphus) #J+1 indexed outcome is no detection
        mu[i,s,1:(J+1)] <- Jprobs[i,s,1:(J+1)] * z[i,t] * w[i,t]
        y[i,s,1:(J+1)] ~ dmultinom(mu[i,s,1:(J+1)]) #multinomial, either detected at one trap j or not detected
      } #s (secondary within primary)
    } #t (primary)
  } #i (individual)
  # Calculate derived population parameters
  for (i in 1:M){
    for (t in 1:n.occasions){
      u[i,t] <- z[i,t] * w[i,t] # Deflated latent state (u)
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
    Nalive[i] <- 1 - Nind[i][which(Nind[i] == 0)]
  } #i
  Nsuper <- sum(Nalive) # Superpopulation size
})

ch13 <- ch0[,primary %in% 1:3,]
ch13 <- ch13[apply(ch13[,,-dim(ch13)[3]], 1, sum) > 0,,]
usage.traps13 <- usage.traps[,primary %in% 1:3]
primary13 <- primary[primary %in% 1:3]
dt13 <-  dt[1:3]

data <- list(y = ch13, 
             mesh = as.matrix(mesh),
             usage.traps = usage.traps13,
             distmat = distmat,
             primary = primary13,
             dt = dt13)

constants <- list(n.prim.occasions = 3,
                  J = 91,
                  M = 700)

inits <- list(
  phi = m$get_par("phi", m = 1, j = 1, s = 1)[1:2], 
  lambda0 = m$get_par("lambda0", m = 1, j = 1, s= 1, k = 1),
  sigma = m$get_par("sigma", m = 1, j = 1, s= 1, k = 1),
  delta = m$state()$delta()[1],
  beta = m$get_par("beta", m = 1, j = 1, s = 1)[1:3], #rdirch(n.prim.occasions, 1), #rdirch for more than one n doesn't work
  rho_ua_Intercept = apply(as.array(1:2), 1, function(x){m$state()$trm(k = (x+1))[2,1]}), #can't remember if its x or x+1
  rho_au_Intercept = apply(as.array(1:2), 1, function(x){m$state()$trm(k = (x+1))[1,2]})
)

Rmodel <- nimbleModel(code = JS_SCR, 
                      constants = constants, 
                      data = data,
                      init = inits
                      )
