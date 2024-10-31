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
  stateprob[(tp+1),1:4] <- stateprob[(tp),1:4] %*% G[1:4,1:4,tp] #membership probability 
  }
   #Psi
  for(t in 1:n.prim.occasions){
        psi[t] <- stateprob[t,2]/sum(stateprob[t,2:3]) #probability of being available given alive
        alive[t] <- stateprob[t,2] + stateprob[t,3] #probability alive
        omega[t] ~ dbeta(5,9)
  }
  
 
  # Likelihood
  for (i in 1:M){
    #sample x and y space, then check if AC is within 
    S[i, 1] ~ dunif(1, upperlimitx)    # x-coord of activity centers
    S[i, 2] ~ dunif(1, upperlimity)    # y coord of activity centers
    hab[i] <- habMat[trunc(S[i, 1]), trunc(S[i, 2])] #for now, just uniformly sampling available mesh, will need to change
    ones[i] ~ dbern(hab[i]) #the ones trick
    Smesh[i] <- whichmesh[trunc(S[i, 1]), trunc(S[i, 2])]
    real[i] ~ dbern(omega[i])
    for (t in 1:n.prim.occasions){
      # State process
      w[i,t] ~ dbern(psi[t]) #available given alive
      z[i,t] ~ dbern(alive[t]) #alive in primary
      # Observation process
      ###can't replace or have flexible size variables
    }#t (primary)
      #Detection at multi-catch trap
      for (s in 1:n.sec.occasions){
        #hus <- array(0, dim = J)
        #exphus <- array(0, dim = J)
        for(j in 1:J){
          #again, this will change when hrc model changes
          hus[j,s,i] <- lambda0 * exp(-( distmat[j, Smesh[i]] *  distmat[j, Smesh[i]])/(2*sigma2)) * usage.traps[j, s]
          exphus[j,s,i] <- exp(-hus[j,s,i])
          Jprobs[i,s,j] <- hus[j,s,i]/sum(hus[1:J,s,i])*(1 - prod(exphus[1:J,s,i]))
        } #j
        Jprobs[i,s,J+1] <- prod(exphus[1:J,s,i]) #J+1 indexed outcome is no detection
        mu[i,s,1:(J+1)] <- Jprobs[i,s,1:(J+1)] * z[i,primary[s]] * w[i,primary[s]] * real[i] #is this right for data augmentation?
        #I think the trap of detection is multinomial, but then the alive, present, and real are bernoulli... do I need a custom distribution?
        y[i,s,1:(J+1)] ~ dmultinom(mu[i,s,1:(J+1)]) #multinomial, either detected at one trap j or not detected
      } #s (secondarys)
  } #i (individual)
  #Derived parameters

  
})
  
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

n.prim.occasions <- 3
n.fake.inds <- 500
M <- dim(ch13)[1] + n.fake.inds
n.sec.occasions <- sum(primary %in% 1:n.prim.occasions)
#subset capture history
ch13 <- ch0[,1:n.sec.occasions,]
ch13 <- ch13[apply(ch13[,,-dim(ch13)[3]], 1, sum) > 0,,] #remove 0 chs
#data augmentation
nodets <- matrix(0, nrow = n.sec.occasions, ncol = dim(ch13)[3])
nodets[,dim(ch13)[3]] <- 1
aug <- aperm(replicate(n.fake.inds, nodets), c(3,1,2))
augch13 <- abind::abind(ch13, aug, along = 1)
dimnames(augch13)[2:3] <- dimnames(ch13)[2:3]
real <- c(rep(1, dim(ch13)[1]), rep(NA, n.fake.inds))
#subset trap usage
usage.traps13 <- usage.traps[,1:n.sec.occasions]
#subset primary code
primary13 <- primary[1:n.sec.occasions]
#subset primary dt
dt13 <-  dt[1:n.prim.occasions]
#convert mesh to raster matrix of 1 and 
habMat <- rasterFromXYZ(cbind(mesh, rep(1, nrow(mesh))))
whichmesh <- rasterFromXYZ(cbind(mesh, 1:nrow(mesh)))

data <- list(y = augch13, 
             real = real,
             habMat = as.matrix(habMat),
             whichmesh = as.matrix(whichmesh),
             usage.traps = usage.traps13,
             distmat = distmat,
             dt = dt13,
             upperlimitx = ncol(habMat),
             upperlimity = nrow(habMat),
             ones = rep(1, M))

constants <- list(n.prim.occasions = 3,
                  n.sec.occasions = 18,
                  J = 91,
                  M = M,
                  primary = primary13,
                  n.mesh = 143
                  )

inits <- list(
  phi = m$get_par("phi", m = 1, j = 1, s = 1)[1:2], 
  lambda0 = m$get_par("lambda0", m = 1, j = 1, s= 1, k = 1),
  sigma = m$get_par("sigma", m = 1, j = 1, s= 1, k = 1),
  delta = m$state()$delta()[1],
  beta = m$get_par("beta", m = 1, j = 1, s = 1)[1:3], #rdirch(n.prim.occasions, 1), #rdirch for more than one n doesn't work
  rho_ua_Intercept = apply(as.array(1:2), 1, function(x){m$state()$trm(k = (x))[2,1]}),
  rho_au_Intercept = apply(as.array(1:2), 1, function(x){m$state()$trm(k = (x))[1,2]})
)

Rmodel <- nimbleModel(code = JS_SCR, 
                      constants = constants, 
                      data = data,
                      init = inits
                      )
