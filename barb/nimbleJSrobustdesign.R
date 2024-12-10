library(nimble)
library(secr)
library(raster)
library(expm)

## Nimble function that calls R matrix eponential from expm package. My base R didnt' have another expm function...
Rexpm <- nimbleRcall(function(x = double(2)){}, Rfun = 'expm',
                     returnType = double(2))

## Main methods function in Nimble that contains all the data to reduce the complexity of the model structure in NIMBLE
## Key is that we are going to put some cleverish caching in here so that we don't have to create too many jprobs in code.
JSguts_nf <- nimbleFunction(
  setup = function(habmat, whichmesh, distmat, dt){
    habmat <- as.matrix(habmat) #just in case stored as df or raster
    habmat[is.na(habmat)] <- 0  ## PVDB: Upated to make sure NA is actually 0.
    whichmesh <- as.matrix(whichmesh)
    J <- nrow(distmat)
    dist2 <- distmat*distmat
    n.sec.occasions <- ncol(trapusage)
    n.prim.occasions <- length(dt) + 1
  },
  # run = function(){},
  methods = list(
    dhab = function(x = double(), #ones data
                     S = double(1), #array of length 2 x and y
                     log = logical(0, default = FALSE)){
      hab <- habmat[trunc(S[1]), trunc(S[2])] #choose from mesh
      if(log) return(log(hab))
      else return(hab)
      returnType(double())
    },
    probdetected = function(lambda = double(), 
                           sigma = double(), 
                           S = double(1)){
      smesh <- whichmesh[trunc(S[2]), trunc(S[1])]  #index of mesh COLUMNS ARE X AND ROWS ARE Y!!!!! *PVDB
      jprobs <- numeric(value = 0, length = J)
      if(is.na(smesh)) return(jprobs) ## PVDB: Made it more robust. Just in case distance isn't calcualted make Jprobs 0.
      
      jprobs <- lambda*exp(-dist2[1:J, smesh]/(2*sigma^2))  ## This could be more efficient. It's where bischoff only computes on nearby traps...
      returnType(double(1))
      return(jprobs)
    },
    transprobs = function(phi = double(1), 
                          b = double(1)){ #note could vectorize to speed up
      G <- array(NA, dim = c(3,3,(n.prim.occasions - 1)))
      for(tp in 1:(n.prim.occasions - 1)){
        G[1:3, 1:3, tp] <- calcGt(phi[1:(n.prim.occasions - 1)], b[1:n.prim.occasions], tp)
      }
      returnType(double(3))
      return(G)
    },
    calcGt = function(phi = double(1), 
                     b = double(1),
                     tp = double()){
      # State membership transition probability 
        #Q is rates
      Q <- matrix(0, 3, 3)
      G <- matrix(0, 3, 3)
      Q[1, 1:3] <- c(NA, NA, 0)
      Q[2, 1:3] <- c(0, log(phi[tp]), -log(phi[tp]))
      Q[3, 1:3] <- c(0, 0, 0)
      #have to fill matrix by column, not row
      G[1, 1:3] <- c(1-b[(tp+1)], b[(tp+1)], 0) #from unborn
      G[2:3, 1] <- c(0,0) #cannot transition back to unborn
      G[2:3,2:3] <- Rexpm((Q[2:3,2:3] * dt[tp])) #matrix exponential for competing states
      returnType(double(2))
      return(G)
    }
  )
)

dcapt_forward_internal <- nimbleFunction(
  setup = function(capthist, trapusage, primary){
    J <- nrow(trapusage)
    nprimary <- max(primary)
    
  }, # Hold data to make model smaller.
  run = function(){},
  methods = list(
    dcapt_marg = function(x = double(), # index of animal 1:M
           real = double(),     # real or not real animal?
           Jprobs = double(1), # probability of detection for each secondary at each trap (single individual)
           G = double(3), # Transition prob matrix
           init = double(0),
           log = logical(0, default=TRUE)){
    ch <- capthist[x,]       
    ## Check if it's not real then prob = 1, unless it was observed, then shit.
    isobs0 <- all(ch == J+1) #all unobserved
    if(real == 0){ # if augmented individual
      ans <- 1 # it must be unobserved
      if(!isobs0){ # if it was observed
        ans <- 0 # that's not possible
      }
      if(log) return(log(ans))
      else return(ans)
    }
    
    ## Now to forward solve.
    ## State 1: unborn
    ## State 2: alive
    ## State 3: dead
    nstates <- dim(G)[1]
    pi <- nimNumeric(value = 0, length = nstates) # prob of state as we go forward
    pobs <- nimNumeric(value = 0, length = nstates) # Pdetect given state
    logL <- 0
    pi[init] <- 1 #init specifies number of state we begin in
        
    for( tp in 1:(nprimary) ){
      idx <- which(primary == tp) #secondaries in primary
      isobs <- any(ch[idx] != (J+1))
      ## Find probability of obs within the primary.
      for( state in 1:nstates ){
        if(state == 1 | state == 3){ #if unborn or dead
          if(isobs) {
            pobs[state] <- 0 #prob 0 of being observed
         } else {
            pobs[state] <- 1 #prob 1 of being unobserved
         }
        } else {
          lp <- 0
          for( j in seq_along(idx) ) {
            occ <- idx[j]
            ## Need to compute detection probs:
            ## Part 1: Not detected but keep on log scale
            keep <- which(trapusage[,occ] > 0)  ## Save a little compute time maybe...
            logpnotdet <- -sum(trapusage[keep,occ]*Jprobs[keep])
            ## If not detected then I'm done. Otherwise calc the prob of being det at just that trap.
            if(ch[occ] == J+1) {
              lp <- lp + logpnotdet
            }else{
              ## If can't be detected, but was then return -infinity.
              if(logpnotdet < 0) {
                lp <- lp + log(Jprobs[ch[occ]]) - log(-logpnotdet) + log(1-exp(logpnotdet)) ## log(pdetect/sum(pdetect) * (1-pnotdet))
              }else{
                lp <- -Inf
              }
            }
          }
          pobs[state] <- exp(lp) #product of probs across secondaries within primary
        }
      }
      pi <- pi * pobs #prob of ch and state
      sumpi <- sum(pi) #summed across states = prob of ch
      logL <- logL + log(sumpi) 
      if (tp < nprimary) { #for all but last primary, times transition probs
        pi <- ((pi %*% G[,,tp])/sumpi)[1,] #
      }
    }
    ans <- sum(pi)
    returnType(double())
    if(log) return(logL)
    else return(exp(logL))
  }
  )
)

#Testing:
# x = double(1), # Observed capture history for one animal
           # primary = double(1), # vector of primary indices for each secondary
           # real = double(),     # real or not real animal?
           # Jprobs = double(2), # probability of detection for each secondary at each trap
           # G = double(2), # Transition prob matrix
           # init = double(0)
# Jprobs <- do.call('rbind', lapply(1:9, FUN = function(x){rdirch(1, rep(1,4))}))
# G <- array(rbind(c(0.3, 0.7, 0), c(0, 0.99, 0.01), c(0, 0, 1)), c(3, 3, 2))
# 
# # debugonce(dcapt_forward)
# dcapt_forward(x=c(4,4,4,1,1,4,4,4,4), primary = c(1,1,1,2,2,3,3,3,3), 
#   real = 1, Jprobs, G, init = 1, log = FALSE)
# Jprobs[,1]*0.99^2
# prod(Jprobs[1:5,1])*0.99*(0.99*prod(Jprobs[6:9,4]) + 0.01)
# 1*0.7*prod(Jprobs[4:5,1])*(0.99*prod(Jprobs[6:9,4]) + 0.01)

#spatial
JS_SCR <- nimbleCode({
  
  # Priors and constraints
  # Detection 
  lambda0 ~ dunif(0, 1) # tilda takes left thing and puts it in x spot
  sigma ~ dunif(0, 10000) #will also need multiple sigmas and lambda0s for covs
  
  # Data Augmentation 
  omega ~ dbeta(5,9) #somewhat informed prior for superpopulation
  
  # Survival and Recruitment Rates
    # Dirichlet prior for entry probabilities
  beta[1:n.prim.occasions] ~ ddirch(alpha[1:n.prim.occasions])

  # Convert entry probs to conditional entry probs
  b[1] <- beta[1]
  for (tp in 2:n.prim.occasions){
    b[tp] <- beta[tp]/prod(1 - b[1:(tp-1)]) 
    #Survival per primary
    phi[(tp-1)] ~ dunif(0, 1)
  } #t
  
  # Convert rates to probabilities
  G[1:3,1:3,1:(n.prim.occasions-1)] <- JSguts$transprobs(phi[1:(n.prim.occasions-1)], b[1:n.prim.occasions])  ## pvdb Index issue with b. Added 1.
  
  probstate1[1:3] <- c(1-b[1], b[1], 0) #create variable for first primary alive prob
  # Likelihood
  for (i in 1:M){
    #Home range center
    S[i, 1] ~ dunif(1, upperlimitx)    # indices on the matrix from 1 to dim + 1
    S[i, 2] ~ dunif(1, upperlimity)    # y coord of activity centers
    ones[i] ~ JSguts$dhab(S[i,1:2]) #the ones trick

    #Data Augmentation
    real[i] ~ dbern(omega) 
    initstate[i] ~ dcat(probstate1[1:3])
    Jprobs[i,1:J] <- JSguts$probdetected(lambda0, sigma, S[i,1:2])
    id[i] ~ dcapt_forward$dcapt_marg(real = real[i], Jprobs = Jprobs[i,1:J],
                                            G = G[1:3,1:3,1:(n.prim.occasions-1)], init = initstate[i])
  }
  #Derived parameters
  ## THIS WILL NEED A FORWARD BACKWARD SOLVE!!!
  ## ******************************************
  #for(t in 1:n.prim.occasions) {
   # N[t] <- sum(z[1:M,t,2]*real[1:M]) # no. alive for each year
    #change to sum of z that equal 2
  #}
  #Nsuper <- sum(real[1:M])
})


#data prep
m <- readRDS("m_sal5.Rds") #put it in repo
n.fake.inds <- 1000
ch <- m$data()$capthist()
traps <- m$data()$traps()
distmat <- m$data()$distances()
mesh <- m$data()$mesh()
primary <- m$data()$primary()
n.prim.occasions <- length(unique(primary))
dt <- diff(m$data()$time())
n.sec.occasions <- sum(primary %in% 1:n.prim.occasions)
trapusage <- usage(traps)

#create capture history that has a row for no detections per occasion
ch0 <- array(data = NA, dim = c(dim(ch)[1:2], dim(ch)[3]+1))
ch0[1:dim(ch)[1], 1:dim(ch)[2], 1:dim(ch)[3]] <- ch
ch0[1:dim(ch)[1], 1:dim(ch)[2], (dim(ch)[3]+1)] <- ifelse(apply(ch, c(1,2), sum)>0, 0, 1)
dimnames(ch0) <- list(dimnames(ch)[[1]], dimnames(ch)[[2]], c(dimnames(ch)[[3]], "no_detection"))

#data augmentation
nodets <- matrix(0, nrow = n.sec.occasions, ncol = dim(ch0)[3])
nodets[,dim(ch0)[3]] <- 1
aug <- aperm(replicate(n.fake.inds, nodets), c(3,1,2))
augch <- abind::abind(ch0, aug, along = 1)
dimnames(augch)[2:3] <- dimnames(ch0)[2:3]


#convert mesh to raster matrix of 1 and 
habmat <- rasterFromXYZ(cbind(mesh, rep(1, nrow(mesh)))) #change from NA to 0 and handle log0
whichmesh <- rasterFromXYZ(cbind(mesh, 1:nrow(mesh)))
real <- c(rep(1, dim(ch0)[1]), rep(NA, n.fake.inds))
M <- length(real)
habmat <- as.matrix(habmat)
whichmesh <- as.matrix(whichmesh)

#reasonable initial states
firstdets <- apply(apply(ch, c(1,2), sum), 1, FUN= function(x){min(which(x > 0))})
guessinitstates <- ifelse(firstdets %in% 1:9, 2,1) 
#probability of recruiting in the first 3 primaries is about equal to
#probability of surviving them

#choose intitial S for mesh near traps
meshrowsneartraps <- which(mesh[, "x"] %in% (traps[,"x"]+400) & #offsets
                             mesh[, "y"] %in% (traps[,"y"]-1600)) 
neartrap <- array(whichmesh %in% meshrowsneartraps, dim=dim(whichmesh))
indsneartrap <- which(neartrap, arr.ind = T, useNames = T)
Sguess <- indsneartrap[sample(1:nrow(indsneartrap), M, replace = T),][,c(2,1)] #columns x and rows y

datay <-  apply(augch, c(1,2), FUN = function(x){which(x > 0)})
                
data <- list(id = 1:M,  # animal ID for indexing data.
             real = real,
             ones = rep(1, M))
#only pass in things the model code needs, not setup code

constants <- list(n.prim.occasions = n.prim.occasions,
                  alpha = rep(1, n.prim.occasions), 
                  J = nrow(traps),
                  M = M,
                  upperlimitx = ncol(habmat)+1,
                  upperlimity = nrow(habmat)+1
)

inits <- list(
  lambda0 = m$get_par("lambda0", m = 1, j = 1, s= 1, k = 1),
  sigma = m$get_par("sigma", m = 1, j = 1, s= 1, k = 1),
  phi = m$get_par("phi", m = 1, j = 1, s = 1), 
  beta = m$get_par("beta", m = 1, j = 1, s = 1), #rdirch(n.prim.occasions, 1), #rdirch for more than one n doesn't work
  initstate = c(guessinitstates, sample(1:3,n.fake.inds, replace = T)),
  S = Sguess
)

JSguts <- JSguts_nf(habmat, whichmesh, distmat, dt)
dcapt_forward <- dcapt_forward_internal(datay, trapusage, primary)


rm(m) ## Remove m as it is big!
rm(ch0)
rm(augch)
rm(ch)
rm(aug)
gc()

startdefineT <- Sys.time()
model <- nimbleModel(code = JS_SCR, 
                      constants = constants, 
                      data = data,
                      init = inits#,
                      # calculate = FALSE#, #avoids calculation step, can do model$calculate later if it was too slow
                      # check = FALSE #won't check to make sure everything's right
)
totaldefineT <- Sys.time() - startdefineT #now takes ~40 seconds

model$simulate("S")
model$simulate("initstate")
model$calculate("Jprobs")
model$calculate("id[1]")  ## check timing in C++
model$Jprobs
model$simulate("initstate")
J <- nrow(traps)
dcapt_forward$dcapt_marg(x = 1, real = 1, Jprobs = model$Jprobs[1,],
                                            G = model$G[1:3,1:3,1:10], init = model$initstate[1])

# test <- compileNimble(dcapt_forward)
# test$dcapt_marg(1, 1, Jprobs = rep(0.25, nrow(traps)), G = NimbleJSmodel$G, c(0,1,0))
# dcapt_forward$dcapt_marg(x=1, real = 1, Jprobs = NimbleJSmodel$Jprobs[1,]+0.05, G = NimbleJSmodel$G, init = 2)

#if you want to try to figure out errors to do with funcitons, try compiling functions by themselves
#compiling the full model will do all of it
#example: compileNimble(JSguts) to check that it compiles, then I can run everythign within in c++ so its faster to check

cmodel <- compileNimble(model)

cmodel$calculate("G")  ## check timing in C++
cmodel$calculate("Jprobs")  ## check timing in C++
cmodel$calculate("id[1]")  ## check timing in C++
cmodel$calculate()  ## RUNS EVERYTHING - posterior log likelihood.

#conf <- configureMCMC #just defines mcmc, samplers for which
# if $setmontitors what to track

# model <- nimbleModel(modelCode, data = data, inits = inits, constants = constants) #, buildDerivs = TRUE) Dont' build derivs.
# conf <- configureMCMC(model)
conf <- configureMCMC(cmodel)
# conf$setMonitors(c("logalpha", "logbeta", "p0", "q", "sigma"))
mcmc <- buildMCMC(conf)
# cmodel <- compileNimble(model)
cmcmc <- compileNimble(mcmc, project = cmodel)
# # mcmc.out <- runMCMC(cmcmc, niter=50000, nburnin=5000, nchains=3, samplesAsCodaMCMC = TRUE)
cmcmc$run(100) #naybe 100 just to see if its slow, can access what everything is in the model
# mvSamples <- cmcmc$mvSamples
# samples <- as.matrix(mvSamples)
# out <- coda::mcmc(samples[-(1:5000),])	# Burn in
# mcmc.sum <- do.call("cbind", summary(out))

# conf$removeSamplers('X')
# for(i in 1:M){
#   conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
#                   type = 'RW_block', silent = TRUE))
# }
# 
# conf$removeSamplers('sigma')
# conf$addSampler(target = 'sigma', 
#                 type = 'slice', silent = TRUE, control = list(adaptive = TRUE, scaleWidth = 0.5))	
# 
# i just want point in space, random walk, i want to block sample acs per animal
#daniel turk, devalpine, efficient block samplers for MCMC large scale bayesian estimation of scr

#has context menu
#data is private$data_$covs(m = 1)
X <- m$design_mats() #setup from gams formulas
pars <- m$par()
for(i in 1:length(pars)){
  pars[[i]] 
}
#jagam will produce smooth code (mgcv function) and will need penalty too (penalty matrix will have prior)
#Marra 2011 paper

#still need to incorporate xy smooth, time smooth beta and phi, and spatial strata
# par_sal <- list(lambda0 ~ trapstratum + trapopen + periodf, #design matrix has row for each trap on each occasion
#                 sigma ~ trapstratum + trapopen + periodf, #same as lambda0
#                 beta ~ s(realtime, k = 6), 
#                 phi ~ s(realtime, k = 3), 
#                 D ~ s(x, y, k = 20) + s(salinity, k = 5))
