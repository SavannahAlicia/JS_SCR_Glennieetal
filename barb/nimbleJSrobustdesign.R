library(nimble)
library(secr)
library(raster)
library(expm)
library(ggplot2)

##Define Functions#-------------------------------------------------------------

## Nimble function that calls R matrix eponential from expm package.
## Will do this call in compiled code.
Rexpm <- nimbleRcall(function(x = double(2)){}, Rfun = 'expm',
                     returnType = double(2))

## Main methods function in Nimble that contains all the data to reduce the complexity of the model structure in NIMBLE
## Key is that we are going to put some cleverish caching in here so that we don't have to create too many jprobs in code.
## Setup includes a matrix of if habitat or not. 
## distmat is distances between each centroid on the mask and all detectors (Distance can be euclidean or non-euclidean).
## dt is the time between each primary occassion.
## dLimit is the maximum distance that an animal could travel from its activity centre. After that pdetect = 0.
JSguts_nf <- nimbleFunction(
  setup = function(habmat, whichmesh, distmat, dt, dLimit = 100){
    habmat <- as.matrix(habmat) #just in case stored as df or raster
    whichmesh <- as.matrix(whichmesh)
    J <- nrow(distmat)
    dist2 <- distmat*distmat
    n.sec.occasions <- ncol(trapusage)
    n.prim.occasions <- length(dt) + 1
    d2Limit <- dLimit^2
  },
  # run = function(){},
  methods = list(
    dhab = function(x = double(), #ones data
                     S = double(1), #array of length 2 x and y
                     log = logical(0, default = FALSE)){
      hab <- habmat[trunc(S[2]), trunc(S[1])] #choose from mesh COLUMNS ARE X AND ROWS ARE Y
      if(log) return(log(hab))
      else return(hab)
      returnType(double(0))
    },
    ## Method to calculate the probability of being detected based on a halfnormal detection function at each detector given activity centre.
    ## returns vector of length traps.
    probdetected = function(lambda = double(), 
                           sigma = double(), 
                           S = double(1)){
      smesh <- whichmesh[trunc(S[2]), trunc(S[1])]  #index of mesh COLUMNS ARE X AND ROWS ARE Y!!!!! *PVDB
      jprobs <- numeric(value = 0, length = J)
      if(is.na(smesh)) return(jprobs)
      keep <- which(dist2[1:J, smesh] <= d2Limit)
      jprobs[keep] <- lambda*exp(-dist2[keep, smesh]/(2*sigma^2))
      returnType(double(1))
      return(jprobs)
    },
    ## Returns transition probabilty matrix.
    transprobs = function(phi = double(1), 
                          b = double(1)){ #note could vectorize to speed up
      G <- array(NA, dim = c(3,3,(n.prim.occasions - 1)))
      for(tp in 1:(n.prim.occasions - 1)){
        G[1:3, 1:3, tp] <- calcGt(phi[1:(n.prim.occasions - 1)], b[1:n.prim.occasions], tp)
      }
      returnType(double(3))
      return(G)
    },
    ## internal method for computing the transition matrix.
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

## distribution function for marginalizing over all detections of a single animal.
## holds all the data for memory issues in NIMBLE and compiling speed.
## capthist - capture history
## trapusage - trap active
## primary - label of primary occassion for each individual (secondary) capture.
dcapt_forward_internal <- nimbleFunction(
  setup = function(capthist, trapusage, primary){
    J <- nrow(trapusage)
    nprimary <- max(primary)
    
  }, # Hold data to make model smaller.
  run = function(){},
  methods = list(
    ## capture history marginal distribution using forward solve. Assumes we are sampling initial state. 
    ## We could remove this and marginalize the extra step in the future.
    dcapt_marg = function(x = double(), # index of animal 1:M
           real = double(),     # real or not real animal?
           Jprobs = double(1), # probability of detection for each secondary at each trap (single individual)
           G = double(3), # Transition prob matrix for each primary
           init = double(0),
           log = logical(0, default=TRUE)){
    ch <- capthist[x,]    #vector of which trap (or no trap) made first detection for occasion    
    ## Check if it's not real then prob = 1, unless it was observed, then shit.
    isobs0 <- all(ch == J+1) #all unobserved
    if(real == 0){ # if augmented individual
      ans <- log(1) # it must be unobserved
      if(!isobs0){ # if it was observed
        ans <- log(0) # that's not possible
      }
    } else { #otherwise its a real individual, in which case detection history depends on state
      ## Now to forward solve.
      ## State 1: unborn
      ## State 2: alive
      ## State 3: dead
      nstates <- dim(G)[1]
      pi <- nimNumeric(value = 0, length = nstates) # prob of state as we go forward
      pobs <- nimNumeric(value = 0, length = nstates) # Pdetect given state
      logL <- 0
      pi[init] <- 1 #init specifies number of state we begin in
      keep <- which(Jprobs > 0)  ## Save a little compute time maybe...

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
          } else { #if alive
            lp <- 0
            for( k in seq_along(idx) ) {
              occ <- idx[k]
              ## Need to compute detection probs:
              ## Part 1: Not detected but keep on log scale
              logpnotdet <- -sum(trapusage[keep,occ]*Jprobs[keep])
              ## If not detected then I'm done. Otherwise calc the prob of being det at just that trap.
              if(ch[occ] == J+1) {
                lp <- lp + logpnotdet
              }else{
                ## If can't be detected, but was then return -infinity.
                if(logpnotdet < 0) { #probability of no detection < 1
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
          if (sumpi != 0){ #will return NA for divide by 0, if all pi 0 just keep 0
          pi <- ((pi %*% G[,,tp])/sumpi)[1,] 
          }
        }
      }
      ans <- logL
    }
    
    returnType(double())
    if(log) return(ans)
    else return(exp(ans))
  }
  )
)

## Nimble code that is kept simple and clean due to data and constants stored in nimble function setup.
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
  G[1:3,1:3,1:(n.prim.occasions-1)] <- JSguts$transprobs(phi[1:(n.prim.occasions-1)], 
                                                         b[1:n.prim.occasions])  ## pvdb Index issue with b. Added 1.
  
  probstate1[1:3] <- c(1-b[1], b[1], 0) #create variable for first primary alive prob
  # Likelihood
  for (i in 1:M){
    #Home range center COLUMNS ARE X AND ROWS ARE Y
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

##Data Setup#-------------------------------------------------------------------
#still need to incorporate xy smooth, time smooth beta and phi, and spatial strata
# par_sal <- list(lambda0 ~ trapstratum + trapopen + periodf, #design matrix has row for each trap on each occasion
#                 sigma ~ trapstratum + trapopen + periodf, #same as lambda0
#                 beta ~ s(realtime, k = 6), 
#                 phi ~ s(realtime, k = 3), 
#                 D ~ s(x, y, k = 20) + s(salinity, k = 5))

m <- readRDS("results/m_sal5.Rds") #put it in repo
n.fake.inds <- 1000
ch <- m$data()$capthist()
#rescale to km
distmat <- m$data()$distances()
distmat <- distmat/1000
mesh <- m$data()$mesh()
yoff <- min(mesh$y)
xoff <- min(mesh$x)
mesh <- data.frame(x = (mesh$x-xoff)/1000, 
                     y = (mesh$y-yoff)/1000)
traps <- m$data()$traps()
trapusage <- usage(traps)
traps <- data.frame(x = (traps$x-xoff)/1000,
                      y = (traps$y-yoff)/1000)
primary <- m$data()$primary()
n.prim.occasions <- length(unique(primary))
dt <- diff(m$data()$time())
n.sec.occasions <- sum(primary %in% 1:n.prim.occasions)



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
habmat[is.na(habmat)] <- 0  ## PVDB: Upated to make sure NA is actually 0.
whichmesh <- rasterFromXYZ(cbind(mesh, 1:nrow(mesh)))
real <- c(rep(1, dim(ch0)[1]), rep(NA, n.fake.inds))
M <- length(real)
n.real.inds <- M - n.fake.inds
habmat <- as.matrix(habmat)
whichmesh <- as.matrix(whichmesh)

#reasonable initial states
firstdets <- apply(apply(ch, c(1,2), sum), 1, FUN= function(x){min(which(x > 0))})
guessinitstates <- ifelse(firstdets %in% 1:9, 2,1) 
#probability of recruiting in the first 3 primaries is about equal to
#probability of surviving them

#choose intitial S for mesh near traps
meshclosetrap <- function(#returns col x and row y of mesh closest to trap
    trapn){ #index of trap
  closestx <- which(abs(mesh[,"x"] - traps[trapn,"x"]) == 
                      min(abs(mesh[,"x"] - traps[trapn,"x"])))
  closesty <- which(abs(mesh[closestx,"y"] - traps[trapn,"y"]) == 
                      min(abs(mesh[closestx,"y"] - traps[trapn,"y"])))
  meshrowneartrap <- closestx[closesty]
  return(which(whichmesh == meshrowneartrap, arr.ind = T, useNames = T)[,c(2,1)] #columns x and rows y
         )
}

datay <-  apply(augch, c(1,2), FUN = function(x){which(x > 0)})
Sguess <- rbind(#real individuals
  do.call(rbind, lapply(as.list(apply(datay[1:n.real.inds,], 1, min)),  meshclosetrap)),
  do.call(rbind, lapply(as.list(sample(1:nrow(traps), n.fake.inds, replace = T)),  meshclosetrap))
)
                
data <- list(id = 1:M,  # animal ID for indexing data.
             real = real,
             ones = rep(1, M),
             initstate = ifelse(apply(as.array(datay[,which(primary ==1)]), 1, min) != 910, 2, NA)#put alive for those detected in first primary
             )
#only pass in things the model code needs, not setup code

constants <- list(n.prim.occasions = n.prim.occasions,
                  alpha = rep(1, n.prim.occasions), 
                  J = nrow(traps),
                  M = M,
                  upperlimitx = ncol(habmat)+1,
                  upperlimity = nrow(habmat)+1
)

realinit <- (!is.na(data$real))*1
realinit[realinit == 1] <- NA
inits <- list(
  lambda0 = m$get_par("lambda0", m = 1, j = 1, s= 1, k = 1),
  sigma = m$get_par("sigma", m = 1, j = 1, s= 1, k = 1)/1000,
  phi = m$get_par("phi", m = 1, j = 1, s = 1), 
  beta = m$get_par("beta", m = 1, j = 1, s = 1), #rdirch(n.prim.occasions, 1), #rdirch for more than one n doesn't work
  initstate = c(guessinitstates, sample(1:2, n.fake.inds, replace = T)), ## InitState can't be 3 Savannah!!
  S = Sguess,
  real = realinit,
  omega = .6
)

JSguts <- JSguts_nf(habmat, whichmesh, distmat, dt, dLimit = 30)
dcapt_forward <- dcapt_forward_internal(datay, trapusage, primary)


##Memory Cleanup#---------------------------------------------------------------
rm(m) ## Remove m as it is big!
rm(ch0)
rm(augch)
rm(ch)
rm(aug)
gc()

##Nimble Object Creation#-------------------------------------------------------
startdefineT <- Sys.time()
model <- nimbleModel(code = JS_SCR, 
                      constants = constants, 
                      data = data,
                      init = inits#,
                      # calculate = FALSE#, #avoids calculation step, can do model$calculate later if it was too slow
                      # check = FALSE #won't check to make sure everything's right
)
totaldefineT <- Sys.time() - startdefineT #now takes ~1 min

startcompileT <- Sys.time()
cmodel <- compileNimble(model)
totalcompileT <- difftime(startcompileT, Sys.time(), "secs") #30 sec

startcalcT <- Sys.time()
cmodel$calculate()
totalcalcT <- difftime(startcalcT, Sys.time(), "secs")

conf <- configureMCMC(cmodel)
conf$setMonitors(c("beta", "phi", "lambda0", "sigma", "S", "omega"))
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = cmodel)

start100mcmcT <- Sys.time()
cmcmc$run(10000, time = TRUE) #took 15 mins, which is 0.155 min per iter, 5 days for 50k
total100mcmcT <- difftime(Sys.time(), start100mcmcT, "secs")

# sum(cmcmc$getTimes())/100*50000/60/60 ## cmcmc$getTimes() returns the sample time for each sampler.

start10kmcmcT <- Sys.time()
mcmc.out <- runMCMC(cmcmc, niter= 30000, nburnin=10000, nchains=2, samplesAsCodaMCMC = TRUE)
total10kmcmcT <- difftime(Sys.time(), start10kmcmcT) #Time difference of 4.922812 days

saveRDS(mcmc.out, file = "mcmcout.Rds")
saveRDS(total10kmcmcT, file = "time.Rds")

mvSamples <- mcmc.out$chain1#cmcmc$mvSamples
samples <- as.matrix(mvSamples)
plotsampletrace <- function(par = "", numits, limits = c(0,max(samplematrix))){
  samplematrix <- samples[,which(grepl(par, colnames(samples)))] #regular expressions to pull parameters for plots
  if(!is.null(dim(samplematrix))){
    sampledf <- data.frame(name = tidyr::pivot_longer(as.data.frame(samplematrix), cols = 1:ncol(samplematrix))$name,
                           value = tidyr::pivot_longer(as.data.frame(samplematrix), cols = 1:ncol(samplematrix))$value,
                           iteration = rep(1:numits, each = ncol(samplematrix)))

  } else {
    sampledf <- data.frame(name = rep(par),
                           value = samplematrix,
                           iteration = 1:numits)
  }
  theplot <- ggplot(data = sampledf) +
    geom_line(mapping = aes(x = iteration, y = value, col = name), alpha = .5) +
    ylim(limits) +
    theme_classic()
  print(theplot)
  return(theplot)
}
# plotsampletrace("S\\[\\d{1,4}\\,\\s1.", 1800)
# S2sample <- samples[,(M+1):(M+M)]
# plotsampletrace("S.*, 1]", 1800)
# plotsampletrace("beta", 1800)
# plotsampletrace("phi", 1800, limits = c(.6,1))
# plotsampletrace("lambda0",1800)
# plotsampletrace("sigma",1800)
# plotsampletrace("omega", 1800, limits = c(.97,1))
# # out <- coda::mcmc(samples[-(1:5000),])	# Burn in
# # mcmc.sum <- do.call("cbind", summary(out))
# 
