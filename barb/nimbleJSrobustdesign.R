library(nimble)
library(secr)
library(raster)
library(expm)

#need functions sum, prod, which, expm::expm, usage,
expmexpm <- function(x){
  expm::expm(x)
}
Rexpm <- nimbleRcall(function(x = double(2)){}, Rfun = 'expmexpm',
                     returnType = double(2))
Rt <- nimbleRcall(function(x = double(2)){}, Rfun = 't',
                  returnType = double(2))

JSguts_nf <- nimbleFunction(
  setup = function(habmat, whichmesh, distmat, trapusage, dt){
    habmat <- as.matrix(habmat) #just in case stored as df or raster
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
    probdetect = function(lambda = double(), 
                           sigma = double(), 
                           S = double(1)){
      smesh <- whichmesh[trunc(S[1]), trunc(S[2])] #index of mesh
      hus0 <- lambda*exp(-dist2[1:J, smesh]/(2*sigma^2))
      jprobs <- matrix(0, nrow = n.sec.occasions, ncol = J+1)
      for( s in 1:n.sec.occasions ){
        hus <- hus0*trapusage[,s]
        jprobs[s, J + 1] <- exp(-sum(hus))
        jprobs[s, 1:J] <- hus/sum(hus) * (1 - jprobs[s, J + 1])
      }
      returnType(double(2))
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
        Q[1:3,1:3] <- c( #note recruitment index starts at 1, but transitions happen between occasions 
          #log(1 - b[(tp+1)])/dt[tp],                   0,          0,
          NA,                                           0,          0,
          #-log(1 - b[(tp+1)])/dt[tp],       log(phi[tp]),          0, #causes errors if b is 1
          NA,                                log(phi[tp]),          0,
          0,                               -log(phi[tp]),          0
        )
        #have to fill matrix by column, not row
        G[1, 1:3] <- c(1-b[(tp+1)], b[(tp+1)], 0) #from unborn
        G[2:3, 1] <- c(0,0) #cannot transition back to unborn
        G[2:3,2:3] <- Rexpm((Q[2:3,2:3] * dt[tp])) #matrix exponential for competing states
        returnType(double(2))
        return(G)
    }#,
    # dmulti_trap = function(x = double(0), #from 1 to J+1 where J+1 is no detection
    #                        prob = double(1),
    #                        z = double(0),
    #                        real = double(0),
    #                        log = logical(0, default = FALSE)){
    #   ans <- 0
    #   if (z == 1 & real == 1){ #if alive and real (note could check without this logic in case logic breaks things)
    #     ans <- log(prob[x])
    #   }
    #   returnType(double(0))
    #   if (log){
    #     return(ans)
    #   } else {
    #     return(exp(ans))
    #   } 
    # }
  )
)





#spatial
JS_SCR <- nimbleCode({
  
  # Priors and constraints
  # Detection 
  lambda0 ~ dunif(0, 1) # tilda takes left thing and puts it in x spot
  #logit(lambda0) = sum(covs)
  #lambda0Intercept ~ dnorm()(think about priors)
  #lambda0trapsstratumSoutheast
  #lambda0trapstratumCentral
  #lambda0trapstratumWest
  #lambda0trapopensheltered
  #lambda0periodf2
  #lambda0periodf3
  #lambda0periodf4
  #lambda0periodf5
  #lambda0periodf6
  #lambda0periodf7
  #lambda0periodf8
  #lambda0periodf9
  #lambda0periodf10
  #lambda0periodf11
  #also how to make sure that sum of possible lambdas is between 0 and 1?
  sigma ~ dunif(0, 10000) #will also need multiple sigmas
  #sigmaIntercept
  #sigmatrapsstratumSoutheast
  #sigmatrapstratumCentral
  #sigmatrapstratumWest
  #sigmatrapopensheltered
  #sigmaperiodf2
  #sigmaperiodf3
  #sigmaperiodf4
  #sigmaperiodf5
  #sigmaperiodf6
  #sigmaperiodf7
  #sigmaperiodf8
  #sigmaperiodf9
  #sigmaperiodf10
  #sigmaperiodf11
  
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
  G[1:3,1:3,1:(n.prim.occasions-1)] <- JSguts$transprobs(phi[1:(n.prim.occasions-1)], b[1:(n.prim.occasions-1)])
  
  probstate1[1:3] <- c(1-b[1], b[1], 0) #create variable for first primary alive prob
  # Likelihood
  for (i in 1:M){
    #Home range center
    S[i, 1] ~ dunif(1, upperlimitx)    # indices on the matrix from 1 to dim + 1
    S[i, 2] ~ dunif(1, upperlimity)    # y coord of activity centers
    ones[i] ~ JSguts$dhab(S[i,1:2]) #the ones trick

    #Data Augmentation
    real[i] ~ dbern(omega) 
    #State process (z is now index of state (1,2,3))
    z[i, 1] ~ dcat(prob = probstate1[1:3])  #first primary, i is alive if it recruits into primary 1
    for (t in 2:n.prim.occasions){ #just the row of G for the state z
      z[i, t] ~ dcat(prob = G[z[i,t-1],1:3,(t-1)])
    }#t (primary)
    #Observation process (multi-catch trap)
    Jprobs[i,1:n.sec.occasions,1:(J+1)] <- JSguts$probdetect(lambda0, sigma, S[i,1:2])
    for (s in 1:n.sec.occasions){
      mu[i,s,1:(J+1)] <- Jprobs[i,s,1:(J+1)] * (z[i, primary[s]] == 2)  * real[i] #is this right for data augmentation?
      #I think the trap of detection is multinomial
       y[i,s] ~ dcat(prob = mu[i,s,1:(J+1)]) #y is which trap detected or J+1 for none
    } #s (secondarys)
  } #i (individual)
  #Derived parameters
  #for(t in 1:n.prim.occasions) {
   # N[t] <- sum(z[1:M,t,2]*real[1:M]) # no. alive for each year
    #change to sum of z that equal 2
  #}
  #Nsuper <- sum(real[1:M])
})


#data prep
m <- readRDS("~/Documents/UniStAndrews/BarBay_OpenSCR/results/m_sal5.Rds")
n.fake.inds <- 1000
ch <- m$data()$capthist()
traps <- m$data()$traps()
distmat <- m$data()$distances()
mesh <- m$data()$mesh()
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
habmat <- rasterFromXYZ(cbind(mesh, rep(1, nrow(mesh))))
whichmesh <- rasterFromXYZ(cbind(mesh, 1:nrow(mesh)))
real <- c(rep(1, dim(ch0)[1]), rep(NA, n.fake.inds))
M <- length(real)

trapusage <- usage(traps)
habmat <- as.matrix(habmat)
whichmesh <- as.matrix(whichmesh)

data <- list(y = apply(augch, c(1,2), FUN = function(x){which(x > 0)}), #index trap (or lack of dets)
             real = real,
             ones = rep(1, M),
             upperlimitx = ncol(habmat)+1,
             upperlimity = nrow(habmat)+1)
#only pass in things the model code needs, not setup code

constants <- list(n.prim.occasions = n.prim.occasions,
                  n.sec.occasions = n.sec.occasions,
                  alpha = rep(1, n.prim.occasions), #make this a constant instead
                  J = nrow(traps),
                  M = M,
                  primary = primary
)

inits <- list(
  lambda0 = m$get_par("lambda0", m = 1, j = 1, s= 1, k = 1),
  sigma = m$get_par("sigma", m = 1, j = 1, s= 1, k = 1),
  phi = m$get_par("phi", m = 1, j = 1, s = 1), 
  beta = m$get_par("beta", m = 1, j = 1, s = 1) #rdirch(n.prim.occasions, 1), #rdirch for more than one n doesn't work
)

JSguts <- JSguts_nf(habmat, whichmesh, distmat, trapusage, dt)

NimbleJSmodel <- nimbleModel(code = JS_SCR, 
                      constants = constants, 
                      data = data,
                      init = inits,
                      calculate = FALSE#, #avoids calculation step, can do model$calculate later if it was too slow
                      #check = false #won't check to make sure everything's right
)

NimbleJSmodel_calculated <- NimbleJSmodel$calculate()
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
