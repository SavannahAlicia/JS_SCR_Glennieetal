library(nimble)
library(secr)
library(raster)

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
  omega ~ dbeta(5,9) #somewhat informed prior
  #survival
  for (tp in 1:(n.prim.occasions-1)){
    phi[tp] ~ dunif(0, 1)
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
  
  # State membership transition probability 
  for (tp in 1:(n.prim.occasions - 1)){
    #Q is rates
    Q[1:3,1:3,tp] <- c( #note recruitment index starts at 1, but transitions happen between occasions 
      log(1 - b[(tp+1)])/dt[tp],                   0,          0,
      -log(1 - b[(tp+1)])/dt[tp],       log(phi[tp]),          0,
      0,                               -log(phi[tp]),          0
    )
    #have to fill matrix by column, not row
    G[1, 1:3, tp] <- c(1-b[(tp+1)], b[(tp+1)], 0) #from unborn
    G[2:3, 1, tp] <- c(0,0,0) #cannot transition back to unborn
    #Qd[1:2,1:2,tp] <- (Q[2:3,2:3,tp] * dt[tp])
    G[2:3,2:3,tp] <- Rexpm((Q[2:3,2:3,tp] * dt[tp])) #matrix exponential for competing states
  }
  probstate1[1:3] <- c(1-b[1], b[1], 0) #create variable for first primary alive prob
  # Likelihood
  for (i in 1:M){
    #Home range center
    S[i, 1] ~ dunif(1, upperlimitx)    # x-coord of activity centers
    S[i, 2] ~ dunif(1, upperlimity)    # y coord of activity centers
    hab[i] <- habMat[trunc(S[i, 1]), trunc(S[i, 2])] #sample x and y space, then check if AC is within mesh
    ones[i] ~ dbern(hab[i]) #the ones trick
    Smesh[i] <- whichmesh[trunc(S[i, 1]), trunc(S[i, 2])] #index of mesh
    #Data Augmentation
    real[i] ~ dbern(omega) 
    #State process
    z[i, 1, 1:3] ~ dmulti(prob = probstate1[1:3], size = 1)  #first primary, i is alive if it recruits into primary 1
    for (t in 2:n.prim.occasions){
      probstate[i,(t-1), 1:3] <- z[i, (t-1), 1:3] %*% G[1:3,1:3,(t-1)]
      z[i, t, 1:3] ~ dmulti(prob = probstate[i,(t-1),1:3], size = 1)
    }#t (primary)
    #Observation process (multi-catch trap)
    for (s in 1:n.sec.occasions){
      for(j in 1:J){
        hus[j,s,i] <- lambda0 * exp(-( distmat[j, Smesh[i]] *  distmat[j, Smesh[i]])/(2*sigma2)) * usage.traps[j, s]
        exphus[j,s,i] <- exp(-hus[j,s,i])
        Jprobs[i,s,j] <- hus[j,s,i]/sum(hus[1:J,s,i])*(1 - prod(exphus[1:J,s,i]))
      } #j
      Jprobs[i,s,J+1] <- prod(exphus[1:J,s,i]) #J+1 indexed outcome is no detection
      mu[i,s,1:(J+1)] <- Jprobs[i,s,1:(J+1)] * z[i, primary[s], 2]  * real[i] #is this right for data augmentation?
      #I think the trap of detection is multinomial, but then the alive, present, and real are bernoulli... do I need a custom distribution?
       y[i,s,1:(J+1)] ~ dmulti(prob = mu[i,s,1:(J+1)], size = 1) #multinomial, either detected at one trap j or not detected
    } #s (secondarys)
  } #i (individual)
  #Derived parameters
  for(t in 1:n.prim.occasions) {
    N[t] <- sum(z[1:M,t,2]*real[1:M]) # no. alive for each year
  }
  Nsuper <- sum(real[1:M])
})


#data prep
m <- readRDS("~/Documents/UniStAndrews/BarBay_OpenSCR/results/m_sal5.Rds")
n.prim.occasions <- 3
n.fake.inds <- 1000
ch <- m$data()$capthist()
traps <- m$data()$traps()
distmat <- m$data()$distances()
mesh <- m$data()$mesh()
primary <- m$data()$primary()
usage.traps <- usage(traps)

#create capture history that has a row for no detections per occasion
ch0 <- array(data = NA, dim = c(dim(ch)[1:2], dim(ch)[3]+1))
ch0[1:dim(ch)[1], 1:dim(ch)[2], 1:dim(ch)[3]] <- ch
ch0[1:dim(ch)[1], 1:dim(ch)[2], (dim(ch)[3]+1)] <- ifelse(apply(ch, c(1,2), sum)>0, 0, 1)
dimnames(ch0) <- list(dimnames(ch)[[1]], dimnames(ch)[[2]], c(dimnames(ch)[[3]], "no_detection"))

#inits from fitted model
beta <- m$get_par("beta", m = 1)
phi <- m$get_par("phi", m = 1)
dt <- diff(m$data()$time())
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
M <- length(real)


data <- list(y = augch13, 
             real = real,
             habMat = as.matrix(habMat),
             whichmesh = as.matrix(whichmesh),
             usage.traps = usage.traps13,
             distmat = distmat,
             dt = dt13,
             upperlimitx = ncol(habMat)+1,
             upperlimity = nrow(habMat)+1,
             ones = rep(1, M))

constants <- list(n.prim.occasions = n.prim.occasions,
                  n.sec.occasions = n.sec.occasions,
                  J = nrow(traps),
                  M = M,
                  primary = primary13
)

inits <- list(
  phi = m$get_par("phi", m = 1, j = 1, s = 1)[1:2], 
  lambda0 = m$get_par("lambda0", m = 1, j = 1, s= 1, k = 1),
  sigma = m$get_par("sigma", m = 1, j = 1, s= 1, k = 1),
  beta = m$get_par("beta", m = 1, j = 1, s = 1)[1:3] #rdirch(n.prim.occasions, 1), #rdirch for more than one n doesn't work
)

Rmodel <- nimbleModel(code = JS_SCR, 
                      constants = constants, 
                      data = data,
                      init = inits
)

#still need to incorporate xy smooth, time smooth beta and phi, and spatial strata
#list(lambda0 ~ trapstratum + trapopen + periodf, 
#     sigma ~ trapstratum + trapopen + periodf, 
#     beta ~ s(realtime, k = k_vals[i,1]), 
#     phi ~ s(realtime, k = k_vals[i,2]), 
#     D ~ s(x, y, k = D_k) + s(salinity, k = sal_k))
