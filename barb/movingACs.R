## Open Population analysis of Barataria Bay 
## Author: Richard Glennie
library(openpopscr)
# set number of cores 
RcppParallel::setThreadOptions(numThreads = 2)

# Data --------------------------------------------------------------------

barb <- readRDS("data/barb.Rds")
start <- list(lambda0 = exp(-5.3), sigma = exp(9.0), 
              beta = 0.5, phi = plogis(1.6), 
              D = exp(1.4),
              sd = 4000)
m_D <- readRDS(paste0("results/m_D_s20.Rds"))
mkval <- readRDS(paste0("results/m_popk63.Rds")) 
par_sal <- list(lambda0 ~ trapstratum + trapopen + periodf, 
                  sigma ~ trapstratum + trapopen + periodf, 
                  beta ~ s(realtime, k = 6), 
                  sd ~ 1,
                  phi ~ s(realtime, k = 3), 
                  D ~ s(x, y, k = 20) + s(salinity, k = 5))
m_sal <- JsTransientModel$new(par_sal, barb, start) 
tmp <- m_sal$par()
tmp$lambda0 <- m_D$par()$lambda0
tmp$sigma <- m_D$par()$sigma 
tmp$phi <- mkval$mle()$phi
tmp$beta <- mkval$mle()$beta
m_sal$set_par(tmp)
m_sal$fit(list(iterlim = 1000, stepmax = 10))
saveRDS(m_sal, file = "results/msal_transient.Rds")  
  
  