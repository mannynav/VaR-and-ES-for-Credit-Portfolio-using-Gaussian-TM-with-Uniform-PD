

#Packages.
library(fitdistrplus)
library(GCPM)

#Size of credit portfolio under study.
PortSize <- length(as.vector(na.omit(as.numeric(lgd$V3))))

#Set correlation.
beta <- 0.20

#Vector to store n-1 summed variables X_1 + X_2 +...+X_{n-1}.
vectS<-c(0)

#Vector of homogenous correlations for each loan.
corrVect<- rep(beta, PortSize)

#Number of simulations.
RUPD <- 1500

#Deterministic Exposure at default of 1.
EAD <- rep(1,PortSize)

#Loss given default from data set.
LGD <- suppressWarnings(as.vector(na.omit(as.numeric(lgd$V3))))

#Vector for uniform default probability among obligors. Can implement stochastic PD as well.
PD<- rep(0.05, PortSize)

#for loop to generate vector of S_n-1 = X_1 + X_2 + ... + X_n-1.
for (i in 1:RUPD) {
  
  #One factor model with F ~ N(0,1).
  F <- rnorm(1)
  
  #Idiosyncratic variables.
  epsilon <- rnorm(PortSize)
  
  #The critical variable X. If X > threshold, no default. If X <= threshold, obligor defaulted.
  X <- sqrt(corrVect)*F + sqrt(1-corrVect)*epsilon

  #For default only model, the threshold is di = qnorm(PDi).
  threshold <- qnorm(PD)
  
  #Apply threshold to critical variable X.
  defaultedLoans <- X < threshold
  
  #Get obligors that have defaulted.
  defaultedLoans[defaultedLoans == 'TRUE'] <- 1
  
  #Vector of obligors that have defaulted.
  lossVariableUPD  <- EAD * LGD * defaultedLoans
  
  #Vector to store sum n-1 of for each simulation i.
  vectS[i] <- sum(lossVariableUPD[1:length(lossVariableUPD)-1])
}

#Sample of portfolio used to fit a distribution to X1,X2,...,Xn.
defaultedLoans[defaultedLoans == 'TRUE'] <- 1
sample <- EAD*LGD*defaultedLoans
hist(sample,breaks = 50,freq = FALSE)

#Adjust sample.
AdjustedSample <- sample[sample>0.001]

#Fit adjusted sample to five distributions.
fit_betaSPD <- fitdist(AdjustedSample,"beta", start = list(shape1 = 0.01, shape2 = 0.01))
fit_paretoSPD <- fitdist(AdjustedSample, distr = "pareto", start = list(shape = 1, scale = 1))
fit_lognormSPD <- fitdist(AdjustedSample, distr = "lnorm")
fit_weibullSPD <- fitdist(AdjustedSample,"weibull",lower = c(0,0))
fit_gammaSPD<-fitdist(AdjustedSample, "gamma")

#Plots of fits.
plot(fit_betaSPD)
plot(fit_paretoSPD)
plot(fit_lognormSPD)
plot(fit_weibullSPD)
plot(fit_gammaSPD)

#Summaries of fits.
summary(fit_betaSPD)
summary(fit_paretoSPD)
summary(fit_lognormSPD)
summary(fit_weibullSPD)
summary(fit_gammaSPD)

#Set total losses.
hist(vectS,breaks = 50, xlab = "Losses",xlim = c(0,300),main = c("Histogram of Loss Distribution"))
summary(vectS)

#This is the conditional MC estimate for each S_(n-1) in vectS
Fncond <- function(x,vectS) {
  s<-c(0) 
  for (i in 1:RUPD) {
    #This distribution function is for the mixture fit. See MixtureFitExample4.9
    s[i] <- pmgp(x-vectS[i], fit_mgp$estimate[1], fit_mgp$estimate[2],fit_mgp$estimate[3], fit_mgp$estimate[4],fit_mgp$estimate[5],fit_mgp$estimate[6])
  }
  return(mean(s))
}

#Quantile.
alpha <- 0.95

#Solve Fncond(1/(1-q)-1,ScaledTotalLossesUPD) = alphaUPD for q.
qu <- uniroot(f = function(q) alpha - Fncond(1/ (1 - q) - 1,vectS), interval = c(0,1))$root

VaR <- 1 /(1 - qu) - 1

#Compute inverse of estimated Fncond, using the ScaledTotalLossesUPD.
FnInsCond <- function(x) {
  sINV<-c(0) 
  for (i in 1:RUPD) {
    sINV[i] <- pmgp(x-vectS[i],  fit_mgp$estimate[1], fit_mgp$estimate[2],fit_mgp$estimate[3], fit_mgp$estimate[4],fit_mgp$estimate[5],fit_mgp$estimate[6])
  }
  return(1-mean(sINV))
}

#Vectorize the inverse distribution function.
FnInvsCondVect <- Vectorize(FnInsCond, "x")

#Integrate the inverse distribution function from the VaR estimate to infinity to get the estimate for m.
m <- integral(FnInvsCondVect,VaR,Inf)

#Compute expected shortfall.
ES - VaR + (m)/(1-alpha)

