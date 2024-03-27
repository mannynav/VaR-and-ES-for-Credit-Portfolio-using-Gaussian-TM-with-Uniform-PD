
library(actuar)
library(mixtools)
library(fitdistrplus)
library(matrixdist)

#Expectation-Maximization for matrix distribution.
x_ini <- ph(structure = "general", dimension = 5)
x_fit <- fit(x_ini, AdjustedSample)
x_fit
p <- seq(0.01, 0.99, 0.01)
sq <- stats::quantile(AdjustedSample, p)
tq <- quan(x_fit, p)
plot(sq, tq)
abline(0, 1, col = "blue")

#####Code snippets used from 202223-Math337 - Financial and Actuarial Mathematics in R ####

#Density function for mixture distribution.
dmgp <- function(x, shapeg1, rateg1, shapeb1, shapeb2,shapeg2, rateg2) {
  dgamma(x, shapeg1, rateg1)*0.245 + dbeta(x, shapeb1, shapeb2)*0.51 + dgamma(x,shapeg2,rateg2)*0.245
}

#Distribution function for mixture distribution.
pmgp <- function(q, shapeg1, rateg1, shapeb1, shapeb2, shapeg2,rateg2) {
  pgamma(q, shapeg1, rateg1)*0.245 + pbeta(q, shapeb1, shapeb2)*0.51 + pgamma(q,shapeg2, rateg2)*0.245
}

#Quantile function for mixture distribution.
qmgp <- function(p, shapeg1, rateg1, shapeb1, shapeb2,shapeg2, rateg2) {
  
  L2 <- function(q, p) {
    (p - pmgp(q, shapeg1, rateg1, shapeb1, shapeb2, shapeg2, rateg2))^2}
  
  sapply(p, function(p) optimize(L2, c(0, 0.055), p = p)$minimum)
}

#Fit the adjusted values for X1,X2,...,Xn to the mixture distribution using MLE.
fit_mgp <- fitdist(AdjustedSample,distr = "mgp",start = 
                     list(shapeg1 = 1, rateg1 = 1, shapeb1 = 1, shapeb2 = 1,shapeg2 = 1, rateg2 = 2),lower = 0)


summary(fit_mgp)
plot(fit_mgp,breaks = 50)
