needed.packages <- c("chron","lubridate","zoo","xts","data.table","highfrequency",
                     "timeSeries","quantmod","ggplot2","dplyr","urca","MASS",
                     "tidyr","diagonals","Matrix")
#lapply(needed.packages, install.packages, character.only = TRUE)
lapply(needed.packages, require, character.only = TRUE)

library(QuoteDynamics)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rm(list=ls())

options(digits.secs = 3)

load(file="IBM_quotes_100ms_2021-01-04.RData")
#load(file="IBM_quotes_100ms_2021-01-07.RData")

library(QuoteDynamics)
ls("package:QuoteDynamics")
# ?quoteDynOptim

q <- na.omit(filled_ts)

# Test, 2000 obs
q <- q[1:5000,]


q <- as.data.table(q)
q[, "tau":= difftime(index,lag(index))]
q <- q[ tau > 0 ]

q[, "ind.N" := (BID.N - lag(BID.N)) | (OFR.N - lag(OFR.N)) != 0]
q[, "ind.T" := (BID.T - lag(BID.T)) | (OFR.T - lag(OFR.T)) != 0]
q[, "ind.K" := (BID.K - lag(BID.K)) | (OFR.K - lag(OFR.K)) != 0]
q[, "ind.P" := (BID.P - lag(BID.P)) | (OFR.P - lag(OFR.P)) != 0]
q[, "ind.Z" := (BID.Z - lag(BID.Z)) | (OFR.Z - lag(OFR.Z)) != 0]

q <- q[ !(ind.N == FALSE & ind.T == FALSE & ind.K == FALSE &
            ind.P == FALSE & ind.Z == FALSE) ]

ts.plot(q$BID.N)
lines(q$OFR.N)
ts.plot(q$BID.T)
lines(q$OFR.T)

spread <- c(mean(q$OFR.N - q$BID.N), mean(q$OFR.T - q$BID.T), mean(q$OFR.K - q$BID.K),
            mean(q$OFR.P - q$BID.P), mean(q$OFR.Z - q$BID.Z))
y <- as.matrix(q[,c("BID.N","OFR.N","BID.T","OFR.T", "BID.K","OFR.K",
                    "BID.P","OFR.P", "BID.Z","OFR.Z")])
tau <- as.matrix(q$tau)
M <- 5
J <- nrow(y)
ind.mat <- as.matrix(q[,c("ind.N","ind.T","ind.K","ind.P","ind.Z")])

rowSums(ind.mat)

# initial values
alpha <- c(0.05,0.05, 0.05,0.05, 0.05,0.05,
           0.05,0.05, 0.05,0.05)
Omega1 <- matrix(c(0.0001,0,
                   0,0.0001), nrow=2, ncol=2, byrow=TRUE)
Omega2 <- matrix(c(0.0001,0,
                   0,0.0001), nrow=2, ncol=2, byrow=TRUE)
Omega3 <- matrix(c(0.0001,0,
                   0,0.0001), nrow=2, ncol=2, byrow=TRUE)
Omega4 <- matrix(c(0.0001,0,
                   0,0.0001), nrow=2, ncol=2, byrow=TRUE)
Omega5 <- matrix(c(0.0001,0,
                   0,0.0001), nrow=2, ncol=2, byrow=TRUE)

sigma <- 0.005
delta1 <- 0.1
delta2 <- -0.1

start <- c(spread, as.vector(alpha), t(chol(Omega1))[lower.tri(t(chol(Omega1)), diag=TRUE)],
           t(chol(Omega2))[lower.tri(t(chol(Omega2)), diag=TRUE)],
           t(chol(Omega3))[lower.tri(t(chol(Omega3)), diag=TRUE)],
           t(chol(Omega4))[lower.tri(t(chol(Omega4)), diag=TRUE)],
           t(chol(Omega5))[lower.tri(t(chol(Omega5)), diag=TRUE)], sigma, delta1, delta2)

# ptm <- proc.time()
# estResults <- QuoteDynamics::quoteDynOptim(start = start, X = y, tau = tau, xtol = 10e-3,
#                                            stop_val = 10e-10, algorithm = "LN_NELDERMEAD", hessian = TRUE,
#                                            step_size = 1e-04, verbose = TRUE)
# proc.time() - ptm
#
# par.est <- estResults$estimate
# lf <- estResults$min_val
# hess <- estResults$hessian
#
# vc <- (1/J)*solve(hess)
#
# cat('\nLog-likelihood function = ',lf)
# cat('\n')
# print( cbind("Param"=round(par.est,4),  "Std.Errors"=sqrt(diag(vc))))














# Selection matrix
S <- diag(2*M)
S1 <- S[1:2,]
S2 <- S[3:4,]

S1 %*% y[1,]
S2 %*% y[1,]




#b <- start

#--------------------------------------------------------------------------
# Multivariate Kalman filter
#--------------------------------------------------------------------------
lnlt <- function(b,y,tau,ind.mat) {

  # dimension of the observations: 2*M
  # dimension of the state variable: 2
  M <- ncol(y)/2
  # Parameters
  spread <- b[1:M]
  ct <- c(-spread[1]/2,spread[1]/2)
  for (i in 2:M) ct <- c(ct, -spread[i]/2,spread[i]/2)
  alpha <- b[(M+1):(3*M)]
  b <- b[-(1:(3*M))]
  C.mat <- matrix(0, nrow=2, ncol=2)
  C.mat[lower.tri(C.mat, diag=TRUE)] <- b[1:3]
  b <- b[-(1:3)]
  Omega <- C.mat %*% t(C.mat)
  GGt <- as.matrix(bdiag(Omega))
  for (i in 2:M) {
    C.mat <- matrix(0, nrow=2, ncol=2)
    C.mat[lower.tri(C.mat, diag=TRUE)] <- b[1:3]
    b <- b[-(1:3)]
    Omega <- C.mat %*% t(C.mat)
    GGt <- as.matrix(bdiag(GGt,Omega))
  }
  dt <- c(0,0)
  Tt <- matrix(c(1,0,
                 0,0), nrow=2, ncol=2)
  sigma <- b[1]
  delta1 <- b[2]
  delta2 <- b[3]

  # Allocate arrays
  t <- nrow(y)
  n <- ncol(y)
  lnl <- rep(0, t)

  # Recursions of the Kalman Filter
  # Initialisation
  st <- rep(0, 2)
  Pt <- diag(2)*1000

  Zt <- cbind(1,alpha*tau[1]^delta2*sigma)

  # Observation
  # conditional mean of yt
  mut <- ct + Zt %*% st
  # conditional variance of yt
  Ft <- Zt %*% Pt %*% t(Zt) + GGt
  # forecast error
  vt <- y[1,] - mut

  lnl[1] <- - 0.5*n*log(2*pi) - 0.5*log(det(Ft)) - 0.5*t(vt) %*% solve(Ft) %*% vt

  # Kalman gain
  Kt <- Pt %*% t(Zt) %*% solve(Ft)
  # latent state
  s0 <- st + Kt %*% vt
  # conditional variance
  P0 <- Pt - Pt %*% t(Zt) %*% t(Kt)

  # Main loop over observations
  for (i in 2:t) {

    # construct selection matrix
    S <- diag(2*M)
    index <- c()
    for (j in 1:M) {
      if (ind.mat[i,j] == FALSE) {index <- c(index, ((2*j-1):(2*j)))}
    }
    if(!is.null(index)){
      S <- S[-index, ]
    }
    # print(index)
    # print(i)
    # print(S)

    Zt <- S %*% cbind(1,alpha*tau[i]^delta2*sigma)
    H <- matrix(c(0,sigma*tau[i]^delta1,
                  0,1), nrow=2, ncol=2, byrow=TRUE)

    # updating
    st <- dt + Tt %*% s0
    Pt <- Tt %*% P0 %*% t(Tt) + H %*% t(H)

    # Observation
    # conditional mean of yt
    mut <- S %*% ct + Zt %*% st
    # conditional variance of yt
    Ft <- Zt %*% Pt %*% t(Zt) + S %*% GGt %*% t(S)
    # forecast error
    vt <- S %*% y[i,] - mut

    lnl[i] <- - 0.5*n*log(2*pi) - 0.5*log(det(Ft)) - 0.5*t(vt) %*% solve(Ft) %*% vt

    # updating
    # Kalman gain
    Kt <- Pt %*% t(Zt) %*% solve(Ft)
    # latent state
    s0 <- st + Kt %*% vt
    # conditional variance
    P0 <- Pt - Pt %*% t(Zt) %*% t(Kt)
  }
  return(lnl)

}

neglog <- function(b,y,tau,ind.mat) {
  lf <- -mean( lnlt(b,y,tau,ind.mat) )
  return(lf)
}

# initial values
alpha <- c(0.05,0.05, 0.05,0.05, 0.05,0.05,
           0.05,0.05, 0.05,0.05)
Omega1 <- matrix(c(0.001,0,
                   0,0.001), nrow=2, ncol=2, byrow=TRUE)
Omega2 <- matrix(c(0.001,0,
                   0,0.001), nrow=2, ncol=2, byrow=TRUE)
Omega3 <- matrix(c(0.001,0,
                   0,0.001), nrow=2, ncol=2, byrow=TRUE)
Omega4 <- matrix(c(0.001,0,
                   0,0.001), nrow=2, ncol=2, byrow=TRUE)
Omega5 <- matrix(c(0.001,0,
                   0,0.001), nrow=2, ncol=2, byrow=TRUE)

sigma <- 0.0005
delta1 <- 0.1
delta2 <- -0.2

start <- c(spread, as.vector(alpha), t(chol(Omega1))[lower.tri(t(chol(Omega1)), diag=TRUE)],
           t(chol(Omega2))[lower.tri(t(chol(Omega2)), diag=TRUE)],
           t(chol(Omega3))[lower.tri(t(chol(Omega3)), diag=TRUE)],
           t(chol(Omega4))[lower.tri(t(chol(Omega4)), diag=TRUE)],
           t(chol(Omega5))[lower.tri(t(chol(Omega5)), diag=TRUE)], sigma, delta1, delta2)
b <- start

# write.table(start, "start.txt", sep = ",", row.names = FALSE, col.names = FALSE)
# write.table(y, "y.txt", sep = ",", row.names = FALSE, col.names = FALSE)
# write.table(tau, "tau.txt", sep = ",", row.names = FALSE, col.names = FALSE)
# write.table(matrix(as.numeric(ind.mat), dim(ind.mat)[1], dim(ind.mat)[2]), "ind_mat.txt", sep = ",", row.names = FALSE, col.names = FALSE)

# con <- file("test_quoteDynOptim.log")
# sink(con, append=TRUE)
# sink(con, append=TRUE, type="message")
#
# Minimum <- QuoteDynamics::quoteDynOptim(start = start, X = y, tau = tau, ident_mat = ind.mat,
#                                         xtol = sqrt(.Machine$double.eps),
#                                         stop_val = sqrt(.Machine$double.eps),
#                                         step_size = sqrt(.Machine$double.eps),
#                                         max_eval = 1,
#                                         algorithm = "LN_COBYLA",
#                                         hessian = TRUE,
#                                         verbose = TRUE)
#
# # Restore output to console
# sink()
# sink(type="message")

PrintAlgorithms()
ptm_cpp_start <- proc.time()
Minimum <- QuoteDynamics::quoteDynOptim(start = start, X = y, tau = tau, ident_mat = ind.mat,
                                        xtol = sqrt(.Machine$double.eps),
                                        stop_val = sqrt(.Machine$double.eps),
                                        step_size = sqrt(.Machine$double.eps),
                                        max_eval = 100*length(start),
                                        algorithm = "LN_NELDERMEAD",
                                        hessian = TRUE,
                                        verbose = TRUE)
ptm_cpp_end <- proc.time()
ptm_cpp_end - ptm_cpp_start

ptm_R_start <- proc.time()
estResults <- optim(start, neglog, y=y, tau=tau, ind.mat=ind.mat, method="Nelder-Mead", hessian=FALSE,
                    control=list(reltol=sqrt(.Machine$double.eps),
                                 abstol=sqrt(.Machine$double.eps),
                                 maxit=100*length(start)))
ptm_R_end <- proc.time()
ptm_R_end - ptm_R_start



mean((start - estResults$par)^2)
mean((start - Minimum$estimate)^2)


# Test random pertubation
start <- start + rnorm(length(start))
ptm_R_start <- proc.time()
estResults <- optim(start, neglog, y=y, tau=tau, ind.mat=ind.mat, method="Nelder-Mead", hessian=FALSE,
                    control=list(reltol=sqrt(.Machine$double.eps), abstol=sqrt(.Machine$double.eps), maxit=100*length(start)))
ptm_R_end <- proc.time()
ptm_R_end - ptm_R_start

PrintAlgorithms()
ptm_cpp_start <- proc.time()
Minimum <- QuoteDynamics::quoteDynOptim(start = start, X = y, tau = tau, ident_mat = ind.mat,
                                        xtol = sqrt(.Machine$double.eps),
                                        stop_val = sqrt(.Machine$double.eps),
                                        step_size = sqrt(.Machine$double.eps),
                                        max_eval = 100 * length(start),
                                        algorithm = "LN_NELDERMEAD",
                                        hessian = FALSE,
                                        verbose = TRUE)
ptm_cpp_end <- proc.time()
ptm_cpp_end - ptm_cpp_start

mean((start - estResults$par)^2)
mean((start - Minimum$estimate)^2)

par.est <- estResults$par

bhat <- par.est

spread.hat <- bhat[1:M]
alpha.hat <- bhat[(M+1):(3*M)]
bhat <- bhat[-(1:(3*M))]

C.mat <- matrix(0, nrow=2, ncol=2)
C.mat[lower.tri(C.mat, diag=TRUE)] <- bhat[1:3]
bhat <- bhat[-(1:3)]
Omega.temp <- C.mat %*% t(C.mat)
Omega <- as.matrix(bdiag(Omega.temp))
for (i in 2:M) {
  C.mat <- matrix(0, nrow=2, ncol=2)
  C.mat[lower.tri(C.mat, diag=TRUE)] <- bhat[1:3]
  bhat <- bhat[-(1:3)]
  Omega.temp <- C.mat %*% t(C.mat)
  Omega <- as.matrix(bdiag(Omega,Omega.temp))
}

sigma.hat <- bhat[1]
delta1.hat <- bhat[2]
delta2.hat <- bhat[3]


# price discovery measures (estimated parameters)

tau <- as.vector(tau)

alpha.tau <- matrix(alpha.hat,nrow=2*M) %*% tau^delta2.hat * sigma.hat
sigma.tau <- sigma.hat * tau^delta1.hat
beta.tau <- rep(1,2*M) + alpha.tau / sigma.tau^2

beta <- function(tau) {
  rep(1,2*M) + matrix(alpha.hat,nrow=2*M) %*% tau^(delta2.hat-delta1.hat)
}
# compare market 1 and 2
plot(beta(1:60)[1,], type="l", ylim=c(1,1.1))
lines(beta(1:60)[3,])

gamma.tau <- solve(Omega.hat) %*% beta.tau %*% diag(sigma.tau^2 / (1 + sigma.tau^2 * diag(t(beta.tau) %*% solve(Omega.hat) %*% beta.tau)))

gamma <- function(tau) {
  if (length(tau) == 1) { x <- sigma.hat^2*tau^(2*delta1.hat) * solve(Omega.hat) %*% beta(tau) /
    as.vector(1 + sigma.hat^2*tau^(2*delta1.hat) * t(beta(tau)) %*% solve(Omega.hat) %*% beta(tau)) } else {
      x <- solve(Omega.hat) %*% beta(tau) %*% diag(sigma.hat^2*tau^(2*delta1.hat)) %*%
        solve(diag(1 + sigma.hat^2*tau^(2*delta1.hat) * diag(t(beta(tau)) %*% solve(Omega.hat) %*% beta(tau))))
    }
  return(x)
}
# compare market 1 and 2
plot(gamma(1:60)[1,], type="l", ylim=c(0,0.4))
lines(gamma(1:60)[3,])

# R2.tau: Sum of IS.tau of all markets
R2.tau <- diag(t(beta.tau) %*% gamma.tau)

# information share for each market
IS.tau.mat <- matrix(NA, nrow=M, ncol=J)
for (j in 1:M) {
  IS.tau.mat[j,] <- diag(t(beta.tau[((2*j-1):(2*j)),]) %*% gamma.tau[((2*j-1):(2*j)),])
}

steps <- c(seq(0.01,1,by=0.01),1:60)
IS1.tau <-  diag(t(beta(steps)[1:2,]) %*% gamma(steps)[1:2,])
IS2.tau <-  diag(t(beta(steps)[3:4,]) %*% gamma(steps)[3:4,])
IS3.tau <-  diag(t(beta(steps)[5:6,]) %*% gamma(steps)[5:6,])
IS4.tau <-  diag(t(beta(steps)[7:8,]) %*% gamma(steps)[7:8,])
IS5.tau <-  diag(t(beta(steps)[9:10,]) %*% gamma(steps)[9:10,])

plot(IS1.tau, ylim=c(0,1), xlab="", ylab="", xaxt="n", type="l")
tickat <- c(50,100,130,160)
labels <- c("500ms","1sec","30sec","60sec")
axis(1, at=tickat, labels = labels)
lines(IS2.tau, lty="dashed", col="blue")
lines(IS3.tau, lty="dashed", col="green")
lines(IS4.tau, lty="dashed", col="red")
lines(IS5.tau, lty="dashed", col="orange")
legend("topright", legend=c("N", "T", "K", "P", "Z"), lty=c(1,2,2,2,2), col=c("black","blue","green","red","orange"), cex=0.8)




