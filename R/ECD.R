##################################################
#         ECD objective function                 #
##################################################
objF_ECD <- function(A, B, w) {
  fk <- log(t(w) %*% A %*% w) + log(t(w) %*% B%*% w) - 2*log(crossprod(w))
  return(fk)
}

##################################################
#    get initial value for ECD algorithm         #
##################################################
ECDini <- function(M, U) {
  # estimating M-envelope contains span(U)
  # where M>0 and is symmetric
  # dimension of the envelope is d
  # based on (M+U) and inv(M)
  p <- dim(U)[2]
  eigM <- eigen(M); eigMU <- eigen(M + U)
  v1 <- eigM$vectors 
  v2 <- eigMU$vectors 
  v <- cbind(v1, v2)
  W0 <- v[, 1]
  Fw0 <- log(t(W0) %*% solve(M+U) %*% W0) + log(t(W0) %*% M %*% W0)
  for (i in 2:(2*p)) {
    W <- v[, i]
    Fw <- log(t(W) %*% solve(M+U) %*% W) + log(t(W) %*% M %*% W)
    if (Fw < Fw0){
      W0 <- W
      Fw0 <- Fw
    }
  }
  return(W0)
}

##################################################
#     optimECD algorithm for solving fk          #
##################################################
optimECD <- function(A, B, w0, maxiter=500, epsilon=1e-08) {
  p <- length(w0)
#  epsilon <- 1e-08
  eigA <- eigen (A + t(A));
  ## already in descending order###
  Gp <- eigA$vectors; dn <- eigA$values
  dn <- diag(dn/2)

  v0 <- crossprod(Gp, w0)
  GBG <- t(Gp) %*% B %*% Gp
  fk <- objF_ECD(dn, GBG, v0)
  v <- v0
  for (iter in 1:maxiter) {
    flg <- 0
    alpha <- 1/(t(v) %*% dn %*% v)
    beta <- 1/(t(v) %*% GBG %*% v)
    delta <- 1/(t(v) %*% v)
    A1 <- as.numeric(alpha)*dn
    B1 <- as.numeric(beta)*GBG
    for (j in 1:p) {
      AB1 <- A1[j, j] + B1[j, j]
      if ((2*delta - AB1) != 0) {
        v[j] <- (crossprod(v, A1[,j]) + crossprod(v, B1[, j]) - AB1*v[j])/(2*delta - AB1)
        flg <- flg + 1
      }
      if (objF_ECD(dn, GBG, v) > (objF_ECD(dn, GBG, v0)) + epsilon) {
        v <- v0
        flg <- flg - 1
      }
    }
    fk1 <- objF_ECD(dn, GBG, v)
    if ((abs (fk - fk1)) < epsilon) break;
    fk <- fk1
  }
  w <- Gp %*% v

  w <- w/norm(w, type = "2")
  return(w)
}


##################################################
#   estimating M-envelope contains span(U)       #
##################################################
ECD1st <- function (M, U, maxiter=500, epsilon=1e-08){
  # estimating M-envelope contains span(U)
  # where M>0 and is symmetric
  # dimension of the envelope is d
  # based on inv(M+U) and (M)
  gamma <- optimECD(M, solve(M + U), ECDini(M, U), maxiter, epsilon=1e-08)
  return(gamma)
}



##################################################
#           ECD algorithm                        #
##################################################
#'@export
ECD <- function(M, U, d, maxiter=500, epsilon=1e-08){
  # estimating M-envelope contains span(U)
  # where M>0 and is symmetric
  # dimension of the envelope is d
  # based on inv(M+U) and (M)
  if ( is.null(maxiter)) {
    maxiter <- 2000
  }
  p <- dim(M)[2]
  Mnew <- M
  Unew <- U
  G <- matrix(0, p, d)
  G0 <- diag(p)
  for (k in 1:d) {
    gk <- ECD1st(Mnew, Unew, maxiter, epsilon)
    G[, k]<- G0 %*% gk
    G0 <- qr.Q(qr(G[, 1:k]),complete=TRUE)[,(k+1):p]
    Mnew <- t(G0) %*% M %*% G0
    Unew <- t(G0) %*% U %*% G0
  }
  Ghat <- G
  return(Ghat)
}
