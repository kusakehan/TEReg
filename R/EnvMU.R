#' @export
EnvMU <- function(M, U, m) {
  p <- dim(U)[1]
  W <- matrix(0, p, (m+1))
  for (k in 1:m) {
    Wk <- W[, 1:k]
    Ek <- M %*% Wk
    temp <- t(Ek) %*% Ek
    QEK <- diag(p) - Ek %*% MASS::ginv(temp) %*% t(Ek)
    W[, (k+1)] <- eigen(QEK %*% U %*% QEK)$vectors[, 1]
  }
  Gamma <- qr.Q(qr(Re(W[, 2:(m+1)])))
  return(Gamma)
}