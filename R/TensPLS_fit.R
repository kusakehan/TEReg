#' @export
#' @import rTensor
#' @importFrom pracma kron sqrtm
#' @importFrom stats cov
TensPLS_fit <- function(Xn, Yn, SigX, u) {
  ss <- dim(Xn)
  len <- length(ss)
  n <- ss[len]
  p <- ss[1:(len-1)]
  m <- length(p)
  SigY <- (n-1)*cov(t(Yn))/n
  
  Sinvhalf <- NULL
  for (i in 1:m) {
    Sinvhalf[[i]] <- pracma::sqrtm(SigX[[i]])$Binv
  }
  
  Sinvhalf[[m+1]] <- pracma::sqrtm(SigY)$Binv
  
  C <- rTensor::ttm(Xn, Yn, m+1)/n
  
  Gamma <- PGamma <- NULL
  for (i in 1:m) {
    M <- SigX[[i]]
    idx <- c(1:(m+1))[-i]
    
    Ck <- rTensor::ttl(C, Sinvhalf[idx], ms = idx)
    
    U <- rTensor::unfold(Ck, row_idx = i, col_idx = idx)@data
    Uk <- U %*% t(U)
    
    Gamma[[i]] <- EnvMU(M, Uk, u[i])
    tmp3 <- t(Gamma[[i]]) %*% SigX[[i]] %*% Gamma[[i]]
    PGamma[[i]] <- Gamma[[i]] %*% solve(tmp3) %*% t(Gamma[[i]]) %*% SigX[[i]]
  }
  
  return(list(Gamma = Gamma, PGamma = PGamma))
}
