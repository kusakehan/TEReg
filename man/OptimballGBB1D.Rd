\name{OptimballGBB1D}
\alias{OptimballGBB1D}
\title{
Estimate the envelope subspace (Feasi 1D)
}
\description{
The 1D algorithm to estimate the envelope subspace with specified dimension based on Wen and Yin (2013).
}
\usage{
OptimballGBB1D(M, U, d, opts=NULL)
}

\arguments{
 \item{M}{M matrix in the envelope objective function. An \eqn{r}-by-\eqn{r} positive semi-definite matrix.}
 \item{U}{U matrix in the envelope objective function. An \eqn{r}-by-\eqn{r} positive semi-definite matrix.}
  \item{d}{Dimension of the envelope. An integer between 0 and r.}
  \item{opts}{Option structure with fields:\cr
            \code{"record = 0"} -- no print out.\cr
            \code{"mxitr"} -- max number of iterations.\cr
            \code{"xtol"} -- stop control for \eqn{||X_k - X_{k-1}||}.\cr
            \code{"gtol"} -- stop control for the projected gradient.\cr
            \code{"ftol"} -- stop control for \eqn{\frac{|F_k - F_{k-1}|}{(1+|F_{k-1}|)}} usually with \code{max{xtol, gtol} > ftol}.\cr

The default values are: \code{"xtol"=1e-08; "gtol"=1e-08; "ftol"=1e-12; "mxitr"=500.}
}
}
\details{
Estimate \code{M}-envelope contains \code{span(U)}
where \code{M > 0} and is symmetric. The
dimension of the envelope is \code{d}.
}
\value{
\item{Ghat}{The orthogonal basis of the envelope subspace.}
}
\references{
Wen, Z., & Yin, W. (2013). A feasible method for optimization with orthogonality constraints. Mathematical Programming, 142(1-2), 397-434.
}
\examples{
##simulate two matrices M and U with an envelope structure#
p <- 20
u <- 5

##randomly generate a semi-orthogonal p-by-u basis matrix (Gamma) for the
##envelope and its orthogonal completion (Gamma0) of dimension p-by-(p-u)
Gamma <- matrix(runif(p*u), p, u)

###make Gamma semi-orthogonal
Gamma <- qr.Q(qr(Gamma))
Gamma0 <- qr.Q(qr(Gamma),complete=TRUE)[,(u+1):p]

## randomly generated symmetric positive definite matrices, M and U, to have
## an exact u-dimensional envelope structure

Phi <- matrix(runif(u^2), u, u)
Phi <- Phi \%*\% t(Phi)
Omega <- matrix(runif(u^2), u, u)
Omega <- Omega \%*\% t(Omega)
Omega0 <- matrix(runif((p-u)^2),p-u,p-u)
Omega0 <- Omega0 \%*\% t(Omega0)
M <- Gamma \%*\% Omega \%*\% t(Gamma) + Gamma0 \%*\% Omega0 \%*\% t(Gamma0)
U <- Gamma \%*\% Phi \%*\% t(Gamma)

# randomly generate symmetric positive definite matrices, Mhat and Uhat, as
# root-n consistent sample estimators for M and U

n=200
X <- mvrnorm(n, mu = rep(0, p), Sigma = M)
Y <- mvrnorm(n, mu = rep(0, p), Sigma = U)
Mhat <- (t(X) \%*\% X)/n
Uhat <- (t(Y) \%*\% Y)/n

Ghat_1D <- OptimballGBB1D(Mhat, Uhat, d=u)
}
