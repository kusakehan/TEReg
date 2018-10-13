\name{ECD}
\alias{ECD}
\title{
ECD algorithm for estimating the envelope subspace
}
\description{
Estimate the envelope subspace with specified dimension based on ECD algorithm that described in Cook, R. D., & Zhang, X. (2018).
}
\usage{
ECD(M, U, d, maxiter=500, epsilon=1e-08)
}

\arguments{
 \item{M}{M matrix in the envelope objective function. An \eqn{r}-by-\eqn{r} positive semi-definite matrix.}
 \item{U}{U matrix in the envelope objective function. An \eqn{r}-by-\eqn{r} positive semi-definite matrix.}
 \item{d}{Envelope dimension. An integer between 0 and \eqn{r}.}
 \item{maxiter}{Maximum number of iterations.}
 \item{epsilon}{Convergence criterion. \eqn{|F_k - F_{k-1}|< \epsilon}, where \eqn{F_k}   is the objective function.}
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
Cook, R. D., & Zhang, X. (2018). Fast envelope algorithms. Statistica Sinica, 28(3), 1179-1197.
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

Ghat_ECD <- ECD(Mhat, Uhat, u)
}
