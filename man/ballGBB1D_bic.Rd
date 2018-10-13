\name{ballGBB1D_bic}
\alias{ballGBB1D_bic}
\title{
Envelope dimension selection based on 1D-BIC
}
\description{
This function selects envelope subspace dimension use 1D-BIC by Zhang, X., & Mai, Q. (2018). The constrained optimization in the 1D algorithm is based on the method proposed by Wen and Yin (2013).
}
\usage{
ballGBB1D_bic(M, U, n, multiD=1, bic_max=10, opts=NULL)
}
\arguments{
\item{M}{M matrix in the envelope objective function. An \eqn{r}-by-\eqn{r} positive semi-definite matrix.}
 \item{U}{U matrix in the envelope objective function. An \eqn{r}-by-\eqn{r} positive semi-definite matrix.}
 \item{n}{Sample size.}
 \item{multiD}{A constant, see Zhang, X., & Mai, Q. (2018), the default value is 1.}
 \item{bic_max}{The maximum dimension to consider, \code{bic_max} is smaller than r, the default value is 10.}
 \item{opts}{Option structure for GBB algorithm. See function \code{OptStiefelGBB}.}
}

\value{
\item{bicval}{The BIC values for different envelope dimensions.}
\item{u}{The dimension selected which corresponds to the smallest BIC values.}
}
\references{
Zhang, X., & Mai, Q. (2018). Model-free envelope dimension selection. Electronic Journal of Statistics, 12(2), 2193-2216.\cr
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

res <- ballGBB1D_bic(Mhat, Uhat, n)

## visualization
plot(1:10, res$bicval, type="o", xlab="Envelope Dimension", ylab="BIC values",
main="Envelope Dimension Selection")
}
