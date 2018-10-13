\name{manifold1D}
\alias{manifold1D}

\title{
Estimate the envelope subspace (ManifoldOptim 1D)
}
\description{
The 1D algorithm to estimate the envelope subspace with specified dimension based on R package "ManifoldOptim".
}
\usage{
manifold1D(M, U, d, params=NULL)
}

\arguments{
  \item{M}{M matrix in the envelope objective function. An \eqn{r}-by-\eqn{r} positive semi-definite matrix.}
 \item{U}{U matrix in the envelope objective function. An \eqn{r}-by-\eqn{r} positive semi-definite matrix.}
  \item{d}{Dimension of the envelope. An integer between 0 and r.}
  \item{params}{Option structure with fields:\cr
            \code{"max_iter"} -- max number of iterations.\cr
            \code{"tol"} -- Tolerance used to assess convergence. See Huang et al (2018) for details on
how this is used.\cr
            \code{"method"} -- Name of optimization method supported by R package \pkg{ManifoldOptim}.
\itemize{
\item{\code{"LRBFGS"}}: Limited-memory RBFGS
\item{\code{"LRTRSR1"}}: Limited-memory RTRSR1
\item{\code{"RBFGS"}}: Riemannian BFGS
\item{\code{"RBroydenFamily"}}: Riemannian Broyden family
\item{\code{"RCG"}}: Riemannian conjugate gradients
\item{\code{"RNewton"}}: Riemannian line-search Newton
\item{\code{"RSD"}}: Riemannian steepest descent
\item{\code{"RTRNewton"}}: Riemannian trust-region Newton
\item{\code{"RTRSD"}}: Riemannian trust-region steepest descent
\item{\code{"RTRSR1"}}: Riemannian trust-region symmetric rank-one update
\item{\code{"RWRBFGS"}}: Riemannian BFGS
}

\code{"check"} -- Should internal manifold object check inputs and print summary message before optimization (TRUE or FALSE).\cr

The default values are: \code{"max_iter"=500; "tol"=1e-08;"method"="RCG"; "check"="False"}.
}
}


\details{
Estimate \code{M}-envelope contains \code{span(U)}
where \code{M > 0} and is symmetric. The
dimension of the envelope is \code{d}.
}
\value{
\item{Gammahat}{The orthogonal basis of the envelope subspace.}
}
\references{
Huang, W., Absil, P. A., Gallivan, K. A., & Hand, P. (2018). ROPTLIB: an object-oriented C++ library for optimization on Riemannian manifolds. ACM Transactions on Mathematical Software (TOMS), 44(4), 43.
}
\note{
\code{ManifoldOptim_module} should be loaded first before using function \code{manifold1D}.
}
\examples{
mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
mani.params <- get.manifold.params(IsCheckParams = TRUE)

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


Ghat_1D <- manifold1D(Mhat, Uhat, u, params=NULL)
}
