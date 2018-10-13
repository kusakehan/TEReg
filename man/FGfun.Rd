\name{FGfun}
\alias{FGfun}
\title{
The Objective function and its gradient}
\description{
The objective function and its gradient value that defined in equation (3.1) of Cook, R. D., & Zhang, X. (2016).\cr
}
\usage{
FGfun(W, M, U)
}
\arguments{
\item{M}{M matrix in the envelope objective function. An \eqn{r}-by-\eqn{r} positive semi-definite matrix.}
 \item{U}{U matrix in the envelope objective function. An \eqn{r}-by-\eqn{r} positive semi-definite matrix.}
  \item{W}{A vector of \eqn{r} by \eqn{d}.}
}
\details{
This is the objective function and its gradient for estimating \code{M}-envelope contains \code{span(U)}, where \code{M > 0} and is symmetric, the
dimension of the envelope is \eqn{d}.
}
\value{
\item{F}{The value of objective function given \code{W}.}
\item{G}{The value of the gradient function given \code{W}.}
}
\references{
Cook, R. D., & Zhang, X. (2016). Algorithms for envelope estimation. Journal of Computational and Graphical Statistics, 25(1), 284-300.
}
