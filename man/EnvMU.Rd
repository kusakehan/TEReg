\name{EnvMU}
\alias{EnvMU}

\title{
Estimate envelope subspace basis
}
\description{
SIMPLS-type algorithm for estimating the \code{M}-envelope contains \code{span(U)} without manifold optimization. This function is a generalization of De Jong, S. (1993) and Cook, R. D., Helland, I. S., & Su, Z. (2013).
}
\usage{
EnvMU(M, U, m)
}
\arguments{
 \item{M}{M matrix in the envelope objective function. An \eqn{r}-by-\eqn{r} positive semi-definite matrix.}
 \item{U}{U matrix in the envelope objective function. An \eqn{r}-by-\eqn{r} positive semi-definite matrix.}
 \item{m}{The envelope dimension.}
}

\value{
\item{Gamma}{The basis of \code{M}-envelope contains \code{span(U)}.}
}
\references{
Cook, R. D., Helland, I. S., & Su, Z. (2013). Envelopes and partial least squares regression. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 75(5), 851-877.

De Jong, S. (1993). SIMPLS: an alternative approach to partial least squares regression. Chemometrics and intelligent laboratory systems, 18(3), 251-263.
}
