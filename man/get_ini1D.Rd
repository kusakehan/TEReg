\name{get_ini1D}
\alias{get_ini1D}
\title{
Initial value for 1D algorithm
}
\description{
This function gives the initial value for 1D algorithm, which is searched from the eigenvalues of \code{M} and \code{M+U}.
}
\usage{
get_ini1D(M, U)
}

\arguments{
 \item{M}{M matrix in the envelope objective function. An \eqn{r}-by-\eqn{r} positive semi-definite matrix.}
 \item{U}{U matrix in the envelope objective function. An \eqn{r}-by-\eqn{r} positive semi-definite matrix.}
}

\value{
\item{W0}{Initial value for the 1D algorihtm}
}

