% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loglogistic4.R
\name{loglogistic4_fn}
\alias{loglogistic4_fn}
\title{4-parameter log-logistic function}
\usage{
loglogistic4_fn(x, theta)
}
\arguments{
\item{x}{numeric vector at which the function is to be evaluated.}

\item{theta}{numeric vector with the four parameters in the form
\code{c(alpha, delta, eta, phi)}.}
}
\value{
Numeric vector of the same length of \code{x} with the values of the
log-logistic function.
}
\description{
Evaluate at a particular set of parameters the 4-parameter log-logistic
function.
}
\details{
The 4-parameter log-logistic function \verb{f(x; theta)} is defined here as

\verb{g(x; theta) = x^eta / (x^eta + phi^eta)}
\verb{f(x; theta) = alpha + delta g(x; theta)}

where \code{x >= 0}, \code{theta = c(alpha, delta, eta, phi)}, \code{eta > 0}, and
\code{phi > 0}.
}
