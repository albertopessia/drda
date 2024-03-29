% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logistic5.R
\name{logistic5_gradient}
\alias{logistic5_gradient}
\alias{logistic5_hessian}
\alias{logistic5_gradient_hessian}
\title{5-parameter logistic function gradient and Hessian}
\usage{
logistic5_gradient(x, theta)

logistic5_hessian(x, theta)

logistic5_gradient_hessian(x, theta)
}
\arguments{
\item{x}{numeric vector at which the function is to be evaluated.}

\item{theta}{numeric vector with the five parameters in the form
\code{c(alpha, delta, eta, phi, nu)}.}
}
\value{
Gradient or Hessian evaluated at the specified point.
}
\description{
Evaluate at a particular set of parameters the gradient and Hessian of the
5-parameter logistic function.
}
\details{
The 5-parameter logistic function \verb{f(x; theta)} is defined here as

\verb{g(x; theta) = 1 / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)}
\verb{f(x; theta) = alpha + delta g(x; theta)}

where \code{theta = c(alpha, delta, eta, phi, nu)}, \code{eta > 0}, and \code{nu > 0}. When
\code{delta} is positive (negative) the curve is monotonically increasing
(decreasing).
}
