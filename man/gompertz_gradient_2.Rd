% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gompertz.R
\name{gompertz_gradient_2}
\alias{gompertz_gradient_2}
\alias{gompertz_hessian_2}
\alias{gompertz_gradient_hessian_2}
\title{Gompertz function gradient and Hessian}
\usage{
gompertz_gradient_2(x, theta)

gompertz_hessian_2(x, theta)

gompertz_gradient_hessian_2(x, theta)
}
\arguments{
\item{x}{numeric vector at which the function is to be evaluated.}

\item{theta}{numeric vector with the four parameters in the form
\code{c(alpha, delta, eta, phi)}.}
}
\value{
Gradient or Hessian of the alternative parameterization evaluated at
the specified point.
}
\description{
Evaluate at a particular set of parameters the gradient and Hessian of the
Gompertz function.
}
\details{
The Gompertz function \verb{f(x; theta)} is defined here as

\verb{g(x; theta) = exp(-exp(-eta * (x - phi)))}
\verb{f(x; theta) = alpha + delta g(x; theta)}

where \code{theta = c(alpha, delta, eta, phi)} and \code{eta > 0}. When \code{delta} is
positive (negative) the curve is monotonically increasing (decreasing).

This set of functions use a different parameterization from
\code{link[drda]{gompertz_gradient}}. To avoid the non-negative
constraints of parameters, the gradient and Hessian computed here are for
the function with \code{eta2 = log(eta)}.

Note that argument \code{theta} is on the original scale and not on the log scale.
}
