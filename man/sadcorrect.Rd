% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sadcorrect.R
\name{sadcorrect}
\alias{sadcorrect}
\title{correct structure errors}
\usage{
sadcorrect(x, xmin = 0.1, log = TRUE, rsm = 0, Nx = NULL,
  Ny = NULL, J = NULL, boundaries = "pad", direction = TRUE)
}
\arguments{
\item{x}{a list of equally sized matrices, the first element is assumed to be the observation}

\item{xmin}{values smaller than \code{xmin} are set to zero}

\item{log}{logical, do you want to log-transfrom the data? (recommended for precipitation)}

\item{rsm}{number of pixels which are linearly smoothed at the edge}

\item{Nx}{size to which the data is extended in x-direction, has to be a whole power of 2}

\item{Ny}{size to which the data is extended in y-direction, has to be a whole power of 2}

\item{J}{largest scale considered}

\item{boundaries}{how to handle the boundary conditions, either "pad", "mirror" or "periodic"}

\item{direction}{if \code{TRUE}, scale and direction are corrected, otherwise only scale}
}
\value{
an object of class \code{sadforecast}
}
\description{
use the inverse 'dtcwt' to correct errors in scale, anisotropy and direction
}
\details{
The algorithm performs the following steps:
\enumerate{
   \item remove values below \code{xmin}
   \item if \code{log=TRUE} log-transform all fields
   \item  set all fields to zero mean, unit variance
   \item apply \code{dtcwt} to all fields
   \item loop over forecasts and scales. If \code{direction=TRUE} loop over the six directions. Multiply forecast energy at each location by the ratio of total observed energy to total forecast energy at that scale (and possibly direction)
   \item apply \code{idtcwt} to all forecasts
   \item reset means and variance of the forecasts to their original values
   \item if \code{log=TRUE} invert the log-transform 
   \item return the list of corrected fields
}
}
\examples{
data(rrain)
ra   <- as.sadforecast( list( rrain[2,1,,], rrain[3,1,,], rrain[3,2,,], rrain[3,3,,] ) )
ra_c <- sadcorrect( ra, rsm=10 )
plot(ra_c)
}
