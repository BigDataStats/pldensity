% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ddpn_eval}
\alias{ddpn_eval}
\title{Model Training}
\usage{
ddpn_eval(model, x, nparticles = 50L)
}
\arguments{
\item{model}{Object of class DDPN}

\item{x}{nxd matrix where each row is an observation}

\item{epochs}{int >=1. If it's larger than one, the data will be recycled. Defaults to 1.}

\item{new_period_every}{int >= 0. When set to zero, all the data corresponds to same time period, and the BAR occurs at the beginning. Otherwise, a BAR update will happen in the frequency specified.}
}
\description{
Each a DDPN model and new data to learn the distribution.
}
\references{
Matthew A. Taddy (2012) Autoregressive Mixture Models for Dynamic Spatial Poisson Processes:
 Application to Tracking Intensity of Violent Crime, Journal of the American Statistical Association,
 105:492, 1403-1417, DOI: 10.1198/jasa.2010.ap09655
}
