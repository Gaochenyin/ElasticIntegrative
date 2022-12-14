% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenData.R
\name{GenData}
\alias{GenData}
\title{Generate data for integrative analysis}
\usage{
GenData(
  seed = NULL,
  beta0 = c(0, 1, 1, 1),
  psi0 = c(0, 1, 1),
  n = 1e+05,
  mean.x = 1,
  m = 2000,
  tlocalpar = 0,
  n.t = NULL,
  equal.weight = T,
  overlap = T
)
}
\arguments{
\item{seed}{random seed for data generation}

\item{beta0}{the coefficients of the linear function for generating Y(0)}

\item{psi0}{the coefficients of the linear function for generating treatment effect}

\item{n}{the size of finite population}

\item{mean.x}{the mean of X for the finite population}

\item{m}{the size of the real-world data (RWD)}

\item{tlocalpar}{the local bias parameter for local alternative}

\item{n.t}{the size of the randomized clinical trial (RCT)}

\item{equal.weight}{the indicator for whether to assign the treatment equally or not}

\item{overlap}{the indicator for assessing the violation of overlap condition}
}
\description{
Generate data for integrative analysis
}
