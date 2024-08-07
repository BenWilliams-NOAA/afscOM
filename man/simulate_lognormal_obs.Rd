% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/obs_utils.R
\name{simulate_lognormal_obs}
\alias{simulate_lognormal_obs}
\title{Simulate observations from a lognormal distribution
#'
A wrapper function around `rlnorm` that generates a
single random observation from a lognormal distribution
centered on a predicted value (`pred`) and given a
level of error (`cv`).}
\usage{
simulate_lognormal_obs(pred, cv)
}
\arguments{
\item{pred}{the predicted value of an observation}

\item{cv}{the real-space coefficient of variation about
the predicted observation value}
}
\description{
Simulate observations from a lognormal distribution
#'
A wrapper function around `rlnorm` that generates a
single random observation from a lognormal distribution
centered on a predicted value (`pred`) and given a
level of error (`cv`).
}
\examples{
\dontrun{
simulate_lognormal_obs(10, 0.20)
}

}
