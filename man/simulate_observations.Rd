% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_observations.R
\name{simulate_observations}
\alias{simulate_observations}
\title{Simulate Observations from Fisheries and Surveys}
\usage{
simulate_observations(naa, waa, selex, faa, zaa, caa, obs_pars, age_error = NA)
}
\arguments{
\item{naa}{numbers-at-age array subset to 1 year (dim [1, nages, nsexes, nregions])}

\item{waa}{weight-at-age array subset to 1 year (dim [1, nages, nsexes, nregions])}

\item{selex}{selectivity-at-age array subset to 1 year (dim [1, nages, nsexes, nregions, nfleets])}

\item{faa}{fishing mortality-at-age array subset to 1 year (dim [1, nages, nsexes, nregions, nfleets])}

\item{zaa}{total mortality-at-age array subset to 1 year (dim [1, nages, nsexes, nregions])}

\item{caa}{total catch-at-age array subset to 1 year (dim [1, nages, nsexes, nregions, nfleets])}

\item{obs_pars}{list of observation process parameters}

\item{age_error}{an optional ageing error matrix}
}
\description{
Simulate observations of a population from common Alaska surveys and
from landed catch data.
}
\examples{
\dontrun{
simulate_observations(naa, waa, selex, faa, zaa, caa, obs_pars, age_error=NA)
}

}
