% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate_dem_params.R
\name{validate_dem_params}
\alias{validate_dem_params}
\title{Validate Demographic Parameter Matrices}
\usage{
validate_dem_params(dem_params, model_dims)
}
\arguments{
\item{dem_params}{a named list containing the demographic rate matrices.
List elements must have the same names as stated above.}

\item{model_dims}{a named list containing expected model dimensions. List
element must be: nyears, nages, nsexes, nregions, and nfleets.}
}
\value{
unmodified dem_params object; if dem_params list is invalid, will
return a descriptive error message
}
\description{
Helper function for validating that all required demographic
parameter matrices have been defined and that they are all of
the correct dimension.
}
\details{
Required demographic parameters are: mort, waa, sexrat, sel,
ret, dmr, and surv_sel

Expected dimensions are: [nyears, nages, nsexes, nregions, nfleets]
}
\examples{
\dontrun{

     model_dims = list(nyears=64, nages=30, nsexes=2, nregions=1, nfleets=2)
     dem_params = list(
         mort = mortality_matrix,
         waa = weight_matrix,
         sexrat = sexratio_matrix,
         sel = selectivity_matrix,
         ret = retention_matrix
         dmr = disard_mortality_matrix,
         surv_sel = survey_selectivity_matrix
     )
     dem_params = validate_dem_params(dem_params, model_dims)

}

}
