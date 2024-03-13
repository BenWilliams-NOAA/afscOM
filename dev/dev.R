# dev code for afscOM

# load ----
library(devtools)
library(pkgdown)

# ignore the development folder from package building
usethis::use_build_ignore("dev")

# dev code
use_package('dplyr')
use_package('TMB')
use_package('abind')
use_package('rlang')
use_package('bbmle')

use_r('data')
use_r('simulate')
# create the package
devtools::document()
devtools::build()

use_mit_license()

usethis::use_vignette('intro')
# add functions
# use_r('fish_age_comp')
# use_r('slx_a50_a95')



remotes::install_github("BenWilliams-NOAA/afscOM")
# create website
# this will create a 'docs' folder where vignettes etc for the associated github website will be held
pkgdown::build_site(examples = FALSE)
