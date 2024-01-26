# dev code for afscOM

# load ----
library(devtools)
library(pkgdown)

# ignore the development folder from package building
usethis::use_build_ignore("dev")

# dev code
use_package('dplyr')
use_package('TMB')

use_mit_license()

# add functions
# use_r('fish_age_comp')
# use_r('slx_a50_a95')

# create the package
devtools::document()
devtools::build()

# create website
# this will create a 'docs' folder where vignettes etc for the associated github website will be held
pkgdown::build_site(examples = FALSE)
