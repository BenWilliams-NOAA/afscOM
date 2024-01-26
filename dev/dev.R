# dev code for afscOM

# load ----
library(devtools)
library(pkgdown)


usethis::use_build_ignore("dev")
# dev code
use_package('dplyr')
use_package('TMB')

use_mit_license()

# add functions
use_r('fish_age_comp')
