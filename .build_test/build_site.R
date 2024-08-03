require(devtools)
require(pkgdown)

setwd("../")

use_readme_rmd()
use_news_md()
#use_vignette("Example") 
devtools::document()

pkgdown::clean_site()
pkgdown::build_site()

