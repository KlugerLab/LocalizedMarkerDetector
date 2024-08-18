library(devtools)
library(pkgdown)

setwd("/banach1/ruiqi/local_marker/LocalizedMarkerDetector")

# Create Template
use_readme_rmd()
# use_news_md()
use_vignette("Example")

# Build Package
usethis::use_mit_license("Ruiqi Li")
devtools::build_readme()
Rcpp::compileAttributes()
devtools::document()
# load all function w/o install
devtools::load_all(".") 


# Render
# rmarkdown::render("vignettes/LMD_demo.Rmd")
pkgdown::clean_site()

# if (!dir.exists("docs/articles")) {
#   dir.create("docs/articles", recursive = TRUE)
# }
library(fs)
dir_copy("articles", "docs/articles", overwrite = TRUE)
# file.copy("articles/LMD_cross_comparison_demo.html", "docs/articles/LMD_cross_comparison_demo.html", overwrite = TRUE)

pkgdown::build_site()


# Install Package from local
devtools::install() # install package

# Install Package from github
devtools::install_github("ruiqi0130/LocalizedMarkerDetector")

