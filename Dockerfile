# Base dockerfile: rocker rstudio
FROM rocker/rstudio:4.1.3

# Set global R options
RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site
ENV RETICULATE_MINICONDA_ENABLED=FALSE


# System libraries required for common R packages
RUN apt-get update && apt-get install -y \
    cmake \
    git \
    libboost-all-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    libfftw3-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libglpk-dev \
    libgmp-dev \
    libgsl-dev \
    libharfbuzz-dev \
    libhdf5-dev \
    libjpeg-dev \
    liblzma-dev \
    libpng-dev \
    libssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libtiff5-dev \
    libxml2-dev \
    libz-dev \
    libzip-dev \
    llvm-10 \
    openjdk-8-jdk \
    pandoc \
    python3-dev \
    python3-pip \
    wget \
    zlib1g-dev \
    libgit2-dev \
    && apt-get clean

RUN R --no-echo --no-restore --no-save -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

# Set working directory inside container
WORKDIR /project
COPY renv.lock renv.lock

RUN mkdir -p renv
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

# Install renv and restore dependencies
RUN R --no-echo --no-restore --no-save -e "renv::restore(confirm = FALSE)"
RUN R --no-echo --no-restore --no-save -e "renv::activate()"

# Optionally: install pkgdown and build documentation
# (Uncomment the lines below if you want the website rendered in the image)
# RUN R -e "install.packages('pkgdown'); pkgdown::build_site()"

# Expose RStudio Server
EXPOSE 8787

# Default command (inherited from rocker/rstudio)
CMD ["/init"]