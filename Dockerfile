FROM rocker/rstudio:4.1.2

LABEL org.opencontainers.image.licenses="GPL-2.0-or-later" \
      org.opencontainers.image.source="https://github.com/rocker-org/rocker-versioned2" \
      org.opencontainers.image.vendor="Rocker Project" \
      org.opencontainers.image.authors="Carl Boettiger <cboettig@ropensci.org>"

RUN echo `pwd`

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libxml2-dev \
    libcairo2-dev \
    libgit2-dev \
    default-libmysqlclient-dev \
    libpq-dev \
    libsasl2-dev \
    libsqlite3-dev \
    libssh2-1-dev \
    libxtst6 \
    libcurl4-openssl-dev \
    xml2 \
    openssl \
    libhdf5-dev \
    libfftw3-dev \
    curl \
    unixodbc-dev

RUN rm -rf /var/lib/apt/lists/*

RUN install2.r --error --skipinstalled \
    tidyverse \
    devtools \
    rmarkdown \
    gert \
    reshape2 \
    changepoint \
    lme4 \
    Formula

RUN R -e 'install.packages("httr")'
RUN R -e 'install.packages("remotes")'
RUN R -e 'options(timeout=9999999)'
RUN R -e 'options(download.file.method = "libcurl")'
RUN R -e 'Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")'
RUN R -e 'remotes::install_github("SoundAg/sounDMR", ref="main")'

RUN rm -rf /tmp/downloaded_packages

