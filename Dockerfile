FROM rocker/tidyverse:latest as base

# Add labels???
LABEL org.opencontainers.image.licenses="GPL-2.0-or-later" \
      org.opencontainers.image.source="https://github.com/rocker-org/rocker-versioned2" \
      org.opencontainers.image.vendor="Rocker Project" \
      org.opencontainers.image.authors="Carl Boettiger <cboettig@ropensci.org>"

# Install bootstrap library.
RUN install2.r --error devtools

# Install the SounDMR package.
RUN R -e 'devtools::install_github("SoundAg/sounDMR")'
