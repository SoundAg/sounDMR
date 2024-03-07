FROM rocker/tidyverse:latest as base

# Add labels???
LABEL org.opencontainers.image.licenses="GPL-2.0-or-later" \
      org.opencontainers.image.source="https://github.com/rocker-org/rocker-versioned2" \
      org.opencontainers.image.vendor="Rocker Project" \
      org.opencontainers.image.authors="Carl Boettiger <cboettig@ropensci.org>"

# Install bootstrap library.
RUN install2.r --error devtools

# Install the SounDMR package.
RUN R -e 'devtools::install_github("jokergoo/GetoptLong")'
RUN R -e 'devtools::install_github("SoundAg/sounDMR")'
COPY splitbychrw.R /opt/
COPY splitbychunkw.R /opt/
COPY getmethylframemultiplew.R /opt/
COPY getmethylframew.R /opt/

# Install pandas
RUN apt-get update && apt-get install -y --no-install-recommends build-essential python3.6 python3-pip python3-setuptools python3-dev
RUN pip3 install pandas
COPY sbchunk.py /opt/
COPY splitbychromosome.py /opt/

