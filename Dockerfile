# Base image https://hub.docker.com/u/rocker/
FROM rocker/shiny:4.2.1

# system libraries of general use
## install debian packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    liblzma-dev \
    libbz2-dev \
    libicu-dev \
	libglpk-dev \
    git \
    systemd
		

## copy necessary files
# delete existing Shiny server examples
RUN rm -rf /srv/shiny-server/*
# copy app folder
COPY / /srv/shiny-server/

# renv.lock file
COPY /renv.lock ./renv.lock

# custom config file
COPY /shiny-server.conf /etc/shiny-server/shiny-server.conf

WORKDIR /srv/shiny-server/

# # install renv & restore packages
# RUN Rscript -e 'install.packages("renv")'
# RUN Rscript -e 'renv::restore()'

## version without renv, needs manual package list, since packages will be installed manually
RUN rm -rf renv*
RUN rm .Rprofile
## install packages manually
RUN bash packages_install.sh


# make sure remotes package is installed
RUN install2.r --error --deps TRUE remotes >> log.txt 2>&1

# install yenpathy from github with specific version
RUN R -e "remotes::install_github('ecohealthalliance/yenpathy@fe19461', build_vignettes = FALSE)"

# expose port
EXPOSE 3838

# run app on container start
CMD ["/usr/bin/shiny-server"]