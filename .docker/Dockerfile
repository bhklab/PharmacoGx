FROM rocker/hadleyverse:latest
MAINTAINER Petr Smirnov <psmirnov2000@gmail.com>

RUN r -e 'source("https://raw.githubusercontent.com/MangoTheCat/remotes/master/install-github.R")$value("mangothecat/remotes")' \
  && r -e 'remotes::install_github("bhklab/PharmacoGx@release")' \
  && rm -rf /tmp/downloaded_packages/ /tmp/*.rds
