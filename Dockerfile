FROM ubuntu:bionic

## Set a default user. Available via runtime flag `--user docker` 
## Add user to 'staff' group, granting them write privileges to /usr/local/lib/R/site.library
## User should also have & own a home directory (for rstudio or linked volumes to work properly). 
RUN useradd docker \
	&& mkdir /home/docker \
	&& chown docker:docker /home/docker \
	&& addgroup docker staff

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		software-properties-common \
                ed \
		less \
		locales \
		vim-tiny \
		wget \
		ca-certificates \
        && add-apt-repository --enable-source --yes "ppa:marutter/rrutter3.5" \
	&& add-apt-repository --enable-source --yes "ppa:marutter/c2d4u3.5" 

## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

## This was not needed before but we need it now
ENV DEBIAN_FRONTEND noninteractive

# Now install R and littler, and create a link for littler in /usr/local/bin
# Default CRAN repo is now set by R itself, and littler knows about it too
# r-cran-docopt is not currently in c2d4u so we install from source
RUN apt-get update \
        && apt-get install -y --no-install-recommends \
                 littler \
 		 r-base \
 		 r-base-dev \
 		 r-recommended \
  	&& ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
 	&& ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
 	&& ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
 	&& ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
 	&& install.r docopt \
 	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
 	&& rm -rf /var/lib/apt/lists/*

CMD ["bash"]

## Install rJava
RUN apt-get -y update && apt-get install -y \
   default-jdk  r-cran-rjava  r-cran-nloptr libssh2-1-dev

RUN yes | apt-get install libv8-dev
RUN yes | apt-get install libxml2-dev
RUN yes | apt-get install libudunits2-dev
RUN yes | apt install libgdal-dev
RUN yes | apt-get install libmagick++-dev

## Install extra R packages using requirements.R
RUN R -e "install.packages('BiocManager'); library('BiocManager')"
RUN R -e "BiocManager::install('Rsamtools')"
RUN R -e "BiocManager::install('GenomicAlignments')"
RUN R -e "install.packages('Biostrings',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('RCurl')"
RUN R -e "BiocManager::install('doParallel')"
RUN R -e "BiocManager::install('ShortRead')"
RUN R -e "BiocManager::install('tidyverse')"

# Create work directory
WORKDIR /temp

# CLEANUP DOCKER FILE SYSTEM
RUN apt-get -qy autoremove

CMD ["/bin/bash"]