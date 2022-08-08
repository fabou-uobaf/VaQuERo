FROM ubuntu:20.04

MAINTAINER famman@cemm.at
ENV LANG en_US.UTF-8
ENV LC_NUMERIC en_US.UTF-8
ENV LC_TIME en_US.UTF-8
ENV LC_MONETARY en_US.UTF-8
ENV LC_PAPER en_US.UTF-8
ENV LC_NAME en_US.UTF-8
ENV LC_ADDRESS en_US.UTF-8
ENV LC_TELEPHONE en_US.UTF-8
ENV LC_MEASUREMENT en_US.UTF-8
ENV LC_IDENTIFICATION en_US.UTF-8
RUN mkdir -p /home/user && cd /home/user
WORKDIR /home/user
ENV HOME /home/user
RUN apt-get update -y
RUN apt-get install -y gnupg2
RUN apt-get install software-properties-common -y
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' -y
RUN apt-get update -y
RUN apt-get install libudunits2-dev libgdal-dev -y
RUN apt-get install git -y
RUN apt-get install r-base r-base-dev -y
RUN git clone https://github.com/fabou-uobaf/VaQuERo.git
RUN Rscript --vanilla VaQuERo/R_package_dependency_install.r
