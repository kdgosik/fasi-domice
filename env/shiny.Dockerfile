FROM registry.gitlab.com/fl84inc/devenv:base
# to be launched from /workspace/comp-bio/

RUN pip install mord \
    colour \
    pyscenic
    
USER root

RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    xtail \
    wget 

# Install Shiny and Shiny Manager
RUN R -e "install.packages('shiny', repos='http://cran.rstudio.com/')"
RUN R -e 'install.packages("shinymanager")'

COPY shiny/scripts /scripts
RUN /scripts/install_shiny_server.sh


# Copy the essential code
COPY TBR/ /app/TBR/

COPY shiny/shiny-server.conf /etc/shiny-server/shiny-server.conf

RUN mkdir -p /mnt/fsx/computational-bio-data

# Remove example app and populate with tbr.R application
RUN rm -rf /srv/shiny-server/*
COPY shiny/app.R /srv/shiny-server/app.R

EXPOSE 3838

CMD ["/scripts/start.sh"]