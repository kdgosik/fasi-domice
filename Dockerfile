FROM rocker/tidyverse

USER root

RUN apt-get update && apt-get install -y libglpk40

RUN wget https://ftp.gnu.org/pub/gnu/libiconv/libiconv-1.15.tar.gz && \
tar -zxvf libiconv-1.15.tar.gz && \
cd libiconv-1.15 && \
./configure --prefix=/usr/local && \
make && \
make install

ENV export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

COPY ./install.R /usr/local/src/install.R

RUN Rscript /usr/local/src/install.R

