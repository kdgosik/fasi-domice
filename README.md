## Overview
Repository for figure and code tracking for "Prevalent cross-cell type QTL trans-regulatory genetic effects impacting innate lymphoid cells in the small intestine".  All relevant code to generate figures are in the figures directory.

## Quickstart Environment Setup

  1. Run docker command `docker run -it -v $PWD:/workspace/fasi-domice -p 8787:8787 kdgosik/main-gitpod bash`
  2. in the vs-code terminal run the command `bin/rstudio.sh` to start the Rstudio environment
  3. enter the set user name: `gitpod` and password: `gitpod`

## Setting up the development environment:
To build image and run from scratch:

  - Install docker
    - Build the docker image, docker build -t kdgosik/main-gitpod:latest .
    - This takes 20-30 mins to build
    - Launch the container using `docker run -it -v $PWD:/workspace/fasi-domice -p 8787:8787 kdgosik/main-gitpod bash`
  - Pull image from Dockerhub and run:
    - docker pull kdgosik/main-gitpod:latest
    - `docker run -it -v $PWD:/workspace/fasi-domice -p 8787:8787 kdgosik/main-gitpod bash`
  - To run demo Rstudio (from within Docker):
    - `cd /workspace/fasi-domice`
    - `bin/rstudio.sh`
    - enter the set user name: `gitpod` and password: `gitpod`
    - got to port 8787 in you browser to access rstudio
  - Then open download-files.R
    - run `setup.R` to create necessary directory locations
    - run `download-files.R` to pull files from the google drive.
