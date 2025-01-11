## Overview
Repository for figure and code tracking for fasi-domice project.  All relevant code to generate figures are in the figures directory.

## Quickstart Environment Setup

  1. Run docker command `docker run -it -v $PWD:/workspace/fasi-domice -p 8787:8787 kdgosik/main-gitpod bash`
  2. in the vs-code terminal run the command `bin/rstudio.sh` to start the Rstudio environment
  3. enter the set user name: `gitpod` and password: `gitpod`

## Setting up the development environment:
To build image and run from scratch:

  - Install docker
    - Build the docker image, docker build -t kdgosik/main-gitpod:latest .
    - This takes 10-15 mins to build
    - Launch the container to go into mgcpy's dev env, docker run -it --rm --name mgcpy-env mgcpy:latest
  - Pull image from Dockerhub and run:
    - docker pull kdgosik/main-gitpod:latest
    - `docker run -it -v $PWD:/workspace/fasi-domice -p 8787:8787 kdgosik/main-gitpod bash`
  - To run demo Rstudio (from within Docker):
    - `bin/rstudio.sh`
    - enter the set user name: `gitpod` and password: `gitpod`
    - got to port 8787 in you browser to access rstudio
  - Then open download-files.R
