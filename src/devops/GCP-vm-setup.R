
# 8 CPU - 32 GB Ram
# n2-standard-8

library(googleComputeEngineR)

## setting up a 13GB RAM instance 
## see gce_list_machinetype() for options of predefined_type
vm <- gce_vm(template = "rstudio",
             name = "rstudio-team",
             username = "kirk", password = "kirk",
             predefined_type = "n2-standard-8")

## wait a bit, login at the IP it gives you

# exec into container to install library
# docker exec -it 6fde45faef0a bash
# 
# update
# apt-get update
# 
# install necessary library
# sudo apt-get install libbz2-dev