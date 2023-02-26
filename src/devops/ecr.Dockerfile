ARG TAG
# This will take an argument regarding which environment to build...
# OPTIONS:
# a) main
# b) train
# c) delft

FROM registry.gitlab.com/fl84inc/comp-bio:${TAG}

COPY . /workspace/comp-bio

WORKDIR /workspace/comp-bio