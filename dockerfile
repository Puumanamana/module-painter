FROM mambaorg/micromamba:git-e2d1bf8-focal-cuda-11.3.1
LABEL author="carisdak@hawaii.edu"

COPY . .

RUN micromamba update --yes --file conda.yaml -n base

WORKDIR /workspace
