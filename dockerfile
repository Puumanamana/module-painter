FROM mambaorg/micromamba:0.24.0
LABEL author="carisdak@hawaii.edu"
ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY --chown=$MAMBA_USER:$MAMBA_USER . .

RUN micromamba install -y -n base -f conda.yaml && \
	micromamba clean --all --yes

RUN rm -rf /tmp/*

WORKDIR /workspace

