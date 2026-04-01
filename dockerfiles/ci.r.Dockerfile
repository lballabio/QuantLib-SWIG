FROM ghcr.io/lballabio/quantlib-swig-devenv:default

LABEL org.opencontainers.image.authors="Luigi Ballabio <luigi.ballabio@gmail.com>"
LABEL description="An environment for QuantLib-SWIG CI builds on Linux"

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y r-base-dev texlive \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

