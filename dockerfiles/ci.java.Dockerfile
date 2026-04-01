FROM ghcr.io/lballabio/quantlib-swig-devenv:threadsafe

LABEL org.opencontainers.image.authors="Luigi Ballabio <luigi.ballabio@gmail.com>"
LABEL description="An environment for QuantLib-SWIG CI builds on Linux"

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y default-jdk \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

