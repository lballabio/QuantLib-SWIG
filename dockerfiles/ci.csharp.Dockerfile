FROM ghcr.io/lballabio/quantlib-swig-devenv:threadsafe

LABEL org.opencontainers.image.authors="Luigi Ballabio <luigi.ballabio@gmail.com>"
LABEL description="An environment for QuantLib-SWIG CI builds on Linux"

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y libicu76 \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN cd /tmp \
 && wget https://dot.net/v1/dotnet-install.sh \
 && bash dotnet-install.sh --install-dir /usr/local/bin/ -c 9.0
