FROM ghcr.io/lballabio/quantlib-swig-devenv:threadsafe
MAINTAINER Luigi Ballabio <luigi.ballabio@gmail.com>
LABEL Description="A development environment for building QuantLib-SWIG on Travis CI"

RUN cd /tmp \
 && wget https://dot.net/v1/dotnet-install.sh \
 && bash dotnet-install.sh --install-dir /usr/local/bin/ -c 6.0
