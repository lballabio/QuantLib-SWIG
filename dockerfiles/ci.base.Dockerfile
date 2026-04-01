ARG tag=rolling
FROM ghcr.io/lballabio/quantlib-devenv:${tag}

LABEL org.opencontainers.image.authors="Luigi Ballabio <luigi.ballabio@gmail.com>"
LABEL description="An environment for QuantLib-SWIG CI builds on Linux"

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y libpcre2-dev \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

ENV swig_version=4.4.1

RUN wget http://downloads.sourceforge.net/project/swig/swig/swig-${swig_version}/swig-${swig_version}.tar.gz \
 && tar xfz swig-${swig_version}.tar.gz \
 && rm swig-${swig_version}.tar.gz \
 && cd swig-${swig_version} \
 && ./configure --prefix=/usr CXXFLAGS=-O3 \
 && make -j 4 && make install \
 && cd .. && rm -rf swig-${swig_version}

RUN mkdir /QuantLib
COPY QuantLib /QuantLib

CMD bash
