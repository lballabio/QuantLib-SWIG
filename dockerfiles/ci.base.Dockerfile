ARG tag=latest
FROM ghcr.io/lballabio/boost:${tag}
MAINTAINER Luigi Ballabio <luigi.ballabio@gmail.com>
LABEL Description="A development environment for building QuantLib-SWIG on Travis CI"

RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y autoconf automake libtool \
                                                      make libpcre2-dev clang git \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

ENV swig_version=4.1.1

RUN wget http://downloads.sourceforge.net/project/swig/swig/swig-${swig_version}/swig-${swig_version}.tar.gz \
 && tar xfz swig-${swig_version}.tar.gz \
 && rm swig-${swig_version}.tar.gz \
 && cd swig-${swig_version} \
 && ./configure --prefix=/usr CXXFLAGS=-O3 \
 && make -j 4 && make install \
 && cd .. && rm -rf swig-${swig_version}

RUN mkdir /QuantLib
COPY QuantLib /QuantLib

RUN rm /usr/lib/libboost_unit_test_framework.so*

CMD bash
