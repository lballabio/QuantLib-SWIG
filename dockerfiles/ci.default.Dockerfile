FROM ghcr.io/lballabio/quantlib-swig-devenv:base
MAINTAINER Luigi Ballabio <luigi.ballabio@gmail.com>
LABEL Description="A development environment for building QuantLib-SWIG on Travis CI"

RUN cd /QuantLib \
 && ./autogen.sh \
 && ./configure --enable-throwing-in-cycles --enable-unity-build --disable-static --disable-test-suite --enable-skip-examples CC=clang CXX=clang++ CXXFLAGS='-O3 -g0' \
 && make -j 4 install \
 && cd .. && rm -rf /QuantLib

RUN ldconfig

CMD bash
