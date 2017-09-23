#!/bin/sh

./configure CXX=clang++ CXXFLAGS='-g -O0 -pedantic -Wall -std=c++11 -I/opt/X11/include -I/opt/local/include' LDFLAGS='-stdlib=libc++ -L/opt/local/lib -L/opt/X11/lib' --prefix=/opt/local/ --no-create --no-recursion
