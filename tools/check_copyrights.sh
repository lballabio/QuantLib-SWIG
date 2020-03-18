#!/bin/bash
#
# Run this from the main QuantLib-SWIG directory

grep -rI --exclude='Makefile*' "Copyright (C)" CSharp Java Python R Scala SWIG | ./tools/collect_copyrights.py

