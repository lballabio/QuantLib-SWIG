# -*- coding: utf-8 -*-
"""
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008 StatPro Italia srl

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <https://www.quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
"""

import os, sys, math, platform, sysconfig
from setuptools import setup, Extension
from setuptools._distutils.ccompiler import get_default_compiler

def is_debug_quantlib():
    return os.getenv("QL_DEBUG", 'False').lower() in ('true', '1', 't')


def define_macros(py_limited_api):

    define_macros = []
    if py_limited_api:
        define_macros += [("Py_LIMITED_API", "0x03080000")]

    compiler = get_default_compiler()

    if compiler == "msvc":
        define_macros += [
            ("__WIN32__", None),
            ("WIN32", None),
            ("NDEBUG", None),
            ("_WINDOWS", None),
            ("NOMINMAX", None),
        ]

    elif compiler == "unix":
        ql_compile_args = os.popen("quantlib-config --cflags").read()[:-1].split()

        define_macros += [
            (arg[2:], None) for arg in ql_compile_args if arg.startswith("-D")
        ]
        define_macros += [("NDEBUG", None)]

    return define_macros


def include_dirs():

    include_dirs = []

    compiler = get_default_compiler()

    if compiler == "msvc":
        try:
            QL_INSTALL_DIR = os.environ["QL_DIR"]
            include_dirs += [QL_INSTALL_DIR]
        except KeyError:
            print("warning: unable to detect QuantLib installation")

        if "INCLUDE" in os.environ:
            include_dirs += [
                d.strip() for d in os.environ["INCLUDE"].split(";") if d.strip()
            ]

    elif compiler == "unix":
        ql_compile_args = os.popen("quantlib-config --cflags").read()[:-1].split()

        include_dirs += [arg[2:] for arg in ql_compile_args if arg.startswith("-I")]

    return include_dirs


def library_dirs():

    library_dirs = []

    compiler = get_default_compiler()

    if compiler == "msvc":
        try:
            QL_INSTALL_DIR = os.environ["QL_DIR"]
            library_dirs += [os.path.join(QL_INSTALL_DIR, "lib")]
        except KeyError:
            print("warning: unable to detect QuantLib installation")

        if "LIB" in os.environ:
            dirs = [dir for dir in os.environ["LIB"].split(";")]
            library_dirs += [d for d in dirs if d.strip()]

    elif compiler == "unix":
        ql_link_args = os.popen("quantlib-config --libs").read()[:-1].split()

        library_dirs += [arg[2:] for arg in ql_link_args if arg.startswith("-L")]

    return library_dirs


def libraries():

    libraries = []

    compiler = get_default_compiler()

    if compiler == "unix":
        ql_link_args = os.popen("quantlib-config --libs").read()[:-1].split()

        libraries += [arg[2:] for arg in ql_link_args if arg.startswith("-l")]

    return libraries


def extra_compile_args():

    extra_compile_args = []

    compiler = get_default_compiler()

    if compiler == "msvc":
        extra_compile_args = ["/GR", "/FD", "/Zm250", "/EHsc", "/bigobj", "/std:c++17"]

        if is_debug_quantlib():
            if "QL_STATIC_RUNTIME" in os.environ:
                extra_compile_args.append("/MTd")
            else:
                extra_compile_args.append("/MDd")
            extra_compile_args.append("/Od")
            extra_compile_args.append("/Ob0")
            extra_compile_args.append("/Zi")
        else:
            if "QL_STATIC_RUNTIME" in os.environ:
                extra_compile_args.append("/MT")
            else:
                extra_compile_args.append("/MD")

    elif compiler == "unix":
        ql_compile_args = os.popen("quantlib-config --cflags").read()[:-1].split()

        extra_compile_args = [
            arg
            for arg in ql_compile_args
            if not arg.startswith("-D")
            if not arg.startswith("-I")
        ] + ["-Wno-unused"]
        if "CXXFLAGS" in os.environ:
            extra_compile_args += os.environ["CXXFLAGS"].split()

    return extra_compile_args


def extra_link_args():

    extra_link_args = []

    compiler = get_default_compiler()

    if compiler == "msvc":

        dbit = round(math.log(sys.maxsize, 2) + 1)
        if dbit == 64:
            machinetype = "/machine:x64"
        else:
            machinetype = "/machine:x86"
        extra_link_args = ["/subsystem:windows", machinetype]
        if is_debug_quantlib():
            extra_link_args += ["/DEBUG"]

    elif compiler == "unix":
        ql_link_args = os.popen("quantlib-config --libs").read()[:-1].split()

        extra_link_args = [
            arg
            for arg in ql_link_args
            if not arg.startswith("-L")
            if not arg.startswith("-l")
        ]
        if "LDFLAGS" in os.environ:
            extra_link_args += os.environ["LDFLAGS"].split()

    return extra_link_args


classifiers = [
    "Development Status :: 6 - Mature",
    "Environment :: Console",
    "Intended Audience :: Developers",
    "Intended Audience :: Education",
    "Intended Audience :: End Users/Desktop",
    "Intended Audience :: Financial and Insurance Industry",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Operating System :: OS Independent",
    "Programming Language :: C++",
    "Programming Language :: Python",
    "Topic :: Office/Business :: Financial",
    "Topic :: Scientific/Engineering",
]

long_description = """
QuantLib (https://www.quantlib.org/) is a free/open-source C++ library
for financial quantitative analysts and developers, aimed at providing
a comprehensive software framework for quantitative finance.
"""


py_limited_api = platform.python_implementation() == "CPython" and not sysconfig.get_config_var("Py_GIL_DISABLED")


with open("./setup.cfg", "w") as f:
    if py_limited_api:	
        f.write('[bdist_wheel]' + os.linesep + 'py_limited_api=cp38' + os.linesep)


setup(
    name="QuantLib",
    version="1.40-rc",
    description="Python bindings for the QuantLib library",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    author="QuantLib Team",
    author_email="quantlib-users@lists.sourceforge.net",
    url="https://www.quantlib.org",
    license="BSD 3-Clause",
    classifiers=classifiers,
    package_dir={"": "src"},
    py_modules=["QuantLib.__init__", "QuantLib.QuantLib"],
    ext_modules=[
        Extension(
            name="QuantLib._QuantLib",
            sources=["src/QuantLib/quantlib_wrap.cpp"],
            py_limited_api=py_limited_api,
            define_macros=define_macros(py_limited_api),
            include_dirs=include_dirs(),
            library_dirs=library_dirs(),
            libraries=libraries(),
            extra_compile_args=extra_compile_args(),
            extra_link_args=extra_link_args(),
        )
    ],
    data_files=[("share/doc/quantlib", ["../LICENSE.TXT"])],
)
