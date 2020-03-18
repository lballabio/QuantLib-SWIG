try:
    from setuptools import setup
except:
    from distutils.core import setup

classifiers = [
    "Development Status :: 5 - Production/Stable",
    "Environment :: Console",
    "Intended Audience :: Developers",
    "Intended Audience :: Science/Research",
    "Intended Audience :: End Users/Desktop",
    "License :: OSI Approved :: BSD License",
    "Natural Language :: English",
    "Programming Language :: C++",
    "Programming Language :: Python",
    "Topic :: Scientific/Engineering",
    "Operating System :: Microsoft :: Windows",
    "Operating System :: POSIX",
    "Operating System :: Unix",
    "Operating System :: MacOS",
]

setup(
    name="QuantLib-Python",
    version="1.18",
    description="Backward-compatible meta-package for the QuantLib module",
    long_description="""
    This module is provided for backward compatibility.
    Use "pip install QuantLib" instead.
    """,
    author="QuantLib Team",
    author_email="quantlib-users@lists.sourceforge.net",
    url="http://quantlib.org",
    license="BSD 3-Clause",
    classifiers=classifiers,
    install_requires=["QuantLib"],
)
