
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005 StatPro Italia srl
 Copyright (C) 2008 Tito Ingargiola
 Copyright (C) 2018, 2019 Matthias Lungwitz

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#ifndef quantlib_montecarlo_tools_i
#define quantlib_montecarlo_tools_i

%include stochasticprocess.i
%include linearalgebra.i
%include randomnumbers.i
%include types.i

%inline %{
Matrix getCovariance(const Array& volatilities, const Matrix& correlations) {
    return QuantLib::getCovariance(volatilities.begin(),
                                   volatilities.end(),
                                   correlations);
}
%}

%{
using QuantLib::Path;
%}

class Path {
    #if defined(SWIGPYTHON)
    %rename(__len__) length;
    #endif
  private:
    Path();
  public:
    Size length() const;
    Real value(Size i) const;
    Real front() const;
    Real back() const;
    Time time(Size i) const;
    %extend {
        #if defined(SWIGPYTHON)
        Real __getitem__(Integer i) {
            Integer size_ = Integer(self->length());
            if (i>=0 && i<size_) {
                return (*self)[i];
            } else if (i<0 && -i<=size_) {
                return (*self)[size_+i];
            } else {
                throw std::out_of_range("path index out of range");
            }
        }
        #endif
    }
};
%template(SamplePath) Sample<Path>;

%{
using QuantLib::PathGenerator;
%}

#if defined(SWIGR)
%rename(nextSample) next;
#endif

template <class GSG>
class PathGenerator {
  public:
    typedef Sample<Path> sample_type;
    PathGenerator(const ext::shared_ptr<StochasticProcess>&,
                      Time length,
                      Size timeSteps,
                      const GSG& generator,
                      bool brownianBridge);
    PathGenerator(const ext::shared_ptr<StochasticProcess>&,
                      const TimeGrid& timeGrid,
                      const GSG& generator,
                      bool brownianBridge);
    const sample_type& next() const;
    const sample_type& antithetic() const;
    Size size() const;
    const TimeGrid& timeGrid() const;
};

%template(GaussianPathGenerator)
    PathGenerator<GaussianRandomSequenceGenerator>;
%template(GaussianSobolPathGenerator)
    PathGenerator<GaussianLowDiscrepancySequenceGenerator>;
%template(InvCumulativeMersenneTwisterPathGenerator)
    PathGenerator<InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>,
                         InverseCumulativeNormal> >;

%{
using QuantLib::MultiPath;
%}

class MultiPath {
    #if defined(SWIGPYTHON)
    %rename(__len__) pathSize;
    #endif
  private:
    MultiPath();
  public:
    Size pathSize() const;
    Size assetNumber() const;
	Path& at(Size j);

    %extend {
        #if defined(SWIGPYTHON)
        const Path& __getitem__(Integer i) {
            Integer assets_ = Integer(self->assetNumber());
            if (i>=0 && i<assets_) {
                return (*self)[i];
            } else if (i<0 && -i<=assets_) {
                return (*self)[assets_+i];
            } else {
                throw std::out_of_range("multi-path index out of range");
            }
        }
        #endif
    }
};

%template(SampleMultiPath) Sample<MultiPath>;

%{
using QuantLib::MultiPathGenerator;
%}

template <class GSG>
class MultiPathGenerator {
  public:
    typedef Sample<MultiPath> sample_type;
    MultiPathGenerator(const ext::shared_ptr<StochasticProcess>&,
                       const TimeGrid& timeGrid,
                       const GSG& generator,
                       bool brownianBridge = false);
    %extend {
      MultiPathGenerator(
                   const ext::shared_ptr<StochasticProcess>& process,
                   const std::vector<Time>& times,
                   const GSG& generator,
                   bool brownianBridge = false) {
          return new MultiPathGenerator<GSG>(process,
                                             TimeGrid(
                                                 times.begin(),
                                                 times.end()),
                                             generator,
                                             brownianBridge);
      }
    }
    const sample_type& next() const;
    const sample_type& antithetic() const;
};

%template(GaussianMultiPathGenerator)
    MultiPathGenerator<GaussianRandomSequenceGenerator>;
%template(GaussianSobolMultiPathGenerator)
    MultiPathGenerator<GaussianLowDiscrepancySequenceGenerator>;

%{
using QuantLib::BrownianBridge;
%}
class BrownianBridge{
public:
  BrownianBridge(Size steps);
  BrownianBridge(const std::vector<Time>& times);
  BrownianBridge(const TimeGrid& timeGrid);

  Size size() const;
  std::vector<Time> times() const;
  std::vector<Real> leftWeight()   const;
  std::vector<Real> rightWeight()  const;
  std::vector<Real> stdDeviation() const;
  %extend{
    std::vector<Real> transform(const std::vector<Real> &input){
      std::vector<Real> outp(input.size());
      $self->transform(input.begin(),input.end(),outp.begin());
      return outp;
    }
    std::vector<unsigned int> bridgeIndex() const{
        return to_vector<unsigned int>($self->bridgeIndex());
    }
    std::vector<unsigned int> leftIndex() const{
        return to_vector<unsigned int>($self->leftIndex());
    }
    std::vector<unsigned int> rightIndex() const{
        return to_vector<unsigned int>($self->rightIndex());
    }
  }
};


#endif
