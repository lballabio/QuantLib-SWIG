
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005 StatPro Italia srl
 Copyright (C) 2008 Tito Ingargiola

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

#if defined(SWIGRUBY)
%mixin Path "Enumerable";
#endif
class Path {
    #if defined(SWIGPYTHON) || defined(SWIGRUBY)
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
        #if defined(SWIGPYTHON) || defined(SWIGRUBY)
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
        #if defined(SWIGRUBY)
        void each() {
            for (Size i=0; i<self->length(); i++)
                rb_yield(rb_float_new((*self)[i]));
        }
        #endif
    }
};

%{
typedef QuantLib::PathGenerator<GaussianRandomSequenceGenerator>
    GaussianPathGenerator;
%}
%template(SamplePath) Sample<Path>;
class GaussianPathGenerator {
  public:
    %extend {
        GaussianPathGenerator(const StochasticProcess1DPtr& process,
                              Time length, Size steps,
                              const GaussianRandomSequenceGenerator& rsg,
                              bool brownianBridge) {
            boost::shared_ptr<StochasticProcess1D> process1d =
                boost::dynamic_pointer_cast<StochasticProcess1D>(process);
            return new GaussianPathGenerator(process1d,length,steps,
                                             rsg,brownianBridge);
        }
    }
    Sample<Path> next() const;
    Sample<Path> antithetic() const;
};

%{
typedef QuantLib::PathGenerator<GaussianLowDiscrepancySequenceGenerator>
    GaussianSobolPathGenerator;
%}
class GaussianSobolPathGenerator {
  public:
    %extend {
        GaussianSobolPathGenerator(
                           const StochasticProcess1DPtr& process,
                           Time length, Size steps,
                           const GaussianLowDiscrepancySequenceGenerator& rsg,
                           bool brownianBridge) {
            boost::shared_ptr<StochasticProcess1D> process1d =
                boost::dynamic_pointer_cast<StochasticProcess1D>(process);
            return new GaussianSobolPathGenerator(process1d,length,steps,
                                                  rsg,brownianBridge);
        }
    }
    Sample<Path> next() const;
    Sample<Path> antithetic() const;
};


%{
using QuantLib::MultiPath;
%}

class MultiPath {
    #if defined(SWIGPYTHON) || defined(SWIGRUBY)
    %rename(__len__)        pathSize;
    #endif
  private:
    MultiPath();
  public:
    Size pathSize() const;
    Size assetNumber() const;
	Path& at(Size j);

    %extend {
        #if defined(SWIGPYTHON) || defined(SWIGRUBY)
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
        #if defined(SWIGRUBY)
        void each_path() {
            for (Size i=0; i<self->assetNumber(); i++)
                rb_yield(SWIG_NewPointerObj(&((*self)[i]),
                                            $descriptor(Path *), 0));
        }
        void each_step() {
            for (Size j=0; j<self->pathSize(); j++) {
                VALUE v = rb_ary_new2(self->assetNumber());
                for (Size i=0; i<self->assetNumber(); i++)
                    rb_ary_store(v,i,rb_float_new((*self)[i][j]));
                rb_yield(v);
            }
        }
        #endif
    }
};

%{
typedef QuantLib::MultiPathGenerator<GaussianRandomSequenceGenerator>
    GaussianMultiPathGenerator;
%}
%template(SampleMultiPath) Sample<MultiPath>;
class GaussianMultiPathGenerator {
  public:
    %extend {
      GaussianMultiPathGenerator(
                   const boost::shared_ptr<StochasticProcess>& process,
                   const std::vector<Time>& times,
                   const GaussianRandomSequenceGenerator& generator,
                   bool brownianBridge = false) {
          return new GaussianMultiPathGenerator(process,
                                                QuantLib::TimeGrid(
                                                    times.begin(),
                                                    times.end()),
                                                generator,
                                                brownianBridge);
      }
    }
    Sample<MultiPath> next() const;
	Sample<MultiPath> antithetic() const;
};

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
    	const std::vector<Size> &tmp = $self->bridgeIndex();
    	std::vector<unsigned int> outp(tmp.size());
    	std::copy(tmp.begin(), tmp.end(), outp.begin());
    	return outp;
    }
    std::vector<unsigned int> leftIndex() const{
    	const std::vector<Size> &tmp = $self->leftIndex();
    	std::vector<unsigned int> outp(tmp.size());
    	std::copy(tmp.begin(), tmp.end(), outp.begin());
    	return outp;
    }
    std::vector<unsigned int> rightIndex() const{
    	const std::vector<Size> &tmp = $self->rightIndex();
    	std::vector<unsigned int> outp(tmp.size());
    	std::copy(tmp.begin(), tmp.end(), outp.begin());
    	return outp;
    }
  }
};


#endif
