
/*
 Copyright (C) 2003 Ferdinando Ametrano
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2016 Gouthaman Balaraman
 Copyright (C) 2019 Matthias Lungwitz

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

#ifndef quantlib_random_numbers_i
#define quantlib_random_numbers_i

%include distributions.i

%{
using QuantLib::Sample;

using QuantLib::LecuyerUniformRng;
using QuantLib::KnuthUniformRng;
using QuantLib::MersenneTwisterUniformRng;
using QuantLib::Xoshiro256StarStarUniformRng;

typedef QuantLib::PseudoRandom::urng_type UniformRandomGenerator;

using QuantLib::CLGaussianRng;
using QuantLib::BoxMullerGaussianRng;
using QuantLib::InverseCumulativeRng;

typedef QuantLib::PseudoRandom::rng_type GaussianRandomGenerator;

using QuantLib::RandomSequenceGenerator;

typedef QuantLib::PseudoRandom::ursg_type UniformRandomSequenceGenerator;
using QuantLib::SobolBrownianGenerator;

using QuantLib::HaltonRsg;
using QuantLib::SobolRsg;
using QuantLib::SobolBrownianBridgeRsg;
using QuantLib::Burley2020SobolRsg;
using QuantLib::Burley2020SobolBrownianBridgeRsg;

typedef QuantLib::LowDiscrepancy::ursg_type
    UniformLowDiscrepancySequenceGenerator;

using QuantLib::InverseCumulativeRsg;

typedef QuantLib::PseudoRandom::rsg_type GaussianRandomSequenceGenerator;
typedef QuantLib::LowDiscrepancy::rsg_type
    GaussianLowDiscrepancySequenceGenerator;
%}

template <class T>
class Sample {
  private:
    Sample();
  public:
    %extend {
        const T& value() { return self->value; }
        Real weight() { return self->weight; }
    }
};

%template(SampleNumber) Sample<Real>;
%template(SampleArray) Sample<Array>;
%template(SampleRealVector) Sample<std::vector<Real> >; 

/************* Uniform number generators *************/

#if defined(SWIGR)
%rename(nextSample) next;
#endif

class LecuyerUniformRng {
  public:
    LecuyerUniformRng(BigInteger seed=0);
    Sample<Real> next() const;
};

class KnuthUniformRng {
  public:
    KnuthUniformRng(BigInteger seed=0);
    Sample<Real> next() const;
};

class MersenneTwisterUniformRng {
  public:
    MersenneTwisterUniformRng(BigInteger seed = 0);
    Sample<Real> next() const;
};

class Xoshiro256StarStarUniformRng {
  public:
    Xoshiro256StarStarUniformRng(BigInteger seed = 0);
    Sample<Real> next() const;
};

class UniformRandomGenerator {
  public:
    UniformRandomGenerator(BigInteger seed=0);
    Sample<Real> next() const;

	%extend {
		// improve performance for direct access. faster version
		Real nextValue() const {
			return (*self).next().value;
		}
	}    
};


/************* Gaussian number generators *************/

template<class RNG> class CLGaussianRng {
  public:
    CLGaussianRng(const RNG& rng);
    Sample<Real> next() const;
};

%template(CentralLimitLecuyerGaussianRng) CLGaussianRng<LecuyerUniformRng>;
%template(CentralLimitKnuthGaussianRng)   CLGaussianRng<KnuthUniformRng>;
%template(CentralLimitMersenneTwisterGaussianRng) CLGaussianRng<MersenneTwisterUniformRng>;
%template(CentralLimitXoshiro256StarStarGaussianRng) CLGaussianRng<Xoshiro256StarStarUniformRng>;

template<class RNG> class BoxMullerGaussianRng {
  public:
    BoxMullerGaussianRng(const RNG& rng);
    Sample<Real> next() const;
};

%template(BoxMullerLecuyerGaussianRng) BoxMullerGaussianRng<LecuyerUniformRng>;
%template(BoxMullerKnuthGaussianRng)   BoxMullerGaussianRng<KnuthUniformRng>;
%template(BoxMullerMersenneTwisterGaussianRng) BoxMullerGaussianRng<MersenneTwisterUniformRng>;
%template(BoxMullerXoshiro256StarStarGaussianRng) BoxMullerGaussianRng<Xoshiro256StarStarUniformRng>;

template<class RNG, class F> class InverseCumulativeRng {
  public:
    InverseCumulativeRng(const RNG& rng);
    Sample<Real> next() const;
};

%template(MoroInvCumulativeLecuyerGaussianRng)
    InverseCumulativeRng<LecuyerUniformRng,MoroInverseCumulativeNormal>;
%template(MoroInvCumulativeKnuthGaussianRng)
    InverseCumulativeRng<KnuthUniformRng,MoroInverseCumulativeNormal>;
%template(MoroInvCumulativeMersenneTwisterGaussianRng)
    InverseCumulativeRng<MersenneTwisterUniformRng,MoroInverseCumulativeNormal>;
%template(MoroInvCumulativeXoshiro256StarStarGaussianRng)
    InverseCumulativeRng<Xoshiro256StarStarUniformRng,MoroInverseCumulativeNormal>;

%template(InvCumulativeLecuyerGaussianRng)
    InverseCumulativeRng<LecuyerUniformRng,InverseCumulativeNormal>;
%template(InvCumulativeKnuthGaussianRng)
    InverseCumulativeRng<KnuthUniformRng,InverseCumulativeNormal>;
%template(InvCumulativeMersenneTwisterGaussianRng)
    InverseCumulativeRng<MersenneTwisterUniformRng,InverseCumulativeNormal>;
%template(InvCumulativeXoshiro256StarStarGaussianRng)
    InverseCumulativeRng<Xoshiro256StarStarUniformRng,InverseCumulativeNormal>;

class GaussianRandomGenerator {
  public:
    GaussianRandomGenerator(const UniformRandomGenerator& rng);
    Sample<Real> next() const;

	%extend {
		// improve performance for direct access, faster version
		Real nextValue() const {
			return (*self).next().value;
		}
	}    
};

/************* Uniform sequence generators *************/


class HaltonRsg {
  public:
    HaltonRsg(Size dimensionality, unsigned long seed = 0,
                  bool randomStart = true, bool randomShift = false);
    const Sample<std::vector<Real> >& nextSequence() const;
    const Sample<std::vector<Real> >& lastSequence() const;
    Size dimension() const;
};

class SobolRsg {
  public:
    enum DirectionIntegers {
            Unit, Jaeckel, SobolLevitan, SobolLevitanLemieux,
            JoeKuoD5, JoeKuoD6, JoeKuoD7,
            Kuo, Kuo2, Kuo3 };
    SobolRsg(Size dimensionality, BigInteger seed=0,
             DirectionIntegers directionIntegers = QuantLib::SobolRsg::Jaeckel);
    const Sample<std::vector<Real> >& nextSequence() const;
    const Sample<std::vector<Real> >& lastSequence() const;
    Size dimension() const;
    void skipTo(Size n);
    %extend{
      std::vector<unsigned int> nextInt32Sequence(){
          return to_vector<unsigned int>($self->nextInt32Sequence());
      }
    }
};

class Burley2020SobolRsg {
  public:
    Burley2020SobolRsg(Size dimensionality,
                       BigInteger seed = 42,
                       SobolRsg::DirectionIntegers directionIntegers = QuantLib::SobolRsg::Jaeckel,
                       BigInteger scrambleSeed = 43);
    const Sample<std::vector<Real> >& nextSequence() const;
    const Sample<std::vector<Real> >& lastSequence() const;
    Size dimension() const;
    %extend{
      std::vector<unsigned int> nextInt32Sequence(){
          return to_vector<unsigned int>($self->nextInt32Sequence());
      }
    }
};


class SobolBrownianBridgeRsg {
  public:
    SobolBrownianBridgeRsg(Size factors, Size steps);
    const Sample<std::vector<Real> >&  nextSequence() const;
    const Sample<std::vector<Real> >&  lastSequence() const;
    Size dimension() const;
};

class Burley2020SobolBrownianBridgeRsg {
  public:
    Burley2020SobolBrownianBridgeRsg(Size factors, Size steps);
    const Sample<std::vector<Real> >&  nextSequence() const;
    const Sample<std::vector<Real> >&  lastSequence() const;
    Size dimension() const;
};


template<class RNG> class RandomSequenceGenerator {
  public:
    RandomSequenceGenerator(Size dimensionality,
                            const RNG& rng);
    RandomSequenceGenerator(Size dimensionality,
                            BigNatural seed = 0);
    const Sample<std::vector<Real> >& nextSequence() const;
    Size dimension() const;
};

%template(LecuyerUniformRsg)
    RandomSequenceGenerator<LecuyerUniformRng>;
%template(KnuthUniformRsg)
    RandomSequenceGenerator<KnuthUniformRng>;
%template(MersenneTwisterUniformRsg)
    RandomSequenceGenerator<MersenneTwisterUniformRng>;
%template(Xoshiro256StarStarUniformRsg)
    RandomSequenceGenerator<Xoshiro256StarStarUniformRng>;

class UniformRandomSequenceGenerator {
  public:
    UniformRandomSequenceGenerator(Size dimensionality,
                                   const UniformRandomGenerator& rng);
    const Sample<std::vector<Real> >& nextSequence() const;
    Size dimension() const;
};

class UniformLowDiscrepancySequenceGenerator {
  public:
    UniformLowDiscrepancySequenceGenerator(
        Size dimensionality,
        BigInteger seed=0,
        SobolRsg::DirectionIntegers directionIntegers = QuantLib::SobolRsg::Jaeckel);
    const Sample<std::vector<Real> >& nextSequence() const;
    Size dimension() const;
};

/************* Gaussian sequence generators *************/

template <class U, class I>
class InverseCumulativeRsg {
  public:
    InverseCumulativeRsg(const U& uniformSequenceGenerator);
    InverseCumulativeRsg(const U& uniformSequenceGenerator,
                             const I& inverseCumulative);
    const Sample<std::vector<Real> >& nextSequence() const;
    Size dimension() const;
};


%template(MoroInvCumulativeLecuyerGaussianRsg)
    InverseCumulativeRsg<RandomSequenceGenerator<LecuyerUniformRng>,
                         MoroInverseCumulativeNormal>;
%template(MoroInvCumulativeKnuthGaussianRsg)
    InverseCumulativeRsg<RandomSequenceGenerator<KnuthUniformRng>,
                         MoroInverseCumulativeNormal>;
%template(MoroInvCumulativeMersenneTwisterGaussianRsg)
    InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>,
                         MoroInverseCumulativeNormal>;
%template(MoroInvCumulativeXoshiro256StarStarGaussianRsg)
    InverseCumulativeRsg<RandomSequenceGenerator<Xoshiro256StarStarUniformRng>,
                         MoroInverseCumulativeNormal>;
%template(MoroInvCumulativeHaltonGaussianRsg)
    InverseCumulativeRsg<HaltonRsg,MoroInverseCumulativeNormal>;
%template(MoroInvCumulativeSobolGaussianRsg)
    InverseCumulativeRsg<SobolRsg,MoroInverseCumulativeNormal>;
%template(MoroInvCumulativeBurley2020SobolGaussianRsg)
    InverseCumulativeRsg<Burley2020SobolRsg,MoroInverseCumulativeNormal>;

%template(InvCumulativeLecuyerGaussianRsg)
    InverseCumulativeRsg<RandomSequenceGenerator<LecuyerUniformRng>,
                         InverseCumulativeNormal>;
%template(InvCumulativeKnuthGaussianRsg)
    InverseCumulativeRsg<RandomSequenceGenerator<KnuthUniformRng>,
                         InverseCumulativeNormal>;
%template(InvCumulativeMersenneTwisterGaussianRsg)
    InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>,
                         InverseCumulativeNormal>;
%template(InvCumulativeXoshiro256StarStarGaussianRsg)
    InverseCumulativeRsg<RandomSequenceGenerator<Xoshiro256StarStarUniformRng>,
                         InverseCumulativeNormal>;
%template(InvCumulativeHaltonGaussianRsg)
    InverseCumulativeRsg<HaltonRsg,InverseCumulativeNormal>;
%template(InvCumulativeSobolGaussianRsg)
    InverseCumulativeRsg<SobolRsg,InverseCumulativeNormal>;
%template(InvCumulativeBurley2020SobolGaussianRsg)
    InverseCumulativeRsg<Burley2020SobolRsg,InverseCumulativeNormal>;

class GaussianRandomSequenceGenerator {
  public:
    GaussianRandomSequenceGenerator(
        const UniformRandomSequenceGenerator& uniformSequenceGenerator);
    const Sample<std::vector<Real> >& nextSequence() const;
    Size dimension() const;
};

class GaussianLowDiscrepancySequenceGenerator {
  public:
    GaussianLowDiscrepancySequenceGenerator(
        const UniformLowDiscrepancySequenceGenerator& u);
    const Sample<std::vector<Real> >& nextSequence() const;
    Size dimension() const;
};



/************* LMM-style sequence generators *************/


%{
using QuantLib::BrownianGenerator;
using QuantLib::MTBrownianGenerator;
using QuantLib::SobolBrownianGenerator;
using QuantLib::BrownianGeneratorFactory;
using QuantLib::MTBrownianGeneratorFactory;
using QuantLib::SobolBrownianGeneratorFactory;
%}

%shared_ptr(BrownianGenerator)
class BrownianGenerator {
  public:
    Real nextStep(std::vector<Real>&);
    Real nextPath();

    Size numberOfFactors() const;
    Size numberOfSteps() const;
  private:
    BrownianGenerator();
};

%shared_ptr(BrownianGeneratorFactory)
class BrownianGeneratorFactory {
  public:
    ext::shared_ptr<BrownianGenerator> create(Size factors,
                                              Size steps) const;
  private:
    BrownianGeneratorFactory();
};

%shared_ptr(MTBrownianGenerator)
class MTBrownianGenerator : public BrownianGenerator {
  public:
    MTBrownianGenerator(Size factors,
                        Size steps,
                        unsigned long seed = 0);
};

%shared_ptr(MTBrownianGeneratorFactory)
class MTBrownianGeneratorFactory : public BrownianGeneratorFactory {
  public:
    MTBrownianGeneratorFactory(unsigned long seed = 0);
};

%shared_ptr(SobolBrownianGenerator)
class SobolBrownianGenerator : public BrownianGenerator {
  public:
    enum Ordering { Factors, Steps, Diagonal };
    SobolBrownianGenerator(Size factors,
                           Size steps,
                           Ordering ordering,
                           unsigned long seed = 0,
                           SobolRsg::DirectionIntegers directionIntegers = SobolRsg::Jaeckel);
};

%shared_ptr(SobolBrownianGeneratorFactory)
class SobolBrownianGeneratorFactory : public BrownianGeneratorFactory {
  public:
    SobolBrownianGeneratorFactory(
                           SobolBrownianGenerator::Ordering ordering,
                           unsigned long seed = 0,
                           SobolRsg::DirectionIntegers directionIntegers = SobolRsg::Jaeckel);
};


#endif
