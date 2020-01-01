/*
 Copyright (C) 2020 Klaus Spanderen

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

#ifndef quantlib_fdm_i
#define quantlib_fdm_i

%include stl.i
%include common.i
%include vectors.i
%include dividends.i
%include stochasticprocess.i

%include boost_shared_ptr.i


// mesher

%{
using QuantLib::Fdm1dMesher;
using QuantLib::FdmBlackScholesMesher;
using QuantLib::Concentrating1dMesher;
using QuantLib::ExponentialJump1dMesher;
using QuantLib::FdmQuantoHelper;
using QuantLib::FdmCEV1dMesher;
using QuantLib::FdmHestonVarianceMesher;
using QuantLib::Uniform1dMesher;
using QuantLib::FdmSimpleProcess1dMesher;
using QuantLib::Predefined1dMesher;
%}


%shared_ptr(Fdm1dMesher)
class Fdm1dMesher {
  public:
    explicit Fdm1dMesher(Size size);
    
    Size size() const;
    Real dplus(Size index) const;
    Real dminus(Size index) const;
    Real location(Size index) const;
    const std::vector<Real>& locations();
};

#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( boost::shared_ptr<Fdm1dMesher> )
#endif
namespace std {
    %template(Fdm1dMesherVector) vector<boost::shared_ptr<Fdm1dMesher> >;
}

%shared_ptr(FdmBlackScholesMesher)
class FdmBlackScholesMesher : public Fdm1dMesher {
  public:
    FdmBlackScholesMesher(
        Size size,
        const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
        Time maturity, Real strike,
        Real xMinConstraint = Null<Real>(),
        Real xMaxConstraint = Null<Real>(),
        Real eps = 0.0001,
        Real scaleFactor = 1.5,
        const std::pair<Real, Real>& cPoint
            = (std::pair<Real, Real>(Null<Real>(), Null<Real>())),
        const std::vector<boost::shared_ptr<Dividend> >& dividendSchedule 
            = std::vector<boost::shared_ptr<Dividend> >(),
        const boost::shared_ptr<FdmQuantoHelper>& fdmQuantoHelper
            = boost::shared_ptr<FdmQuantoHelper>(),
        Real spotAdjustment = 0.0);

    static boost::shared_ptr<GeneralizedBlackScholesProcess> processHelper(
         const Handle<Quote>& s0,
         const Handle<YieldTermStructure>& rTS,
         const Handle<YieldTermStructure>& qTS,
         Volatility vol);
};

#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( boost::tuple<Real, Real, bool> )
#endif

%template(Concentrating1dMesherPoint) boost::tuple<Real, Real, bool>;
%template(Concentrating1dMesherPointVector) std::vector<boost::tuple<Real, Real, bool> >;


%shared_ptr(Concentrating1dMesher)
class Concentrating1dMesher : public Fdm1dMesher {
  public:
    Concentrating1dMesher(
        Real start, Real end, Size size,
        const std::pair<Real, Real>& cPoints
                 = (std::pair<Real, Real>(Null<Real>(), Null<Real>())),
        const bool requireCPoint = false);

    Concentrating1dMesher(
        Real start, Real end, Size size,
        const std::vector<boost::tuple<Real, Real, bool> >& cPoints,
        Real tol = 1e-8);
};

%shared_ptr(ExponentialJump1dMesher)
class ExponentialJump1dMesher : public Fdm1dMesher {
   public:
     ExponentialJump1dMesher(Size steps, Real beta, Real jumpIntensity, 
                             Real eta, Real eps = 1e-3);
};

%shared_ptr(FdmCEV1dMesher)
class FdmCEV1dMesher : public Fdm1dMesher {
  public:
    FdmCEV1dMesher(
        Size size,
        Real f0, Real alpha, Real beta,
        Time maturity,
        Real eps = 0.0001,
        Real scaleFactor = 1.5,
        const std::pair<Real, Real>& cPoint
            = (std::pair<Real, Real>(Null<Real>(), Null<Real>())));
};

%shared_ptr(FdmHestonVarianceMesher)
class FdmHestonVarianceMesher : public Fdm1dMesher {
  public:
    FdmHestonVarianceMesher(
        Size size,
        const boost::shared_ptr<HestonProcess> & process,
        Time maturity, Size tAvgSteps = 10, Real epsilon = 0.0001);

    Real volaEstimate() const;
};

%shared_ptr(FdmSimpleProcess1dMesher)
class FdmSimpleProcess1dMesher : public Fdm1dMesher {
  public:
      FdmSimpleProcess1dMesher(
        Size size,
        const boost::shared_ptr<StochasticProcess1D>& process,
        Time maturity, Size tAvgSteps = 10, Real epsilon = 0.0001,
        Real mandatoryPoint = Null<Real>());
};

%shared_ptr(Uniform1dMesher)
class Uniform1dMesher : public Fdm1dMesher {
  public:
    Uniform1dMesher(Real start, Real end, Size size);
};

%shared_ptr(Predefined1dMesher)
class Predefined1dMesher : public Fdm1dMesher {
  public:
    explicit Predefined1dMesher(const std::vector<Real>& x);
};


%{
using QuantLib::FdmLinearOpIterator;
using QuantLib::FdmLinearOpLayout;
using QuantLib::FdmMesher;
using QuantLib::FdmMesherComposite;
%}


class FdmLinearOpIterator {
  public:
  
#if defined(SWIGPYTHON)
    %extend {
        static FdmLinearOpIterator create(
            const std::vector<unsigned int>& dim,
            const std::vector<unsigned int>& coordinates, Size index) {
            
            std::vector<Size> _dim(dim.size());
            std::copy(dim.begin(), dim.end(), _dim.begin());
            std::vector<Size> _coordinates(coordinates.size());
            std::copy(coordinates.begin(), coordinates.end(), 
                _coordinates.begin());
            
            return FdmLinearOpIterator(_dim, _coordinates, index);
         }
         
         std::vector<unsigned int> coordinates() {
            const std::vector<Size>& c = $self->coordinates();
             std::vector<unsigned int> tmp(c.size());
             std::copy(c.begin(), c.end(), tmp.begin());
             
             return tmp;
         }
     }
#else  
    FdmLinearOpIterator(
        const std::vector<Size>& dim,
        const std::vector<Size>& coordinates, Size index);

    const std::vector<Size>& coordinates();
#endif

    %extend {
        void increment() {
            ++(*$self);
        }
        bool notEqual(const FdmLinearOpIterator& iterator) {
            return $self->operator!=(iterator);
        }
    }    
    Size index() const;
};

%shared_ptr(FdmLinearOpLayout)
class FdmLinearOpLayout {
  public:
  
#if defined(SWIGPYTHON)
    %extend {
        static boost::shared_ptr<FdmLinearOpLayout> create(
            const std::vector<unsigned int>& dim) {
            std::vector<Size> _dim(dim.size());
            std::copy(dim.begin(), dim.end(), _dim.begin());
            
            return boost::make_shared<FdmLinearOpLayout>(_dim);
        }
        
        Size index(const std::vector<unsigned int>& coordinates) const {
            std::vector<Size> tmp(coordinates.size());
            std::copy(coordinates.begin(), coordinates.end(), tmp.begin());
            
            return $self->index(tmp);
        }
        
        const std::vector<unsigned int> spacing() {        
            std::vector<unsigned int> tmp($self->spacing().size());
            std::copy($self->spacing().begin(), $self->spacing().end(),
                tmp.begin());
                
             return tmp;
        }
        const std::vector<unsigned int> dim() const {
            std::vector<unsigned int> tmp($self->dim().size());
            std::copy($self->dim().begin(), $self->dim().end(),
                tmp.begin());
                
             return tmp;
        }        
    }
#else  
    explicit FdmLinearOpLayout(const std::vector<Size>& dim);

    const std::vector<Size>& spacing();
    const std::vector<Size>& dim() const;

    Size index(const std::vector<Size>& coordinates) const;
    
#endif
    
    FdmLinearOpIterator begin() const;
    FdmLinearOpIterator end() const;
    
    Size size() const;
    
    Size neighbourhood(const FdmLinearOpIterator& iterator,
                       Size i, Integer offset) const;

    Size neighbourhood(const FdmLinearOpIterator& iterator,
                       Size i1, Integer offset1,
                       Size i2, Integer offset2) const;

    %extend {
        FdmLinearOpIterator iter_neighbourhood(
            const FdmLinearOpIterator& iterator, Size i, Integer offset) const {
            
            return $self->iter_neighbourhood(iterator, i, offset);
        }
    }
    
#if defined(SWIGPYTHON)
    private:
        FdmLinearOpLayout();    
#endif    
};


%shared_ptr(FdmMesher)
class FdmMesher {
  private:
    FdmMesher();
};

%shared_ptr(FdmMesherComposite)
class FdmMesherComposite : public FdmMesher {
  public:
    FdmMesherComposite(
        const boost::shared_ptr<FdmLinearOpLayout>& layout,
        const std::vector<boost::shared_ptr<Fdm1dMesher> > & mesher);

    // convenient constructors
    explicit FdmMesherComposite(
        const std::vector<boost::shared_ptr<Fdm1dMesher> > & mesher);
    explicit FdmMesherComposite(
        const boost::shared_ptr<Fdm1dMesher>& mesher);
    FdmMesherComposite(const boost::shared_ptr<Fdm1dMesher>& m1,
                       const boost::shared_ptr<Fdm1dMesher>& m2);
    FdmMesherComposite(const boost::shared_ptr<Fdm1dMesher>& m1,
                       const boost::shared_ptr<Fdm1dMesher>& m2,
                       const boost::shared_ptr<Fdm1dMesher>& m3);
    FdmMesherComposite(const boost::shared_ptr<Fdm1dMesher>& m1,
                       const boost::shared_ptr<Fdm1dMesher>& m2,
                       const boost::shared_ptr<Fdm1dMesher>& m3,
                       const boost::shared_ptr<Fdm1dMesher>& m4);


    Real dplus(const FdmLinearOpIterator& iter, Size direction) const;
    Real dminus(const FdmLinearOpIterator& iter, Size direction) const;
    Real location(const FdmLinearOpIterator& iter, Size direction) const;
    %extend {
        Array locations(Size direction) const {
            return $self->locations(direction);
        }
    }

    const std::vector<boost::shared_ptr<Fdm1dMesher> >&
        getFdm1dMeshers() const;
};


// fdm operators

%{
using QuantLib::FdmLinearOp;
using QuantLib::FdmLinearOpComposite;
%}


%shared_ptr(FdmLinearOp)
class FdmLinearOp {
  public:
    virtual ~FdmLinearOp();
    virtual Disposable<Array> apply(const Array& r) const = 0;
};

%shared_ptr(FdmLinearOpComposite)
class FdmLinearOpComposite : public FdmLinearOp {
  public:
    virtual Size size() const = 0;
    virtual void setTime(Time t1, Time t2) = 0;
    virtual Disposable<Array> apply_mixed(const Array& r) const = 0;    
    virtual Disposable<Array> 
        apply_direction(Size direction, const Array& r) const = 0;
    virtual Disposable<Array> 
        solve_splitting(Size direction, const Array& r, Real s) const = 0;
    virtual Disposable<Array> 
        preconditioner(const Array& r, Real s) const = 0;
};


#if defined(SWIGPYTHON)
%{
class FdmLinearOpCompositeProxy : public FdmLinearOpComposite {
  public:
      FdmLinearOpCompositeProxy(PyObject* callback) : callback_(callback) {
        Py_XINCREF(callback_);
    }
    
    FdmLinearOpCompositeProxy& operator=(const FdmLinearOpCompositeProxy& f) {
        if ((this != &f) && (callback_ != f.callback_)) {
            Py_XDECREF(callback_);
            callback_ = f.callback_;
            Py_XINCREF(callback_);
        }
        return *this;
    }
    
    ~FdmLinearOpCompositeProxy() {
        Py_XDECREF(callback_);
    }
        
    Size size() const {
        PyObject* pyResult = PyObject_CallMethod(callback_,"size", NULL);
        
        QL_ENSURE(pyResult != NULL,
                  "failed to call size() on Python object");
                  
        Size result = PyInt_AsLong(pyResult);
        Py_XDECREF(pyResult);
        
        return result;    
    }

    void setTime(Time t1, Time t2) {
        PyObject* pyResult 
            = PyObject_CallMethod(callback_,"setTime","dd", t1, t2);
            
        QL_ENSURE(pyResult != NULL,
                  "failed to call setTime() on Python object");
                                    
        Py_XDECREF(pyResult);
    }

    Disposable<Array> apply(const Array& r) const {
        return apply(r, "apply");        
    }

    Disposable<Array> apply_mixed(const Array& r) const {
        return apply(r, "apply_mixed");        
    }
    
    Disposable<Array> apply_direction(Size direction, const Array& r) const {
        PyObject* pyArray 
            = SWIG_NewPointerObj(SWIG_as_voidptr(&r), SWIGTYPE_p_Array, 0);
            
        PyObject* pyResult 
            = PyObject_CallMethod(callback_, "apply_direction", "kO", 
                (unsigned long)(direction), pyArray);
            
        Py_XDECREF(pyArray); 
            
        return extractArray(pyResult, "apply_direction");        
    }
    
    Disposable<Array> solve_splitting(
        Size direction, const Array& r, Real s) const {

        PyObject* pyArray 
            = SWIG_NewPointerObj(SWIG_as_voidptr(&r), SWIGTYPE_p_Array, 0);
            
        PyObject* pyResult 
            = PyObject_CallMethod(callback_, "solve_splitting", "kOd", 
                (unsigned long)(direction), pyArray, s);
            
        Py_XDECREF(pyArray); 
            
        return extractArray(pyResult, "solve_splitting");        
    }
    
    Disposable<Array> preconditioner(const Array& r, Real s) const {
        PyObject* pyArray 
            = SWIG_NewPointerObj(SWIG_as_voidptr(&r), SWIGTYPE_p_Array, 0);
            
        PyObject* pyResult 
            = PyObject_CallMethod(callback_, "preconditioner", "Od",pyArray, s);
            
        Py_XDECREF(pyArray); 
            
        return extractArray(pyResult, "preconditioner");        
    }

  private:
      Disposable<Array> extractArray(
          PyObject* pyResult, const std::string& methodName) const {
          
        QL_ENSURE(pyResult != NULL,
                  "failed to call " + methodName + " on Python object");

        QL_ENSURE(pyResult != Py_None, methodName + " returned None");
            
        Array* ptr;            
        const int err = SWIG_ConvertPtr(
            pyResult, (void **) &ptr, SWIGTYPE_p_Array,    SWIG_POINTER_EXCEPTION);

        if (err != 0) {
            Py_XDECREF(pyResult);
            QL_FAIL("return type must be of type QuantLib Array in " 
                + methodName);
        }
        
          Array tmp(*ptr);          
        Py_XDECREF(pyResult);
         
          return tmp;
      }
      
    Disposable<Array> apply(
        const Array& r, const std::string& methodName) const {

        PyObject* pyArray 
            = SWIG_NewPointerObj(SWIG_as_voidptr(&r), SWIGTYPE_p_Array, 0);
            
        PyObject* pyResult 
            = PyObject_CallMethod(callback_, methodName.c_str(), "O", pyArray);
            
        Py_XDECREF(pyArray); 
        
        return extractArray(pyResult, methodName);        
    }
        
    PyObject* callback_;    
};
%}

%shared_ptr(FdmLinearOpCompositeProxy)
class FdmLinearOpCompositeProxy {
  public:
    FdmLinearOpCompositeProxy(PyObject* callback);
    
    Size size() const;
    void setTime(Time t1, Time t2);
      
    Disposable<Array> apply(const Array& r) const;
    Disposable<Array> apply_mixed(const Array& r) const;    
    Disposable<Array> apply_direction(Size direction, const Array& r) const;
    Disposable<Array> solve_splitting(Size direction, const Array& r, Real s) const;
    Disposable<Array> preconditioner(const Array& r, Real s) const;
};

#elif defined(SWIGJAVA) || defined(SWIGCSHARP)

%{
class FdmLinearOpCompositeDelegate {
  public:
      virtual ~FdmLinearOpCompositeDelegate() {}
      
    virtual Size size() const {
        QL_FAIL("implementation of FdmLinearOpCompositeDelegate.size is missing");        
    }
    
    virtual void setTime(Time t1, Time t2) {
        QL_FAIL("implementation of FdmLinearOpCompositeDelegate.setTime is missing");    
    }
      
    virtual Array apply(const Array& r) const {
        QL_FAIL("implementation of FdmLinearOpCompositeDelegate.apply is missing");    
    }
    
    virtual Array apply_mixed(const Array& r) const {
        QL_FAIL("implementation of FdmLinearOpCompositeDelegate.apply_mixed is missing");    
    }    
    
    virtual Array apply_direction(Size direction, const Array& r) const {
        QL_FAIL("implementation of FdmLinearOpCompositeDelegate.apply_direction is missing");    
    }
    
    virtual Array solve_splitting(Size direction, const Array& r, Real s) const {
        QL_FAIL("implementation of FdmLinearOpCompositeDelegate.solve_splitting is missing");        
    }    
    
    virtual Array preconditioner(const Array& r, Real dt) const {
        return solve_splitting(0, r, dt);
    }
};

class FdmLinearOpCompositeProxy : public FdmLinearOpComposite {
  public:
      FdmLinearOpCompositeProxy(FdmLinearOpCompositeDelegate* delegate)
      : delegate_(delegate) {}
      
      Size size() const { return delegate_->size(); }
    void setTime(Time t1, Time t2) { delegate_->setTime(t1, t2); }
    
    Disposable<Array> apply(const Array& r) const {
        Array retVal = delegate_->apply(r);
        return retVal;
    }
    Disposable<Array> apply_mixed(const Array& r) const {
        Array retVal = delegate_->apply_mixed(r);
        return retVal;
    }        
    Disposable<Array> apply_direction(Size direction, const Array& r) const {
        Array retVal = delegate_->apply_direction(direction, r);
        return retVal;
    }
    Disposable<Array> solve_splitting(
        Size direction, const Array& r, Real s) const {
        Array retVal = delegate_->solve_splitting(direction, r, s);
        return retVal;
    }
    Disposable<Array> preconditioner(const Array& r, Real s) const {
        Array retVal = delegate_->preconditioner(r, s);
        return retVal;
    }
               
  private:
      FdmLinearOpCompositeDelegate* const delegate_; 
};
%}

%shared_ptr(FdmLinearOpCompositeProxy)
class FdmLinearOpCompositeProxy : public FdmLinearOpComposite {
  public:
      FdmLinearOpCompositeProxy(FdmLinearOpCompositeDelegate* delegate);
      
      Size size() const;
    void setTime(Time t1, Time t2);
    
    Disposable<Array> apply(const Array& r) const;
    Disposable<Array> apply_mixed(const Array& r) const;
    Disposable<Array> apply_direction(Size direction, const Array& r) const;
    Disposable<Array> solve_splitting(
        Size direction, const Array& r, Real s) const;
    Disposable<Array> preconditioner(const Array& r, Real s) const;
               
  private:
      FdmLinearOpCompositeDelegate* const delegate_; 
};


%feature("director") FdmLinearOpCompositeDelegate;

class FdmLinearOpCompositeDelegate {
  public:
      virtual ~FdmLinearOpCompositeDelegate();
      
    virtual Size size() const;
    virtual void setTime(Time t1, Time t2);
      
    virtual Array apply(const Array& r) const;
    virtual Array apply_mixed(const Array& r) const;    
    virtual Array apply_direction(Size direction, const Array& r) const;
    virtual Array solve_splitting(
        Size direction, const Array& r, Real s) const;    
    virtual Array preconditioner(const Array& r, Real dt) const;
};
    
#endif

%{
using QuantLib::FdmDiscountDirichletBoundary;
using QuantLib::FdmDirichletBoundary;
using QuantLib::BoundaryCondition;
using QuantLib::FdmBlackScholesOp;
using QuantLib::Fdm2dBlackScholesOp;
using QuantLib::FdmBatesOp;
using QuantLib::FdmCEVOp;
using QuantLib::FdmG2Op;
using QuantLib::FdmHestonHullWhiteOp;
using QuantLib::FdmHestonOp;
using QuantLib::FdmHullWhiteOp;
using QuantLib::FdmLocalVolFwdOp;
using QuantLib::FdmOrnsteinUhlenbeckOp;
using QuantLib::FdmSabrOp;
using QuantLib::OrnsteinUhlenbeckProcess;

typedef std::vector<boost::shared_ptr<BoundaryCondition<FdmLinearOp> > > FdmBoundaryConditionSet;
%}



%shared_ptr(BoundaryCondition<FdmLinearOp>);

template <class Operator>
class BoundaryCondition {
  public:
      enum Side { None, Upper, Lower };
        
    virtual ~BoundaryCondition();
    virtual void applyBeforeApplying(Operator&) const = 0;
    virtual void applyAfterApplying(Array&) const = 0;
    virtual void applyBeforeSolving(Operator&, Array& rhs) const = 0;
    virtual void applyAfterSolving(Array&) const = 0;
    virtual void setTime(Time t) = 0;
};

typedef std::vector<boost::shared_ptr<BoundaryCondition<FdmLinearOp> > > FdmBoundaryConditionSet;

%template(BoundaryConditionFdmLinearOp) BoundaryCondition<FdmLinearOp>; 

#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( boost::shared_ptr<BoundaryCondition<FdmLinearOp> > )
#endif

%template(FdmBoundaryConditionSet) std::vector<boost::shared_ptr<BoundaryCondition<FdmLinearOp> > >;

%shared_ptr(FdmDirichletBoundary)
class FdmDirichletBoundary : public BoundaryCondition<FdmLinearOp> {
  public:
    typedef BoundaryCondition<FdmLinearOp>::Side Side;

    FdmDirichletBoundary(const boost::shared_ptr<FdmMesher>& mesher,
                         Real valueOnBoundary, Size direction, Side side);

    void applyBeforeApplying(FdmLinearOp&) const;
    void applyBeforeSolving(FdmLinearOp&, Array&) const;
    void applyAfterApplying(Array&) const;
    void applyAfterSolving(Array&) const;
    void setTime(Time);
    
    Real applyAfterApplying(Real x, Real value) const;
};

%shared_ptr(FdmDiscountDirichletBoundary)
class FdmDiscountDirichletBoundary
        : public BoundaryCondition<FdmLinearOp> {
  public:
    typedef BoundaryCondition<FdmLinearOp>::Side Side;

    FdmDiscountDirichletBoundary(
        const boost::shared_ptr<FdmMesher>& mesher,
        const boost::shared_ptr<YieldTermStructure>& rTS,
        Time maturityTime,
        Real valueOnBoundary,
        Size direction, Side side);

    void setTime(Time);
    void applyBeforeApplying(FdmLinearOp&) const;
    void applyBeforeSolving(FdmLinearOp&, Array&) const;
    void applyAfterApplying(Array&) const;
    void applyAfterSolving(Array&) const;
};


%define DeclareOperator(OperatorName, constructor)
%shared_ptr(OperatorName)
class OperatorName : public FdmLinearOpComposite {
  public:
    OperatorName(constructor);
    
    Size size() const;
    void setTime(Time t1, Time t2);

    Disposable<Array> apply(const Array& r) const;
    Disposable<Array> apply_mixed(const Array& r) const;
    Disposable<Array> apply_direction(Size direction, const Array& r) const;
    Disposable<Array>
        solve_splitting(Size direction, const Array& r, Real s) const;
    Disposable<Array> preconditioner(const Array& r, Real s) const;    
};
%enddef

%define COMMA ,
%enddef

DeclareOperator(FdmBatesOp, 
        const boost::shared_ptr<FdmMesher>& mesher COMMA
        const boost::shared_ptr<BatesProcess>& batesProcess COMMA
        const FdmBoundaryConditionSet& bcSet COMMA
        Size integroIntegrationOrder COMMA
        const boost::shared_ptr<FdmQuantoHelper>& quantoHelper
                                    = boost::shared_ptr<FdmQuantoHelper>()                                    
)

DeclareOperator(FdmBlackScholesOp, 
        const boost::shared_ptr<FdmMesher>& mesher COMMA
        const boost::shared_ptr<GeneralizedBlackScholesProcess>& process COMMA
        Real strike COMMA
        bool localVol = false COMMA
        Real illegalLocalVolOverwrite = -Null<Real>() COMMA
        Size direction = 0 COMMA
        const boost::shared_ptr<FdmQuantoHelper>& quantoHelper
            = boost::shared_ptr<FdmQuantoHelper>()
)

DeclareOperator(Fdm2dBlackScholesOp, 
        const boost::shared_ptr<FdmMesher>& mesher COMMA
        const boost::shared_ptr<GeneralizedBlackScholesProcess>& p1 COMMA
        const boost::shared_ptr<GeneralizedBlackScholesProcess>& p2 COMMA
        Real correlation COMMA
        Time maturity COMMA
        bool localVol = false COMMA
        Real illegalLocalVolOverwrite = -Null<Real>()
)        

DeclareOperator(FdmCEVOp, 
    const boost::shared_ptr<FdmMesher>& mesher COMMA
    const boost::shared_ptr<YieldTermStructure>& rTS COMMA
    Real f0 COMMA Real alpha COMMA Real beta COMMA
    Size direction
)

DeclareOperator(FdmG2Op,
    const boost::shared_ptr<FdmMesher>& mesher COMMA
    const boost::shared_ptr<G2>& model COMMA
    Size direction1 COMMA Size direction2
)

DeclareOperator(FdmHestonHullWhiteOp,
    const boost::shared_ptr<FdmMesher>& mesher COMMA
    const boost::shared_ptr<HestonProcess>& hestonProcess COMMA
    const boost::shared_ptr<HullWhiteProcess>& hwProcess COMMA
    Real equityShortRateCorrelation
)

DeclareOperator(FdmHestonOp,
    const boost::shared_ptr<FdmMesher>& mesher COMMA
    const boost::shared_ptr<HestonProcess>& hestonProcess COMMA
    const boost::shared_ptr<FdmQuantoHelper>& quantoHelper
        = boost::shared_ptr<FdmQuantoHelper>() COMMA
    const boost::shared_ptr<LocalVolTermStructure>& leverageFct
        = boost::shared_ptr<LocalVolTermStructure>()
)

DeclareOperator(FdmHullWhiteOp,
    const boost::shared_ptr<FdmMesher>& mesher COMMA
    const boost::shared_ptr<HullWhite>& model COMMA
    Size direction
)

DeclareOperator(FdmLocalVolFwdOp,
    const boost::shared_ptr<FdmMesher>& mesher COMMA
    const boost::shared_ptr<Quote>& spot COMMA
    const boost::shared_ptr<YieldTermStructure>& rTS COMMA
    const boost::shared_ptr<YieldTermStructure>& qTS COMMA
    const boost::shared_ptr<LocalVolTermStructure>& localVol COMMA
    Size direction = 0
)

DeclareOperator(FdmOrnsteinUhlenbeckOp,
    const boost::shared_ptr<FdmMesher>& mesher COMMA
    const boost::shared_ptr<OrnsteinUhlenbeckProcess>& p COMMA
    const boost::shared_ptr<YieldTermStructure>& rTS COMMA
    Size direction = 0
)

DeclareOperator(FdmSabrOp,
    const boost::shared_ptr<FdmMesher>& mesher COMMA
    const boost::shared_ptr<YieldTermStructure>& rTS COMMA
    Real f0 COMMA
    Real alpha COMMA
    Real beta COMMA
    Real nu COMMA
    Real rho
)

%{
using QuantLib::TripleBandLinearOp;
using QuantLib::FirstDerivativeOp;
using QuantLib::SecondDerivativeOp;
using QuantLib::NinePointLinearOp;

%}

%shared_ptr(TripleBandLinearOp)
class TripleBandLinearOp : public FdmLinearOp {
  public:
    TripleBandLinearOp(Size direction,
                       const boost::shared_ptr<FdmMesher>& mesher);

    Disposable<Array> apply(const Array& r) const;
    Disposable<Array> solve_splitting(const Array& r, Real a,
                                      Real b = 1.0) const;

    Disposable<TripleBandLinearOp> mult(const Array& u) const;
    Disposable<TripleBandLinearOp> multR(const Array& u) const;
    Disposable<TripleBandLinearOp> add(const TripleBandLinearOp& m) const;
    Disposable<TripleBandLinearOp> add(const Array& u) const;

    void axpyb(const Array& a, const TripleBandLinearOp& x,
               const TripleBandLinearOp& y, const Array& b);

    void swap(TripleBandLinearOp& m);
};

%shared_ptr(Disposable<TripleBandLinearOp>)
%template(DisposableTripleBandLinearOp) Disposable<TripleBandLinearOp>;

%shared_ptr(FirstDerivativeOp)
class FirstDerivativeOp : public TripleBandLinearOp {
  public:
    FirstDerivativeOp(Size direction,
                      const boost::shared_ptr<FdmMesher>& mesher);
                          
    Disposable<Array> apply(const Array& r) const;
    Disposable<Array> solve_splitting(
        const Array& r, Real a, Real b = 1.0) const;                     
};

%shared_ptr(SecondDerivativeOp)
class SecondDerivativeOp : public TripleBandLinearOp {
  public:
    SecondDerivativeOp(Size direction,
        const boost::shared_ptr<FdmMesher>& mesher);

    Disposable<Array> apply(const Array& r) const;
    Disposable<Array> solve_splitting(
        const Array& r, Real a, Real b = 1.0) const;                     
};

%shared_ptr(NinePointLinearOp)
class NinePointLinearOp : public FdmLinearOp {
  public:
    NinePointLinearOp(Size d0, Size d1,
        const boost::shared_ptr<FdmMesher>& mesher);

    Disposable<Array> apply(const Array& r) const;
};

%shared_ptr(Disposable<NinePointLinearOp>)
%template(DisposableNinePointLinearOp) Disposable<NinePointLinearOp>;


// fdm schemes

%{
using QuantLib::CraigSneydScheme;
using QuantLib::CrankNicolsonScheme;
using QuantLib::ImplicitEulerScheme;
using QuantLib::DouglasScheme;
using QuantLib::ExplicitEulerScheme;
using QuantLib::HundsdorferScheme;
using QuantLib::MethodOfLinesScheme;
using QuantLib::ModifiedCraigSneydScheme;
%}

%shared_ptr(CraigSneydScheme)
class CraigSneydScheme  {
  public:
    CraigSneydScheme(Real theta, Real mu,
        const boost::shared_ptr<FdmLinearOpComposite> & map,
        const FdmBoundaryConditionSet& bcSet = FdmBoundaryConditionSet());

    void step(Array& a, Time t);
    void setStep(Time dt);
};

%shared_ptr(ImplicitEulerScheme)
class ImplicitEulerScheme {
  public:
    enum SolverType { BiCGstab, GMRES };

    ImplicitEulerScheme(
        const boost::shared_ptr<FdmLinearOpComposite>& map,
        const FdmBoundaryConditionSet& bcSet = FdmBoundaryConditionSet(),
        Real relTol = 1e-8,
        SolverType solverType = BiCGstab);

    void step(Array& a, Time t);
    void setStep(Time dt);

    Size numberOfIterations() const;
};

%shared_ptr(CrankNicolsonScheme)
class CrankNicolsonScheme  {
  public:
    CrankNicolsonScheme(
        Real theta,
        const boost::shared_ptr<FdmLinearOpComposite>& map,
        const FdmBoundaryConditionSet& bcSet = FdmBoundaryConditionSet(),
        Real relTol = 1e-8,
        ImplicitEulerScheme::SolverType solverType
            = ImplicitEulerScheme::BiCGstab);

    void step(Array& a, Time t);
    void setStep(Time dt);

    Size numberOfIterations() const;
};

%shared_ptr(DouglasScheme)
class DouglasScheme  {
  public:
    DouglasScheme(Real theta,
        const boost::shared_ptr<FdmLinearOpComposite> & map,
        const FdmBoundaryConditionSet& bcSet = FdmBoundaryConditionSet());

    void step(Array& a, Time t);
    void setStep(Time dt);
};

%shared_ptr(ExplicitEulerScheme)
class ExplicitEulerScheme  {
  public:
    ExplicitEulerScheme(
        const boost::shared_ptr<FdmLinearOpComposite>& map,
        const FdmBoundaryConditionSet& bcSet = FdmBoundaryConditionSet());

    void step(Array& a, Time t);
    void setStep(Time dt);
};

%shared_ptr(HundsdorferScheme)
class HundsdorferScheme  {
  public:
    HundsdorferScheme(Real theta, Real mu,
        const boost::shared_ptr<FdmLinearOpComposite> & map,
        const FdmBoundaryConditionSet& bcSet = FdmBoundaryConditionSet());

    void step(Array& a, Time t);
    void setStep(Time dt);
};

%shared_ptr(MethodOfLinesScheme)
class MethodOfLinesScheme  {
  public:
    MethodOfLinesScheme(
        const Real eps, const Real relInitStepSize,
        const boost::shared_ptr<FdmLinearOpComposite>& map,
        const FdmBoundaryConditionSet& bcSet = FdmBoundaryConditionSet());

    void step(Array& a, Time t);
    void setStep(Time dt);
};

%shared_ptr(ModifiedCraigSneydScheme)
class ModifiedCraigSneydScheme  {
  public:
    ModifiedCraigSneydScheme(Real theta, Real mu,
        const boost::shared_ptr<FdmLinearOpComposite> & map,
        const FdmBoundaryConditionSet& bcSet = FdmBoundaryConditionSet());

    void step(Array& a, Time t);
    void setStep(Time dt);
};



#endif