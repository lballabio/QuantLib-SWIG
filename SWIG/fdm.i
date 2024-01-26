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
%include tuple.i
%include functions.i
%include options.i
%include basketoptions.i
%include dividends.i
%include settings.i
%include shortratemodels.i


// mesher

%{
using QuantLib::Fdm1dMesher;
using QuantLib::FdmBlackScholesMesher;
using QuantLib::Concentrating1dMesher;
using QuantLib::ExponentialJump1dMesher;
using QuantLib::FdmQuantoHelper;
using QuantLib::FdmCEV1dMesher;
using QuantLib::FdmHestonVarianceMesher;
using QuantLib::FdmHestonLocalVolatilityVarianceMesher;
using QuantLib::Uniform1dMesher;
using QuantLib::FdmSimpleProcess1dMesher;
using QuantLib::Predefined1dMesher;
using QuantLib::Glued1dMesher;
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
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<Fdm1dMesher> )
#endif
namespace std {
    %template(Fdm1dMesherVector) vector<ext::shared_ptr<Fdm1dMesher> >;
}

%shared_ptr(FdmBlackScholesMesher)
class FdmBlackScholesMesher : public Fdm1dMesher {
  public:
    #if defined(SWIGPYTHON)
    %feature("kwargs") FdmBlackScholesMesher;
    #endif

    FdmBlackScholesMesher(
        Size size,
        const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
        Time maturity, Real strike,
        doubleOrNull xMinConstraint = Null<Real>(),
        doubleOrNull xMaxConstraint = Null<Real>(),
        Real eps = 0.0001,
        Real scaleFactor = 1.5,
        const std::pair<Real, Real>& cPoint
            = (std::pair<Real, Real>(Null<Real>(), Null<Real>())),
        const std::vector<ext::shared_ptr<Dividend> >& dividendSchedule
            = std::vector<ext::shared_ptr<Dividend> >(),
        const ext::shared_ptr<FdmQuantoHelper>& fdmQuantoHelper
            = ext::shared_ptr<FdmQuantoHelper>(),
        Real spotAdjustment = 0.0);

    static ext::shared_ptr<GeneralizedBlackScholesProcess> processHelper(
         const Handle<Quote>& s0,
         const Handle<YieldTermStructure>& rTS,
         const Handle<YieldTermStructure>& qTS,
         Volatility vol);
};

%template(Concentrating1dMesherPoint) ext::tuple<Real, Real, bool>;
%template(Concentrating1dMesherPointVector) std::vector<ext::tuple<Real, Real, bool> >;


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
        const std::vector<ext::tuple<Real, Real, bool> >& cPoints,
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
    #if defined(SWIGPYTHON)
    %feature("kwargs") FdmCEV1dMesher;
    #endif

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
        const ext::shared_ptr<HestonProcess> & process,
        Time maturity, Size tAvgSteps = 10, Real epsilon = 0.0001);

    Real volaEstimate() const;
};

%shared_ptr(FdmHestonLocalVolatilityVarianceMesher)
class FdmHestonLocalVolatilityVarianceMesher : public Fdm1dMesher {
  public:
    FdmHestonLocalVolatilityVarianceMesher(
        Size size,
        const ext::shared_ptr<HestonProcess>& process,
        const ext::shared_ptr<LocalVolTermStructure>& leverageFct,
        Time maturity, Size tAvgSteps = 10, Real epsilon = 0.0001);

    Real volaEstimate() const;
};


%shared_ptr(FdmSimpleProcess1dMesher)
class FdmSimpleProcess1dMesher : public Fdm1dMesher {
  public:
      FdmSimpleProcess1dMesher(
        Size size,
        const ext::shared_ptr<StochasticProcess1D>& process,
        Time maturity, Size tAvgSteps = 10, Real epsilon = 0.0001,
        doubleOrNull mandatoryPoint = Null<Real>());
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

%shared_ptr(Glued1dMesher)
class Glued1dMesher : public Fdm1dMesher {
  public:
    Glued1dMesher(
        const Fdm1dMesher& leftMesher,
        const Fdm1dMesher& rightMesher);
};


%{
using QuantLib::FdmLinearOpIterator;
using QuantLib::FdmLinearOpLayout;
using QuantLib::FdmMesher;
using QuantLib::FdmMesherComposite;
%}


class FdmLinearOpIterator {
  public:
    %extend {
        FdmLinearOpIterator(const std::vector<unsigned int>& dim) {
            return new FdmLinearOpIterator(to_vector<Size>(dim));
        }    
        FdmLinearOpIterator(const std::vector<unsigned int>& dim,
                            const std::vector<unsigned int>& coordinates,
                            Size index) {
            return new FdmLinearOpIterator(to_vector<Size>(dim),
                                           to_vector<Size>(coordinates),
                                           index);
        }
        std::vector<unsigned int> coordinates() {
            return to_vector<unsigned int>($self->coordinates());
        }
        void increment() {
            ++(*$self);
        }
        bool notEqual(const FdmLinearOpIterator& iterator) {
            return self->operator!=(iterator);
        }
    }

    Size index() const;
};

%shared_ptr(FdmLinearOpLayout)
class FdmLinearOpLayout {
  public:
    %extend {
        FdmLinearOpLayout(const std::vector<unsigned int>& dim) {
            return new FdmLinearOpLayout(to_vector<Size>(dim));
        }

        std::vector<unsigned int> spacing() {
            return to_vector<unsigned int>($self->spacing());
        }

        std::vector<unsigned int> dim() const {
            return to_vector<unsigned int>($self->dim());
        }

        Size index(const std::vector<unsigned int>& coordinates) const {
            return $self->index(to_vector<Size>(coordinates));
        }
    }

    FdmLinearOpIterator begin() const;
    FdmLinearOpIterator end() const;

    Size size() const;

    Size neighbourhood(const FdmLinearOpIterator& iterator,
                       Size i, Integer offset) const;

    Size neighbourhood(const FdmLinearOpIterator& iterator,
                       Size i1, Integer offset1,
                       Size i2, Integer offset2) const;

    FdmLinearOpIterator iter_neighbourhood(
        const FdmLinearOpIterator& iterator, Size i, Integer offset) const;
};


%shared_ptr(FdmMesher)
class FdmMesher {
  private:
    FdmMesher();
  public:
    Real dplus(const FdmLinearOpIterator& iter, Size direction) const;
    Real dminus(const FdmLinearOpIterator& iter, Size direction) const;
    Real location(const FdmLinearOpIterator& iter, Size direction) const;
    Array locations(Size direction) const;
    ext::shared_ptr<FdmLinearOpLayout> layout() const;
};

%shared_ptr(FdmMesherComposite)
class FdmMesherComposite : public FdmMesher {
  public:
    FdmMesherComposite(
        const ext::shared_ptr<FdmLinearOpLayout>& layout,
        const std::vector<ext::shared_ptr<Fdm1dMesher> > & mesher);

    // convenient constructors
    explicit FdmMesherComposite(
        const std::vector<ext::shared_ptr<Fdm1dMesher> > & mesher);
    explicit FdmMesherComposite(
        const ext::shared_ptr<Fdm1dMesher>& mesher);
    FdmMesherComposite(const ext::shared_ptr<Fdm1dMesher>& m1,
                       const ext::shared_ptr<Fdm1dMesher>& m2);
    FdmMesherComposite(const ext::shared_ptr<Fdm1dMesher>& m1,
                       const ext::shared_ptr<Fdm1dMesher>& m2,
                       const ext::shared_ptr<Fdm1dMesher>& m3);
    FdmMesherComposite(const ext::shared_ptr<Fdm1dMesher>& m1,
                       const ext::shared_ptr<Fdm1dMesher>& m2,
                       const ext::shared_ptr<Fdm1dMesher>& m3,
                       const ext::shared_ptr<Fdm1dMesher>& m4);

    const std::vector<ext::shared_ptr<Fdm1dMesher> >&
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
    virtual Array apply(const Array& r) const;

  private:
    FdmLinearOp();
};

%{
class SparseMatrix {
  public:
    std::vector<unsigned int> row_idx, col_idx;
    std::vector<Real> data;    
};
%}

%shared_ptr(SparseMatrix)
class SparseMatrix {
  public:
    std::vector<unsigned int> row_idx, col_idx;
    std::vector<Real> data;    
};

%shared_ptr(FdmLinearOpComposite)
class FdmLinearOpComposite : public FdmLinearOp {
  public:    
    virtual Size size() const;
    virtual void setTime(Time t1, Time t2);

    virtual Array apply_mixed(const Array& r) const;
    virtual Array apply_direction(Size direction, const Array& r) const;
    virtual Array solve_splitting(Size direction, const Array& r, Real s) const;
    virtual Array preconditioner(const Array& r, Real s) const;

    %extend {
        ext::shared_ptr<SparseMatrix> to_sparse_matrix() const {
            
            ext::shared_ptr<SparseMatrix> a = ext::make_shared<SparseMatrix>();
            
            const QuantLib::SparseMatrix m = self->toMatrix(); 
            
            Size entries(0);
            for (auto iter1 = m.begin1(); iter1 != m.end1(); ++iter1)
                entries+=std::distance(iter1.begin(), iter1.end());
        
            a->row_idx.reserve(entries);
            a->col_idx.reserve(entries);
            a->data.reserve(entries);
            
            for (auto iter1 = m.begin1(); iter1 != m.end1(); ++iter1)
                for (auto iter2 = iter1.begin(); iter2 != iter1.end(); ++iter2) {
                    a->row_idx.push_back(iter1.index1());
                    a->col_idx.push_back(iter2.index2());
                    a->data.push_back(*iter2);
                }
    
            return a;
        }        
    }
  private:
      FdmLinearOpComposite();
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
    
    FdmLinearOpCompositeProxy(const FdmLinearOpCompositeProxy& p) 
    : callback_(p.callback_) {
        Py_XINCREF(callback_);
    }

    ~FdmLinearOpCompositeProxy() {
        Py_XDECREF(callback_);
    }

    Size size() const {
        PyObject* pyResult = PyObject_CallMethod(callback_,"size", NULL);

        QL_ENSURE(pyResult != NULL,
                  "failed to call size() on Python object");
        QL_ENSURE(PyLong_Check(pyResult), "size() is not an int");

        Size result = PyLong_AsLong(pyResult);
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

    Array apply(const Array& r) const {
        return apply(r, "apply");        
    }

    Array apply_mixed(const Array& r) const {
        return apply(r, "apply_mixed");        
    }

    Array apply_direction(Size direction, const Array& r) const {
        PyObject* pyArray = SWIG_NewPointerObj(
            SWIG_as_voidptr(&r), SWIGTYPE_p_Array, 0);
            
        PyObject* pyResult 
            = PyObject_CallMethod(callback_, "apply_direction", "kO", 
                (unsigned long)(direction), pyArray);
            
        Py_XDECREF(pyArray); 
            
        return extractArray(pyResult, "apply_direction");        
    }

    Array solve_splitting(Size direction, const Array& r, Real s) const {
        PyObject* pyArray = SWIG_NewPointerObj(
            SWIG_as_voidptr(&r), SWIGTYPE_p_Array, 0);
            
        PyObject* pyResult 
            = PyObject_CallMethod(callback_, "solve_splitting", "kOd", 
                (unsigned long)(direction), pyArray, s);
            
        Py_XDECREF(pyArray); 
            
        return extractArray(pyResult, "solve_splitting");        
    }

    Array preconditioner(const Array& r, Real s) const {
        PyObject* pyArray = SWIG_NewPointerObj(
            SWIG_as_voidptr(&r), SWIGTYPE_p_Array, 0);
            
        PyObject* pyResult 
            = PyObject_CallMethod(callback_, "preconditioner", "Od",pyArray, s);
            
        Py_XDECREF(pyArray); 
            
        return extractArray(pyResult, "preconditioner");        
    }

  private:
    Array apply(const Array& r, const char* methodName) const {
        PyObject* pyArray = SWIG_NewPointerObj(
            SWIG_as_voidptr(&r), SWIGTYPE_p_Array, 0);

        PyObject* pyResult 
            = PyObject_CallMethod(callback_, methodName, "O", pyArray);

        Py_DECREF(pyArray);
        return extractArray(pyResult, methodName);        
    }

  private:        
    PyObject* callback_;    
};
%}

%shared_ptr(FdmLinearOpCompositeProxy)
class FdmLinearOpCompositeProxy : public FdmLinearOpComposite {
  public:
    FdmLinearOpCompositeProxy(PyObject* callback);

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

    Array apply(const Array& r) const {
        Array retVal = delegate_->apply(r);
        return retVal;
    }
    Array apply_mixed(const Array& r) const {
        Array retVal = delegate_->apply_mixed(r);
        return retVal;
    }        
    Array apply_direction(Size direction, const Array& r) const {
        Array retVal = delegate_->apply_direction(direction, r);
        return retVal;
    }
    Array solve_splitting(Size direction, const Array& r, Real s) const {
        Array retVal = delegate_->solve_splitting(direction, r, s);
        return retVal;
    }
    Array preconditioner(const Array& r, Real s) const {
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
using QuantLib::FdmTimeDepDirichletBoundary;
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
using QuantLib::FdmZabrOp;
using QuantLib::FdmDupire1dOp;
using QuantLib::FdmBlackScholesFwdOp;
using QuantLib::FdmHestonFwdOp;
using QuantLib::FdmSquareRootFwdOp;

typedef BoundaryCondition<FdmLinearOp> FdmBoundaryCondition;
typedef std::vector<ext::shared_ptr<FdmBoundaryCondition> > FdmBoundaryConditionSet;
%}

%shared_ptr(FdmBoundaryCondition);
class FdmBoundaryCondition {
   %rename(NoSide) None;

  public:    
    enum Side { None, Upper, Lower }; 

    virtual void applyBeforeApplying(FdmLinearOp&) const;
    virtual void applyAfterApplying(Array&) const;
    virtual void applyBeforeSolving(FdmLinearOp&, Array& rhs) const;
    virtual void applyAfterSolving(Array&) const;
    virtual void setTime(Time t);

  private:
    FdmBoundaryCondition();
};


typedef std::vector<ext::shared_ptr<FdmBoundaryCondition> > FdmBoundaryConditionSet;

#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<FdmBoundaryCondition> )
#endif

%template(FdmBoundaryConditionSet) std::vector<ext::shared_ptr<FdmBoundaryCondition> >;

%shared_ptr(FdmDirichletBoundary)
class FdmDirichletBoundary : public FdmBoundaryCondition {
  public:
    typedef FdmBoundaryCondition::Side Side;

    FdmDirichletBoundary(const ext::shared_ptr<FdmMesher>& mesher,
                         Real valueOnBoundary, Size direction, Side side);

	void applyAfterApplying(Array&) const;
    Real applyAfterApplying(Real x, Real value) const;
};

%shared_ptr(FdmDiscountDirichletBoundary)
class FdmDiscountDirichletBoundary : public FdmBoundaryCondition {
  public:
    typedef FdmBoundaryCondition::Side Side;

    FdmDiscountDirichletBoundary(
        const ext::shared_ptr<FdmMesher>& mesher,
        const ext::shared_ptr<YieldTermStructure>& rTS,
        Time maturityTime,
        Real valueOnBoundary,
        Size direction, Side side);
};

#if defined(SWIGPYTHON) || defined(SWIGJAVA) || defined(SWIGCSHARP)
%shared_ptr(FdmTimeDepDirichletBoundary)
class FdmTimeDepDirichletBoundary : public FdmBoundaryCondition {
  public:
    typedef FdmBoundaryCondition::Side Side;

    %extend {
#if defined(SWIGPYTHON)
        FdmTimeDepDirichletBoundary(
            const ext::shared_ptr<FdmMesher>& mesher,
            PyObject* function,
            Size direction, Side side) {

            const ext::function<Real(Real)> f = UnaryFunction(function);
            return new FdmTimeDepDirichletBoundary(
                mesher, f, direction, side);
        }
#elif defined(SWIGJAVA) || defined(SWIGCSHARP)
        FdmTimeDepDirichletBoundary(
            const ext::shared_ptr<FdmMesher>& mesher,
            UnaryFunctionDelegate* function,
            Size direction, Side side) {

            const ext::function<Real(Real)> f = UnaryFunction(function);
            return new FdmTimeDepDirichletBoundary(
                mesher, f, direction, side);        
         }
#endif
    }
};
#endif


%shared_ptr(FdmBatesOp)
class FdmBatesOp : public FdmLinearOpComposite {
  public:
    FdmBatesOp(
        const ext::shared_ptr<FdmMesher>& mesher,
        const ext::shared_ptr<BatesProcess>& batesProcess,
        const FdmBoundaryConditionSet& bcSet,
        Size integroIntegrationOrder,
        const ext::shared_ptr<FdmQuantoHelper>& quantoHelper
                                    = ext::shared_ptr<FdmQuantoHelper>());
};

%shared_ptr(FdmBlackScholesOp)
class FdmBlackScholesOp : public FdmLinearOpComposite {
  public:
    FdmBlackScholesOp(
        const ext::shared_ptr<FdmMesher>& mesher,
        const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
        Real strike,
        bool localVol = false,
        doubleOrNull illegalLocalVolOverwrite = -Null<Real>(),
        Size direction = 0,
        const ext::shared_ptr<FdmQuantoHelper>& quantoHelper
            = ext::shared_ptr<FdmQuantoHelper>());
};

%shared_ptr(Fdm2dBlackScholesOp)
class Fdm2dBlackScholesOp : public FdmLinearOpComposite {
  public:  
    Fdm2dBlackScholesOp( 
        const ext::shared_ptr<FdmMesher>& mesher,
        const ext::shared_ptr<GeneralizedBlackScholesProcess>& p1,
        const ext::shared_ptr<GeneralizedBlackScholesProcess>& p2,
        Real correlation,
        Time maturity,
        bool localVol = false,
        doubleOrNull illegalLocalVolOverwrite = -Null<Real>());
};        

%shared_ptr(FdmCEVOp)
class FdmCEVOp : public FdmLinearOpComposite {
  public:
      FdmCEVOp(
        const ext::shared_ptr<FdmMesher>& mesher,
        const ext::shared_ptr<YieldTermStructure>& rTS,
        Real f0, Real alpha, Real beta,
        Size direction);
};

%shared_ptr(FdmG2Op)
class FdmG2Op : public FdmLinearOpComposite {
  public:
    FdmG2Op(
        const ext::shared_ptr<FdmMesher>& mesher,
        const ext::shared_ptr<G2>& model,
        Size direction1, Size direction2);
};

%shared_ptr(FdmHestonHullWhiteOp)
class FdmHestonHullWhiteOp : public FdmLinearOpComposite {
  public:
    FdmHestonHullWhiteOp(
        const ext::shared_ptr<FdmMesher>& mesher,
        const ext::shared_ptr<HestonProcess>& hestonProcess,
        const ext::shared_ptr<HullWhiteProcess>& hwProcess,
        Real equityShortRateCorrelation);
};

%shared_ptr(FdmHestonOp)
class FdmHestonOp : public FdmLinearOpComposite {
  public:
    FdmHestonOp(
	    const ext::shared_ptr<FdmMesher>& mesher,
	    const ext::shared_ptr<HestonProcess>& hestonProcess,
	    const ext::shared_ptr<FdmQuantoHelper>& quantoHelper
	        = ext::shared_ptr<FdmQuantoHelper>(),
	    const ext::shared_ptr<LocalVolTermStructure>& leverageFct
	        = ext::shared_ptr<LocalVolTermStructure>());
};

%shared_ptr(FdmHullWhiteOp)
class FdmHullWhiteOp : public FdmLinearOpComposite {
  public:
    FdmHullWhiteOp(
        const ext::shared_ptr<FdmMesher>& mesher,
        const ext::shared_ptr<HullWhite>& model,
        Size direction);
};

%shared_ptr(FdmLocalVolFwdOp)
class FdmLocalVolFwdOp : public FdmLinearOpComposite {
  public:
      FdmLocalVolFwdOp(
        const ext::shared_ptr<FdmMesher>& mesher,
        const ext::shared_ptr<Quote>& spot,
        const ext::shared_ptr<YieldTermStructure>& rTS,
        const ext::shared_ptr<YieldTermStructure>& qTS,
        const ext::shared_ptr<LocalVolTermStructure>& localVol,
        Size direction = 0);
};

%shared_ptr(FdmOrnsteinUhlenbeckOp)
class FdmOrnsteinUhlenbeckOp : public FdmLinearOpComposite {
  public:
    FdmOrnsteinUhlenbeckOp(
        const ext::shared_ptr<FdmMesher>& mesher,
        const ext::shared_ptr<OrnsteinUhlenbeckProcess>& p,
        const ext::shared_ptr<YieldTermStructure>& rTS,
        Size direction = 0);
};

%shared_ptr(FdmSabrOp)
class FdmSabrOp : public FdmLinearOpComposite {
  public:
      FdmSabrOp(
        const ext::shared_ptr<FdmMesher>& mesher,
        const ext::shared_ptr<YieldTermStructure>& rTS,
        Real f0,
        Real alpha,
        Real beta,
        Real nu,
        Real rho);
};

%shared_ptr(FdmZabrOp)
class FdmZabrOp : public FdmLinearOpComposite {
  public:
    FdmZabrOp(
        const ext::shared_ptr<FdmMesher> & mesher,
        const Real beta,
        const Real nu,
        const Real rho, 
        const Real gamma);
};

%shared_ptr(FdmDupire1dOp)
class FdmDupire1dOp : public FdmLinearOpComposite {
  public:
    FdmDupire1dOp(
        const ext::shared_ptr<FdmMesher> & mesher,
        const Array &localVolatility);
};

%shared_ptr(FdmBlackScholesFwdOp)
class FdmBlackScholesFwdOp : public FdmLinearOpComposite {
  public:
    FdmBlackScholesFwdOp(
        const ext::shared_ptr<FdmMesher>& mesher,
        const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
        Real strike,
        bool localVol = false,
        Real illegalLocalVolOverwrite = -Null<Real>(),
        Size direction = 0);
};

%shared_ptr(FdmSquareRootFwdOp)
class FdmSquareRootFwdOp : public FdmLinearOpComposite {
  public:
    enum TransformationType { Plain, Power, Log };

    FdmSquareRootFwdOp(
        const ext::shared_ptr<FdmMesher>& mesher,
        Real kappa, Real theta, Real sigma,
        Size direction,
        TransformationType type = Plain);
};

%shared_ptr(FdmHestonFwdOp)
class FdmHestonFwdOp : public FdmLinearOpComposite {
  public:
    FdmHestonFwdOp(
        const ext::shared_ptr<FdmMesher>& mesher,
        const ext::shared_ptr<HestonProcess>& process,
        FdmSquareRootFwdOp::TransformationType type 
            = FdmSquareRootFwdOp::Plain,
        const ext::shared_ptr<LocalVolTermStructure> & leverageFct
            = ext::shared_ptr<LocalVolTermStructure>());
};

%{
using QuantLib::TripleBandLinearOp;
using QuantLib::FirstDerivativeOp;
using QuantLib::SecondDerivativeOp;
using QuantLib::NinePointLinearOp;
using QuantLib::SecondOrderMixedDerivativeOp;
using QuantLib::NthOrderDerivativeOp;
%}

%shared_ptr(TripleBandLinearOp)
class TripleBandLinearOp : public FdmLinearOp {
  public:
    TripleBandLinearOp(Size direction,
                       const ext::shared_ptr<FdmMesher>& mesher);

    Array apply(const Array& r) const;
    Array solve_splitting(const Array& r, Real a, Real b = 1.0) const;

    TripleBandLinearOp mult(const Array& u) const;
    TripleBandLinearOp multR(const Array& u) const;
    TripleBandLinearOp add(const TripleBandLinearOp& m) const;
    TripleBandLinearOp add(const Array& u) const;

    void axpyb(const Array& a, const TripleBandLinearOp& x,
               const TripleBandLinearOp& y, const Array& b);
    void swap(TripleBandLinearOp& m);
};


%shared_ptr(FirstDerivativeOp)
class FirstDerivativeOp : public TripleBandLinearOp {
  public:
    FirstDerivativeOp(Size direction,
                      const ext::shared_ptr<FdmMesher>& mesher);
};

%shared_ptr(SecondDerivativeOp)
class SecondDerivativeOp : public TripleBandLinearOp {
  public:
    SecondDerivativeOp(Size direction,
        const ext::shared_ptr<FdmMesher>& mesher);
};

%shared_ptr(NinePointLinearOp)
class NinePointLinearOp : public FdmLinearOp {
  public:
    NinePointLinearOp(Size d0, Size d1,
        const ext::shared_ptr<FdmMesher>& mesher);
};

%shared_ptr(SecondOrderMixedDerivativeOp)
class SecondOrderMixedDerivativeOp : public NinePointLinearOp {
public:
    SecondOrderMixedDerivativeOp(
        Size d0, Size d1, 
        const ext::shared_ptr<FdmMesher>& mesher);
};

%shared_ptr(NthOrderDerivativeOp)
class NthOrderDerivativeOp : public FdmLinearOp {
  public:
    NthOrderDerivativeOp(
        Size direction, Size order, Integer nPoints,
        const ext::shared_ptr<FdmMesher>& mesher);
};


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
        const ext::shared_ptr<FdmLinearOpComposite> & map,
        const FdmBoundaryConditionSet& bcSet = FdmBoundaryConditionSet());

    void step(Array& a, Time t);
    void setStep(Time dt);
};

%shared_ptr(ImplicitEulerScheme)
class ImplicitEulerScheme {
  public:
    enum SolverType { BiCGstab, GMRES };

    #if defined(SWIGPYTHON)
        %feature("kwargs") ImplicitEulerScheme;
    #endif

    ImplicitEulerScheme(
        const ext::shared_ptr<FdmLinearOpComposite>& map,
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
    #if defined(SWIGPYTHON)
        %feature("kwargs") CrankNicolsonScheme;
    #endif

    CrankNicolsonScheme(
        Real theta,
        const ext::shared_ptr<FdmLinearOpComposite>& map,
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
        const ext::shared_ptr<FdmLinearOpComposite> & map,
        const FdmBoundaryConditionSet& bcSet = FdmBoundaryConditionSet());

    void step(Array& a, Time t);
    void setStep(Time dt);
};

%shared_ptr(ExplicitEulerScheme)
class ExplicitEulerScheme  {
  public:
    ExplicitEulerScheme(
        const ext::shared_ptr<FdmLinearOpComposite>& map,
        const FdmBoundaryConditionSet& bcSet = FdmBoundaryConditionSet());

    void step(Array& a, Time t);
    void setStep(Time dt);
};

%shared_ptr(HundsdorferScheme)
class HundsdorferScheme  {
  public:
    HundsdorferScheme(Real theta, Real mu,
        const ext::shared_ptr<FdmLinearOpComposite> & map,
        const FdmBoundaryConditionSet& bcSet = FdmBoundaryConditionSet());

    void step(Array& a, Time t);
    void setStep(Time dt);
};

%shared_ptr(MethodOfLinesScheme)
class MethodOfLinesScheme  {
  public:
    MethodOfLinesScheme(
        const Real eps, const Real relInitStepSize,
        const ext::shared_ptr<FdmLinearOpComposite>& map,
        const FdmBoundaryConditionSet& bcSet = FdmBoundaryConditionSet());

    void step(Array& a, Time t);
    void setStep(Time dt);
};

%shared_ptr(ModifiedCraigSneydScheme)
class ModifiedCraigSneydScheme  {
  public:
    ModifiedCraigSneydScheme(Real theta, Real mu,
        const ext::shared_ptr<FdmLinearOpComposite> & map,
        const FdmBoundaryConditionSet& bcSet = FdmBoundaryConditionSet());

    void step(Array& a, Time t);
    void setStep(Time dt);
};


// step condition

%{
using QuantLib::StepCondition;
using QuantLib::FdmSnapshotCondition;
using QuantLib::FdmAmericanStepCondition;
using QuantLib::FdmArithmeticAverageCondition;
using QuantLib::FdmSimpleSwingCondition;
using QuantLib::FdmBermudanStepCondition;
using QuantLib::FdmSimpleStorageCondition;
using QuantLib::FdmSimpleSwingCondition;
using QuantLib::FdmDividendHandler;
using QuantLib::FdmInnerValueCalculator;
using QuantLib::FdmCellAveragingInnerValue;
using QuantLib::FdmLogInnerValue;
using QuantLib::FdmLogBasketInnerValue;
using QuantLib::FdmZeroInnerValue;
using QuantLib::FdmAffineModelSwapInnerValue;
using QuantLib::FdmStepConditionComposite;
%}

%shared_ptr(StepCondition<Array>);

template <class array_type>
class StepCondition {
  public:
    virtual void applyTo(array_type& a, Time t) const;

  private:
    StepCondition();
};

%template(FdmStepCondition) StepCondition<Array>;

#if defined(SWIGPYTHON)
%{
class FdmStepConditionProxy : public StepCondition<Array> {
  public:
    FdmStepConditionProxy(PyObject* callback) : callback_(callback) {
        Py_XINCREF(callback_);
    }
    
    FdmStepConditionProxy(const FdmStepConditionProxy& p) 
    : callback_(p.callback_) {
        Py_XINCREF(callback_);
    }

    FdmStepConditionProxy& operator=(const FdmStepConditionProxy& f) {
        if ((this != &f) && (callback_ != f.callback_)) {
            Py_XDECREF(callback_);
            callback_ = f.callback_;
            Py_XINCREF(callback_);
        }
        return *this;
    }

    ~FdmStepConditionProxy() {
        Py_XDECREF(callback_);
    }

    void applyTo(Array& a, Time t) const {
        PyObject* pyArray = SWIG_NewPointerObj(
            SWIG_as_voidptr(&a), SWIGTYPE_p_Array, 0);
            
        PyObject* pyResult 
            = PyObject_CallMethod(callback_, "applyTo", "Od",pyArray, t);

        Py_XDECREF(pyArray);
    }
    
  private:       
    PyObject* callback_;    
};
%}

%shared_ptr(FdmStepConditionProxy)
class FdmStepConditionProxy : public StepCondition<Array> {
  public:
    FdmStepConditionProxy(PyObject* callback);
};

#elif defined(SWIGJAVA) || defined(SWIGCSHARP)

%{
class FdmStepConditionDelegate {
  public:
    virtual ~FdmStepConditionDelegate() {}

    virtual void applyTo(Array& a, Time t) const {
        QL_FAIL("implementation of FdmStepCondition.applyTo is missing");        
    }
};

class FdmStepConditionProxy : public StepCondition<Array> {
  public:
    FdmStepConditionProxy(FdmStepConditionDelegate* delegate)
    : delegate_(delegate) {}

    void applyTo(Array& a, Time t) const {
        delegate_->applyTo(a, t);
    }
    
  private:  
      FdmStepConditionDelegate* const delegate_; 
};
%}

%shared_ptr(FdmStepConditionProxy)
class FdmStepConditionProxy : public StepCondition<Array> {
  public:
    FdmStepConditionProxy(FdmStepConditionDelegate* delegate);
};


%feature("director") FdmStepConditionDelegate;

class FdmStepConditionDelegate {
  public:      
    virtual ~FdmStepConditionDelegate();
    virtual void applyTo(Array& a, Time t) const;
};

#endif


%shared_ptr(FdmInnerValueCalculator)
class FdmInnerValueCalculator {
  public:
    virtual Real innerValue(const FdmLinearOpIterator& iter, Time t);
    virtual Real avgInnerValue(const FdmLinearOpIterator& iter, Time t);

  private:
    FdmInnerValueCalculator();
};

#if defined(SWIGPYTHON)
%{
class FdmInnerValueCalculatorProxy : public FdmInnerValueCalculator {
  public:
    FdmInnerValueCalculatorProxy(PyObject* callback) : callback_(callback) {
        Py_XINCREF(callback_);
    }
    
    FdmInnerValueCalculatorProxy(const FdmInnerValueCalculatorProxy& p) 
    : callback_(p.callback_) {
        Py_XINCREF(callback_);
    }

    FdmInnerValueCalculatorProxy& operator=(const FdmInnerValueCalculatorProxy& f) {
        if ((this != &f) && (callback_ != f.callback_)) {
            Py_XDECREF(callback_);
            callback_ = f.callback_;
            Py_XINCREF(callback_);
        }
        return *this;
    }

    ~FdmInnerValueCalculatorProxy() {
        Py_XDECREF(callback_);
    }

    Real innerValue(const FdmLinearOpIterator& iter, Time t) {
        return getValue(iter, t, "innerValue");
    }

    Real avgInnerValue(const FdmLinearOpIterator& iter, Time t) {
        return getValue(iter, t, "avgInnerValue");
    }
    
  private: 
      Real getValue(const FdmLinearOpIterator& iter, Time t, const char* methodName) {
        PyObject* pyIter = SWIG_NewPointerObj(
            SWIG_as_voidptr(&iter), SWIGTYPE_p_FdmLinearOpIterator, 0);

        PyObject* pyResult 
            = PyObject_CallMethod(callback_, methodName, "Od", pyIter, t);

        Py_DECREF(pyIter);

        QL_ENSURE(pyResult != NULL, "failed to call innerValue function on Python object");

        const Real result = PyFloat_AsDouble(pyResult);

        Py_XDECREF(pyResult);

        return result;
      }
            
    PyObject* callback_;    
};
%}

%shared_ptr(FdmInnerValueCalculatorProxy)
class FdmInnerValueCalculatorProxy : public FdmInnerValueCalculator {
  public:
    FdmInnerValueCalculatorProxy(PyObject* callback);
};

#elif defined(SWIGJAVA) || defined(SWIGCSHARP)

%{
class FdmInnerValueCalculatorDelegate {
  public:
    virtual ~FdmInnerValueCalculatorDelegate() {}
      
    virtual Real innerValue(const FdmLinearOpIterator& iter, Time t) {
        QL_FAIL("implementation of FdmInnerValueCalculatorDelegate.innerValue is missing");        
    }
    virtual Real avgInnerValue(const FdmLinearOpIterator& iter, Time t) {
        QL_FAIL("implementation of FdmInnerValueCalculatorDelegate.avgInnerValue is missing");            
    }
};

class FdmInnerValueCalculatorProxy : public FdmInnerValueCalculator {
  public:
    FdmInnerValueCalculatorProxy(FdmInnerValueCalculatorDelegate* delegate)
    : delegate_(delegate) {}

    Real innerValue(const FdmLinearOpIterator& iter, Time t) {
        return delegate_->innerValue(iter, t);
    }
    Real avgInnerValue(const FdmLinearOpIterator& iter, Time t) {
        return delegate_->avgInnerValue(iter, t);
    }
    
  private:  
      FdmInnerValueCalculatorDelegate* const delegate_; 
};
%}

%shared_ptr(FdmInnerValueCalculatorProxy)
class FdmInnerValueCalculatorProxy : public FdmInnerValueCalculator {
  public:
    FdmInnerValueCalculatorProxy(FdmInnerValueCalculatorDelegate* delegate);
};


%feature("director") FdmInnerValueCalculatorDelegate;

class FdmInnerValueCalculatorDelegate {
  public:
    virtual ~FdmInnerValueCalculatorDelegate();

    virtual Real innerValue(const FdmLinearOpIterator& iter, Time t);
    virtual Real avgInnerValue(const FdmLinearOpIterator& iter, Time t);
};
#endif


%shared_ptr(FdmCellAveragingInnerValue)
class FdmCellAveragingInnerValue : public FdmInnerValueCalculator {
  public:

#if defined(SWIGPYTHON)
    %extend {
        FdmCellAveragingInnerValue(
            const ext::shared_ptr<Payoff>& payoff,
            const ext::shared_ptr<FdmMesher>& mesher,
            Size direction,
            PyObject* gridMapping) {

                UnaryFunction f(gridMapping);
                return new FdmCellAveragingInnerValue(payoff, mesher, direction, f);
        }
    }
#elif defined(SWIGJAVA) || defined(SWIGCSHARP)
    %extend {
        FdmCellAveragingInnerValue(
            const ext::shared_ptr<Payoff>& payoff,
            const ext::shared_ptr<FdmMesher>& mesher,
            Size direction,        
            UnaryFunctionDelegate* gridMapping) {
            
                UnaryFunction f(gridMapping);
                return new FdmCellAveragingInnerValue(payoff, mesher, direction, f);            
        }
    }
#endif

    %extend {
        FdmCellAveragingInnerValue(
            const ext::shared_ptr<Payoff>& payoff,
            const ext::shared_ptr<FdmMesher>& mesher,
            Size direction) {
                return new FdmCellAveragingInnerValue(payoff, mesher, direction);            
        }
    }
};


%shared_ptr(FdmLogInnerValue)
class FdmLogInnerValue : public FdmCellAveragingInnerValue {
  public:
    FdmLogInnerValue(const ext::shared_ptr<Payoff>& payoff,
                     const ext::shared_ptr<FdmMesher>& mesher,
                     Size direction);
};


%shared_ptr(FdmLogBasketInnerValue)
class FdmLogBasketInnerValue : public FdmInnerValueCalculator {
  public:
    FdmLogBasketInnerValue(const ext::shared_ptr<BasketPayoff>& payoff,
                           const ext::shared_ptr<FdmMesher>& mesher);
};

%shared_ptr(FdmZeroInnerValue)
class FdmZeroInnerValue : public FdmInnerValueCalculator {
  public:
    FdmZeroInnerValue();
};


%shared_ptr(FdmAffineModelSwapInnerValue<G2>)
%shared_ptr(FdmAffineModelSwapInnerValue<HullWhite>)


#if !defined(SWIGR)

#if !defined(SWIGJAVA)
%template(TimeToDateMap) std::map<Time, Date>;
#endif

template <class ModelType>
class FdmAffineModelSwapInnerValue : public FdmInnerValueCalculator {
  public:
#if defined(SWIGJAVA)
    %extend {
        FdmAffineModelSwapInnerValue(
            const ext::shared_ptr<ModelType>& disModel,
            const ext::shared_ptr<ModelType>& fwdModel,
            const ext::shared_ptr<VanillaSwap>& swap,
            const std::vector<Time>& exerciseTimes,
            const std::vector<Date>& exerciseDates,
            const ext::shared_ptr<FdmMesher>& mesher,
            Size direction) {

            QL_REQUIRE(exerciseTimes.size() == exerciseDates.size(),
                "different exercise dates and times length");

            std::map<Time, Date> t2d;
            for (Size i=0; i < exerciseTimes.size(); ++i) 
                t2d[ exerciseTimes[i] ] = exerciseDates[i];

            return new FdmAffineModelSwapInnerValue<ModelType>(
                disModel, fwdModel, swap, t2d, mesher, direction);
        }
    }
#else
    FdmAffineModelSwapInnerValue(
        const ext::shared_ptr<ModelType>& disModel,
        const ext::shared_ptr<ModelType>& fwdModel,
        const ext::shared_ptr<VanillaSwap>& swap,
        const std::map<Time, Date>& exerciseDates,
        const ext::shared_ptr<FdmMesher>& mesher,
        Size direction);
#endif
};

%template(FdmAffineG2ModelSwapInnerValue) FdmAffineModelSwapInnerValue<G2>;
%template(FdmAffineHullWhiteModelSwapInnerValue) FdmAffineModelSwapInnerValue<HullWhite>;
#endif

%shared_ptr(FdmSnapshotCondition)
class FdmSnapshotCondition : public StepCondition<Array> {
public:
    explicit FdmSnapshotCondition(Time t);

    Time getTime() const;       
    const Array& getValues() const;
};

#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<StepCondition<Array> > )
#endif

%template(FdmStepConditionVector) std::vector<ext::shared_ptr<StepCondition<Array> > > ;
 
%shared_ptr(FdmStepConditionComposite)
class FdmStepConditionComposite : public StepCondition<Array> {
public:
    typedef std::vector<ext::shared_ptr<StepCondition<Array> > > Conditions;
    %extend {
        FdmStepConditionComposite(
            const std::vector<Time> & stoppingTimes,
            const std::vector<ext::shared_ptr<StepCondition<Array> > > & conditions) {
            return new FdmStepConditionComposite(
                std::list<std::vector<Time> >(1, stoppingTimes), 
                std::list<ext::shared_ptr<StepCondition<Array> > >(
                    conditions.begin(), conditions.end()));
        }
    }

    const std::vector<Time>& stoppingTimes() const;
    const std::vector<ext::shared_ptr<StepCondition<Array> > > & conditions() const;

    static ext::shared_ptr<FdmStepConditionComposite> joinConditions(
                const ext::shared_ptr<FdmSnapshotCondition>& c1,
                const ext::shared_ptr<FdmStepConditionComposite>& c2);

    static ext::shared_ptr<FdmStepConditionComposite> vanillaComposite(
         const std::vector<ext::shared_ptr<Dividend> >& schedule,
         const ext::shared_ptr<Exercise>& exercise,
         const ext::shared_ptr<FdmMesher>& mesher,
         const ext::shared_ptr<FdmInnerValueCalculator>& calculator,
         const Date& refDate,
         const DayCounter& dayCounter);
};


%shared_ptr(FdmAmericanStepCondition)
class FdmAmericanStepCondition : public StepCondition<Array> {
  public:
    FdmAmericanStepCondition(
        const ext::shared_ptr<FdmMesher> & mesher,
        const ext::shared_ptr<FdmInnerValueCalculator> & calculator);
};

%shared_ptr(FdmArithmeticAverageCondition)
class FdmArithmeticAverageCondition : public StepCondition<Array> {
  public:
    FdmArithmeticAverageCondition(
        const std::vector<Time> & averageTimes,
        Real, Size pastFixings,
        const ext::shared_ptr<FdmMesher> & mesher,
        Size equityDirection);
};

%shared_ptr(FdmBermudanStepCondition)
class FdmBermudanStepCondition : public StepCondition<Array> {
  public:
    FdmBermudanStepCondition(
        const std::vector<Date> & exerciseDates,
        const Date& referenceDate,
        const DayCounter& dayCounter,
        const ext::shared_ptr<FdmMesher> & mesher,
        const ext::shared_ptr<FdmInnerValueCalculator> & calculator);

    const std::vector<Time>& exerciseTimes() const;
};

%shared_ptr(FdmSimpleStorageCondition)
class FdmSimpleStorageCondition : public StepCondition<Array> {
  public:
    FdmSimpleStorageCondition(
        const std::vector<Time> & exerciseTimes,
        const ext::shared_ptr<FdmMesher>& mesher,
        const ext::shared_ptr<FdmInnerValueCalculator>& calculator,
        Real changeRate);
};

%shared_ptr(FdmSimpleSwingCondition)
class FdmSimpleSwingCondition : public StepCondition<Array> {
  public:
      FdmSimpleSwingCondition(
              const std::vector<Time> & exerciseTimes,
              const ext::shared_ptr<FdmMesher>& mesher,
              const ext::shared_ptr<FdmInnerValueCalculator>& calculator,
              Size swingDirection,
              Size minExercises = 0);
};

%shared_ptr(FdmDividendHandler)
class FdmDividendHandler : public StepCondition<Array> {
  public:
    FdmDividendHandler(const std::vector<ext::shared_ptr<Dividend> >& schedule,
                       const ext::shared_ptr<FdmMesher>& mesher,
                       const Date& referenceDate,
                       const DayCounter& dayCounter,
                       Size equityDirection);

    const std::vector<Time>& dividendTimes() const;
    const std::vector<Date>& dividendDates() const;
    const std::vector<Real>& dividends() const;
};


// solver

%{
using QuantLib::FdmSolverDesc;
using QuantLib::Fdm1DimSolver;
using QuantLib::FdmBackwardSolver;
using QuantLib::Fdm2dBlackScholesSolver;
using QuantLib::Fdm2DimSolver;
using QuantLib::Fdm3DimSolver;
using QuantLib::FdmG2Solver;
using QuantLib::FdmHestonHullWhiteSolver;
using QuantLib::FdmHestonSolver;
using QuantLib::FdmHullWhiteSolver;
using QuantLib::FdmNdimSolver;
%}


struct FdmSolverDesc {
  public:
    %extend {
        FdmSolverDesc(
            const ext::shared_ptr<FdmMesher>& mesher,
            const FdmBoundaryConditionSet& bcSet,
            const ext::shared_ptr<FdmStepConditionComposite>& condition,
            const ext::shared_ptr<FdmInnerValueCalculator>& calculator,
            Time maturity,
            Size timeSteps,
            Size dampingSteps) {
            
            const FdmSolverDesc desc = { 
                mesher, bcSet, condition, calculator, 
                maturity, timeSteps, dampingSteps };
            
            return new FdmSolverDesc(desc);            
        }
        
        ext::shared_ptr<FdmMesher> getMesher() const { return self->mesher; }
        const FdmBoundaryConditionSet& getBcSet() const { return self->bcSet; }
        ext::shared_ptr<FdmStepConditionComposite>
            getStepConditions() const { return self->condition; }
        ext::shared_ptr<FdmInnerValueCalculator>
            getCalculator() const { return self->calculator; }
        Time getMaturity() const { return self->maturity; }
        Size getTimeSteps() const { return self->timeSteps; }
        Size getDampingSteps() const { return self->dampingSteps; }        
    }
};

%shared_ptr(Fdm1DimSolver)
class Fdm1DimSolver {
  public:
    Fdm1DimSolver(const FdmSolverDesc& solverDesc,
                  const FdmSchemeDesc& schemeDesc,
                  const ext::shared_ptr<FdmLinearOpComposite>& op);

    Real interpolateAt(Real x) const;
    Real thetaAt(Real x) const;

    Real derivativeX(Real x) const;
    Real derivativeXX(Real x) const;
};

%shared_ptr(FdmBackwardSolver)
class FdmBackwardSolver {
  public:    
    FdmBackwardSolver(
      const ext::shared_ptr<FdmLinearOpComposite>& map,
      const FdmBoundaryConditionSet& bcSet,
      const ext::shared_ptr<FdmStepConditionComposite> condition,
      const FdmSchemeDesc& schemeDesc);

    void rollback(Array& a, Time from, Time to,
                  Size steps, Size dampingSteps);
};


%shared_ptr(Fdm2dBlackScholesSolver)
class Fdm2dBlackScholesSolver {
  public:
    #if defined(SWIGPYTHON)
    %feature("kwargs") Fdm2dBlackScholesSolver;
    #endif
    
    %extend {    
        Fdm2dBlackScholesSolver(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& p1,
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& p2,
            const Real correlation,
            const FdmSolverDesc& solverDesc,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer(),
            bool localVol = false,
            Real illegalLocalVolOverwrite = -Null<Real>()) {
                return new Fdm2dBlackScholesSolver(
                    Handle<GeneralizedBlackScholesProcess>(p1), 
                    Handle<GeneralizedBlackScholesProcess>(p2), 
                    correlation, solverDesc, schemeDesc, 
                    localVol, illegalLocalVolOverwrite); 
        }
    }
    
    Real valueAt(Real x, Real y) const;
    Real thetaAt(Real x, Real y) const;

    Real deltaXat(Real x, Real y) const;
    Real deltaYat(Real x, Real y) const;
    Real gammaXat(Real x, Real y) const;
    Real gammaYat(Real x, Real y) const;
    Real gammaXYat(Real x, Real y) const;
};


%shared_ptr(Fdm2DimSolver)
class Fdm2DimSolver {
  public:
    Fdm2DimSolver(const FdmSolverDesc& solverDesc,
                  const FdmSchemeDesc& schemeDesc,
                  const ext::shared_ptr<FdmLinearOpComposite>& op);

    Real interpolateAt(Real x, Real y) const;
    Real thetaAt(Real x, Real y) const;

    Real derivativeX(Real x, Real y) const;
    Real derivativeY(Real x, Real y) const;
    Real derivativeXX(Real x, Real y) const;
    Real derivativeYY(Real x, Real y) const;
    Real derivativeXY(Real x, Real y) const;
};


%shared_ptr(Fdm3DimSolver)
class Fdm3DimSolver {
  public:
    Fdm3DimSolver(const FdmSolverDesc& solverDesc,
                  const FdmSchemeDesc& schemeDesc,
                  const ext::shared_ptr<FdmLinearOpComposite>& op);

    void performCalculations() const;

    Real interpolateAt(Real x, Real y, Rate z) const;
    Real thetaAt(Real x, Real y, Rate z) const;
};


%shared_ptr(FdmG2Solver)
class FdmG2Solver {
  public:
    %extend {
        FdmG2Solver(
            const ext::shared_ptr<G2>& model,
            const FdmSolverDesc& solverDesc,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer()) {
                return new FdmG2Solver(Handle<G2>(model), solverDesc, schemeDesc);            
        }
    }
    Real valueAt(Real x, Real y) const;
};


%shared_ptr(FdmHestonHullWhiteSolver)
class FdmHestonHullWhiteSolver {
  public:
    %extend {
        FdmHestonHullWhiteSolver(
            const ext::shared_ptr<HestonProcess>& hestonProcess,
            const ext::shared_ptr<HullWhiteProcess>& hwProcess,
            Rate corrEquityShortRate,
            const FdmSolverDesc& solverDesc,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer()) {
                return new FdmHestonHullWhiteSolver(
                    Handle<HestonProcess>(hestonProcess),
                    Handle<HullWhiteProcess>(hwProcess),
                    corrEquityShortRate, solverDesc, schemeDesc);                    
        }
    }

    Real valueAt(Real s, Real v, Rate r) const;
    Real thetaAt(Real s, Real v, Rate r) const;

    Real deltaAt(Real s, Real v, Rate r, Real eps) const;
    Real gammaAt(Real s, Real v, Rate r, Real eps) const;
};


%shared_ptr(FdmHestonSolver)
class FdmHestonSolver {
  public:
    #if defined(SWIGPYTHON)
    %feature("kwargs") FdmHestonSolver;
    #endif

    %extend {
        FdmHestonSolver(
            const ext::shared_ptr<HestonProcess>& process,
            const FdmSolverDesc& solverDesc,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer(),
            const ext::shared_ptr<FdmQuantoHelper>& quantoHelper
                = ext::shared_ptr<FdmQuantoHelper>(),
            const ext::shared_ptr<LocalVolTermStructure>& leverageFct
                = ext::shared_ptr<LocalVolTermStructure>()) {

                return new FdmHestonSolver(
                    Handle<HestonProcess>(process),
                    solverDesc, schemeDesc, 
                    Handle<FdmQuantoHelper>(quantoHelper), 
                    leverageFct);
        }
    }

    Real valueAt(Real s, Real v) const;
    Real thetaAt(Real s, Real v) const;

    Real deltaAt(Real s, Real v) const;
    Real gammaAt(Real s, Real v) const;

    Real meanVarianceDeltaAt(Real s, Real v) const;
    Real meanVarianceGammaAt(Real s, Real v) const;
};


%shared_ptr(FdmHullWhiteSolver)
class FdmHullWhiteSolver {
  public:
    %extend {
        FdmHullWhiteSolver(
            const ext::shared_ptr<HullWhite>& model,
            const FdmSolverDesc& solverDesc,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer()) {
                return new FdmHullWhiteSolver(
                    Handle<HullWhite>(model), solverDesc, schemeDesc);
        }
    }
    Real valueAt(Real r) const;
};


%shared_ptr(FdmNdimSolver<4>);
%shared_ptr(FdmNdimSolver<5>);
%shared_ptr(FdmNdimSolver<6>);

template <Size N>
class FdmNdimSolver {
  public:
    FdmNdimSolver(const FdmSolverDesc& solverDesc,
                  const FdmSchemeDesc& schemeDesc,
                  const ext::shared_ptr<FdmLinearOpComposite>& op);

    Real interpolateAt(const std::vector<Real>& x) const;
    Real thetaAt(const std::vector<Real>& x) const;
};


%template(Fdm4dimSolver) FdmNdimSolver<4>;
%template(Fdm5dimSolver) FdmNdimSolver<5>;
%template(Fdm6dimSolver) FdmNdimSolver<6>;


// utilities

%{
using QuantLib::FdmIndicesOnBoundary;
using QuantLib::RiskNeutralDensityCalculator;
using QuantLib::BSMRNDCalculator;
using QuantLib::CEVRNDCalculator;
using QuantLib::GBSMRNDCalculator;
using QuantLib::HestonRNDCalculator;
using QuantLib::LocalVolRNDCalculator;
using QuantLib::SquareRootProcessRNDCalculator;
%}

%shared_ptr(FdmIndicesOnBoundary)
class FdmIndicesOnBoundary {
  public:
    FdmIndicesOnBoundary(const ext::shared_ptr<FdmLinearOpLayout>& l,
                          Size direction, FdmDirichletBoundary::Side side);

    %extend {
        std::vector<unsigned int> getIndices() const {
            return to_vector<unsigned int>($self->getIndices());
        }
    }
};


%shared_ptr(RiskNeutralDensityCalculator)
class RiskNeutralDensityCalculator {
  public:
    virtual Real pdf(Real x, Time t) const;
    virtual Real cdf(Real x, Time t) const;
    virtual Real invcdf(Real p, Time t) const;

  private:
    RiskNeutralDensityCalculator();
};

%shared_ptr(BSMRNDCalculator)
class BSMRNDCalculator : public RiskNeutralDensityCalculator {
  public:
    explicit BSMRNDCalculator(
        const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
};

%shared_ptr(CEVRNDCalculator)
class CEVRNDCalculator : public RiskNeutralDensityCalculator {
  public:
    CEVRNDCalculator(Real f0, Real alpha, Real beta);

    Real massAtZero(Time t) const;
};

%shared_ptr(GBSMRNDCalculator)
class GBSMRNDCalculator : public RiskNeutralDensityCalculator {
public:
    explicit GBSMRNDCalculator(
        const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
};

%shared_ptr(HestonRNDCalculator)
class HestonRNDCalculator : public RiskNeutralDensityCalculator {
public:
    HestonRNDCalculator(
        const ext::shared_ptr<HestonProcess>& hestonProcess,
        Real integrationEps= 1e-6,
        Size maxIntegrationIterations = 10000ul);
};


%shared_ptr(LocalVolRNDCalculator)
class LocalVolRNDCalculator : public RiskNeutralDensityCalculator {
  public:
#if defined(SWIGPYTHON)
%feature("kwargs") LocalVolRNDCalculator;
#endif

    LocalVolRNDCalculator(
        const ext::shared_ptr<Quote>& spot,
        const ext::shared_ptr<YieldTermStructure>& rTS,
        const ext::shared_ptr<YieldTermStructure>& qTS,
        const ext::shared_ptr<LocalVolTermStructure>& localVol,
        Size xGrid = 101, Size tGrid = 51,
        Real x0Density = 0.1,
        Real localVolProbEps = 1e-6,
        Size maxIter = 10000,
        Time gaussianStepSize = -Null<Time>());

    ext::shared_ptr<Fdm1dMesher> mesher(Time t) const;
    %extend {
        std::vector<unsigned int> rescaleTimeSteps() const {
            return to_vector<unsigned int>($self->rescaleTimeSteps());
        }
    }
};

%shared_ptr(SquareRootProcessRNDCalculator)
class SquareRootProcessRNDCalculator : public RiskNeutralDensityCalculator {
  public:
    SquareRootProcessRNDCalculator(
        Real v0, Real kappa, Real theta, Real sigma);

    Real stationary_pdf(Real v) const;
    Real stationary_cdf(Real v) const;
    Real stationary_invcdf(Real q) const;
};


#endif
