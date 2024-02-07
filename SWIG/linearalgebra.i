
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005 StatPro Italia srl
 Copyright (C) 2005 Dominic Thuillier

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

#ifndef quantlib_linear_algebra_i
#define quantlib_linear_algebra_i

%include common.i
%include types.i
%include stl.i

%{
using QuantLib::Array;
using QuantLib::Matrix;
using QuantLib::SampledCurve;
%}

%define QL_TYPECHECK_ARRAY       4210    %enddef
%define QL_TYPECHECK_MATRIX      4220    %enddef

#if defined(SWIGPYTHON)
%{
bool ArrayFromSequence(PyObject* source, Array* target) {
    if (PyTuple_Check(source) || PyList_Check(source)) {
        Size size = (PyTuple_Check(source) ?
                     PyTuple_Size(source) :
                     PyList_Size(source));
        *target = Array(size);
        for (Size i=0; i<size; i++) {
            PyObject* o = PySequence_GetItem(source,i);
            if (PyFloat_Check(o)) {
                (*target)[i] = PyFloat_AsDouble(o);
                Py_DECREF(o);
            } else if (PyLong_Check(o)) {
                (*target)[i] = PyLong_AsDouble(o);
                Py_DECREF(o);
            } else {
                Py_DECREF(o);
                return false;
            }
        }
        return true;
    } else {
        return false;
    }
}
%}

%typemap(in) Array (Array* v, void *argp, int res = 0) {
    if (ArrayFromSequence($input,&$1)) {
        ;
    } else {
        // copied from SWIGTYPE typemap -- might need updating for newer SWIG
        res = SWIG_ConvertPtr($input, &argp, $&descriptor, %convertptr_flags);
        if (!SWIG_IsOK(res)) {
            %argument_fail(res, "$type", $symname, $argnum);
        }
        if (!argp) {
            %argument_nullref("$type", $symname, $argnum);
        } else {
            $1 = *(%reinterpret_cast(argp, $&ltype));
        }
    }
};
%typemap(in) const Array& (Array temp, void *argp = 0, int res = 0) {
    if (ArrayFromSequence($input,&temp)) {
        $1 = &temp;
    } else {
        // copied from SWIGTYPE typemap -- might need updating for newer SWIG
        res = SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags);
        if (!SWIG_IsOK(res)) {
            %argument_fail(res, "$type", $symname, $argnum);
        }
        if (!argp) { %argument_nullref("$type", $symname, $argnum); }
        $1 = %reinterpret_cast(argp, $ltype);
    }
};
%typecheck(QL_TYPECHECK_ARRAY) Array {
    /* native sequence? */
    if (PyTuple_Check($input) || PyList_Check($input)) {
        Size size = PySequence_Size($input);
        if (size == 0) {
            $1 = 1;
        } else {
            PyObject* o = PySequence_GetItem($input,0);
            if (PyNumber_Check(o))
                $1 = 1;
            else
                $1 = 0;
            Py_DECREF(o);
        }
    } else {
        // copied from SWIGTYPE typemap -- might need updating for newer SWIG
        void *vptr = 0;
        int res = SWIG_ConvertPtr($input, &vptr, $&descriptor, SWIG_POINTER_NO_NULL);
        $1 = SWIG_CheckState(res);
    }
}
%typecheck(QL_TYPECHECK_ARRAY) const Array & {
    /* native sequence? */
    if (PyTuple_Check($input) || PyList_Check($input)) {
        Size size = PySequence_Size($input);
        if (size == 0) {
            $1 = 1;
        } else {
            PyObject* o = PySequence_GetItem($input,0);
            if (PyNumber_Check(o))
                $1 = 1;
            else
                $1 = 0;
            Py_DECREF(o);
        }
    } else {
        // copied from SWIGTYPE typemap -- might need updating for newer SWIG
        void *vptr = 0;
        int res = SWIG_ConvertPtr($input, &vptr, $descriptor, SWIG_POINTER_NO_NULL);
        $1 = SWIG_CheckState(res);
    }
}



%typemap(in) Matrix (Matrix* m, void *argp, int res = 0) {
    if (PyTuple_Check($input) || PyList_Check($input)) {
        Size rows, cols;
        rows = (PyTuple_Check($input) ?
                PyTuple_Size($input) :
                PyList_Size($input));
        if (rows > 0) {
            // look ahead
            PyObject* o = PySequence_GetItem($input,0);
            if (PyTuple_Check(o) || PyList_Check(o)) {
                cols = (PyTuple_Check(o) ?
                        PyTuple_Size(o) :
                        PyList_Size(o));
                Py_DECREF(o);
            } else {
                PyErr_SetString(PyExc_TypeError, "Matrix expected");
                Py_DECREF(o);
                SWIG_fail;
            }
        } else {
            cols = 0;
        }
        $1 = Matrix(rows,cols);
        for (Size i=0; i<rows; i++) {
            PyObject* o = PySequence_GetItem($input,i);
            if (PyTuple_Check(o) || PyList_Check(o)) {
                Size items = (PyTuple_Check(o) ?
                              PyTuple_Size(o) :
                              PyList_Size(o));
                if (items != cols) {
                    PyErr_SetString(PyExc_TypeError,
                        "Matrix must have equal-length rows");
                    Py_DECREF(o);
                    SWIG_fail;
                }
                for (Size j=0; j<cols; j++) {
                    PyObject* d = PySequence_GetItem(o,j);
                    if (PyFloat_Check(d)) {
                        $1[i][j] = PyFloat_AsDouble(d);
                        Py_DECREF(d);
                    } else if (PyLong_Check(d)) {
                        $1[i][j] = PyLong_AsDouble(d);
                        Py_DECREF(d);
                    } else {
                        PyErr_SetString(PyExc_TypeError,"doubles expected");
                        Py_DECREF(d);
                        Py_DECREF(o);
                        SWIG_fail;
                    }
                }
                Py_DECREF(o);
            } else {
                PyErr_SetString(PyExc_TypeError, "Matrix expected");
                Py_DECREF(o);
                SWIG_fail;
            }
        }
    } else {
        // copied from SWIGTYPE typemap -- might need updating for newer SWIG
        res = SWIG_ConvertPtr($input, &argp, $&descriptor, %convertptr_flags);
        if (!SWIG_IsOK(res)) {
            %argument_fail(res, "$type", $symname, $argnum);
        }
        if (!argp) {
            %argument_nullref("$type", $symname, $argnum);
        } else {
            $1 = *(%reinterpret_cast(argp, $&ltype));
        }
    }
};
%typemap(in) const Matrix & (Matrix temp, void *argp = 0, int res = 0) {
    if (PyTuple_Check($input) || PyList_Check($input)) {
        Size rows, cols;
        rows = (PyTuple_Check($input) ?
                PyTuple_Size($input) :
                PyList_Size($input));
        if (rows > 0) {
            // look ahead
            PyObject* o = PySequence_GetItem($input,0);
            if (PyTuple_Check(o) || PyList_Check(o)) {
                cols = (PyTuple_Check(o) ?
                        PyTuple_Size(o) :
                        PyList_Size(o));
                Py_DECREF(o);
            } else {
                PyErr_SetString(PyExc_TypeError, "Matrix expected");
                Py_DECREF(o);
                SWIG_fail;
            }
        } else {
            cols = 0;
        }

        temp = Matrix(rows,cols);
        for (Size i=0; i<rows; i++) {
            PyObject* o = PySequence_GetItem($input,i);
            if (PyTuple_Check(o) || PyList_Check(o)) {
                Size items = (PyTuple_Check(o) ?
                                        PyTuple_Size(o) :
                                        PyList_Size(o));
                if (items != cols) {
                    PyErr_SetString(PyExc_TypeError,
                        "Matrix must have equal-length rows");
                    Py_DECREF(o);
                    SWIG_fail;
                }
                for (Size j=0; j<cols; j++) {
                    PyObject* d = PySequence_GetItem(o,j);
                    if (PyFloat_Check(d)) {
                        temp[i][j] = PyFloat_AsDouble(d);
                        Py_DECREF(d);
                    } else if (PyLong_Check(d)) {
                        temp[i][j] = PyLong_AsDouble(d);
                        Py_DECREF(d);
                    } else {
                        PyErr_SetString(PyExc_TypeError,"doubles expected");
                        Py_DECREF(d);
                        Py_DECREF(o);
                        SWIG_fail;
                    }
                }
                Py_DECREF(o);
            } else {
                PyErr_SetString(PyExc_TypeError, "Matrix expected");
                Py_DECREF(o);
                SWIG_fail;
            }
        }
        $1 = &temp;
    } else {
        // copied from SWIGTYPE typemap -- might need updating for newer SWIG
        res = SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags);
        if (!SWIG_IsOK(res)) {
            %argument_fail(res, "$type", $symname, $argnum);
        }
        if (!argp) { %argument_nullref("$type", $symname, $argnum); }
        $1 = %reinterpret_cast(argp, $ltype);
    }
};
%typecheck(QL_TYPECHECK_MATRIX) Matrix {
    /* native sequence? */
    if (PyTuple_Check($input) || PyList_Check($input)) {
        $1 = 1;
    /* wrapped Matrix? */
    } else {
        // copied from SWIGTYPE typemap -- might need updating for newer SWIG
        void *vptr = 0;
        int res = SWIG_ConvertPtr($input, &vptr, $&descriptor, SWIG_POINTER_NO_NULL);
        $1 = SWIG_CheckState(res);
    }
}
%typecheck(QL_TYPECHECK_MATRIX) const Matrix & {
    /* native sequence? */
    if (PyTuple_Check($input) || PyList_Check($input)) {
        $1 = 1;
    /* wrapped Matrix? */
    } else {
        // copied from SWIGTYPE typemap -- might need updating for newer SWIG
        void *vptr = 0;
        int res = SWIG_ConvertPtr($input, &vptr, $descriptor, SWIG_POINTER_NO_NULL);
        $1 = SWIG_CheckState(res);
    }
}
#endif

#if defined(SWIGR)
swigr_list_converter(Array,_p_Array,numeric)
%Rruntime %{
setMethod('print', '_p_Matrix',
function(x) print(as.matrix(x)))

setMethod("as.matrix", "_p_Matrix",
function(x) matrix(data=as.numeric(x$dataVector),
        nrow=x$rows(), ncol=x$columns()))

setMethod("print", "_p_SampledCurve",
function(x) print(as.data.frame(x))
)

setMethod("as.data.frame", "_p_SampledCurve",
function(x,row.names,optional)
data.frame("grid"=as(x$grid(), "numeric"),
"values"=as(x$values(), "numeric")))

setMethod("plot", "_p_SampledCurve",
function(x,y) plot(as.data.frame(x)))

%}
#endif

#if defined(SWIGR)
%Rruntime %{
setMethod("+", c("_p_Array", "_p_Array"),
    function(e1,e2) Array___add__(e1,e2))
setMethod("-", c("_p_Array", "_p_Array"),
    function(e1,e2) Array___sub__(e1,e2))
setMethod("*", c("_p_Array", "_p_Array"),
    function(e1,e2) Array___mul__(e1,e2))
setMethod("*", c("_p_Array", "numeric"),
    function(e1,e2) Array___mul__(e1,e2))
setMethod("/", c("_p_Array", "numeric"),
    function(e1,e2) Array___div__(e1,e2))    
%}
#endif


#if defined(SWIGCSHARP)
%rename(QlArray) Array;
#endif
class Array {
    #if defined(SWIGPYTHON)
    %rename(__len__)   size;
    #endif
  public:
    Array();
    Array(Size n, Real fill = 0.0);
    Array(const Array&);
    Size size() const;
    %extend {
        std::string __str__() {
            std::ostringstream out;
            out << *self;
            return out.str();
        }
        #if defined(SWIGPYTHON) || defined(SWIGR)
        Array __add__(const Array& a) {
            return Array(*self+a);
        }
        Array __sub__(const Array& a) {
            return Array(*self-a);
        }
        Array __mul__(Real a) {
            return Array(*self*a);
        }
        Real __mul__(const Array& a) {
            return QuantLib::DotProduct(*self,a);
        }
        Array __mul__(const Matrix& a) {
            return *self*a;
        }
        Array __div__(Real a) {
            return Array(*self/a);
        }
        #endif
        #if defined(SWIGPYTHON)
        Array __rmul__(Real a) {
            return Array(*self*a);
        }
        Array __getslice__(Integer i, Integer j) {
            Integer size_ = static_cast<Integer>(self->size());
            if (i<0)
                i = size_+i;
            if (j<0)
                j = size_+j;
            i = std::max(0,i);
            j = std::min(size_,j);
            Array tmp(j-i);
            std::copy(self->begin()+i,self->begin()+j,tmp.begin());
            return tmp;
        }
        void __setslice__(Integer i, Integer j, const Array& rhs) {
            Integer size_ = static_cast<Integer>(self->size());
            if (i<0)
                i = size_+i;
            if (j<0)
                j = size_+j;
            i = std::max(0,i);
            j = std::min(size_,j);
            QL_ENSURE(static_cast<Integer>(rhs.size()) == j-i,
                      "arrays are not resizable");
            std::copy(rhs.begin(),rhs.end(),self->begin()+i);
        }
        bool __bool__() {
            return (self->size() != 0);
        }
        #endif
        #if defined(SWIGPYTHON)
        Real __getitem__(Integer i) {
            Integer size_ = static_cast<Integer>(self->size());
            if (i>=0 && i<size_) {
                return (*self)[i];
            } else if (i<0 && -i<=size_) {
                return (*self)[size_+i];
            } else {
                throw std::out_of_range("array index out of range");
            }
        }
        void __setitem__(Integer i, Real x) {
            Integer size_ = static_cast<Integer>(self->size());
            if (i>=0 && i<size_) {
                (*self)[i] = x;
            } else if (i<0 && -i<=size_) {
                (*self)[size_+i] = x;
            } else {
                throw std::out_of_range("array index out of range");
            }
        }
        #elif defined(SWIGR)
        Real __getitem__(Integer i) {
            Integer size_ = static_cast<Integer>(self->size());
            if (i>=0 && i<size_) {
                return (*self)[i];
            } else {
                throw std::out_of_range("array index out of range");
            }
        }
        void __setitem__(Integer i, Real x) {
            Integer size_ = static_cast<Integer>(self->size());
            if (i>=0 && i<size_) {
                (*self)[i] = x;
            } else {
                throw std::out_of_range("array index out of range");
            }
        }
        #elif defined(SWIGCSHARP) || defined(SWIGJAVA)
        Real get(Size i) {
            if (i<self->size())
                return (*self)[i];
            else
                throw std::out_of_range("array index out of range");
        }
        void set(Size i, Real x) {
            if (i<self->size())
                (*self)[i] = x;
            else
                throw std::out_of_range("array index out of range");
        }
        #endif
    }
};


// matrix class
%{
using QuantLib::outerProduct;
using QuantLib::transpose;
using QuantLib::SVD;
%}

#if defined(SWIGPYTHON)
%{
class MatrixRow {
    Matrix::row_iterator begin_;
    Integer columns_;
  public:
    MatrixRow(Matrix::row_iterator begin, Size columns) : begin_(begin), columns_((Integer)columns) {}
    Real __getitem__(Integer i) {
        if (i >= 0 && i < columns_)
            return begin_[i];
        else if (i < 0 && -i <= columns_)
            return begin_[columns_+i];
        else
            throw std::out_of_range("matrix indexes out of range");
    }
    void __setitem__(Integer i, Real x) {
        if (i >= 0 && i < columns_)
            begin_[i] = x;
        else if (i < 0 && -i <= columns_)
            begin_[columns_+i] = x;
        else
            throw std::out_of_range("matrix indexes out of range");
    }
};
%}

class MatrixRow {
    MatrixRow();
  public:
    Real __getitem__(Integer i);
    void __setitem__(Integer i, Real x);
};
#endif

class Matrix {
  public:
    Matrix();
    Matrix(Size rows, Size columns, Real fill = 0.0);
    Matrix(const Matrix&);
    Size rows() const;
    Size columns() const;
    %extend {
        std::string __str__() {
            std::ostringstream out;
            out << *self;
            return out.str();
        }
        #if defined(SWIGPYTHON)
        Matrix __add__(const Matrix& m) {
            return *self+m;
        }
        Matrix __sub__(const Matrix& m) {
            return *self-m;
        }
        Matrix __mul__(Real x) {
            return *self*x;
        }
        Array __mul__(const Array& x) {
            return *self*x;
        }
        Matrix __mul__(const Matrix& x) {
            return *self*x;
        }
        Matrix __div__(Real x) {
            return *self/x;
        }
        #endif
        #if defined(SWIGPYTHON)
        MatrixRow __getitem__(Integer i) {
            Integer rows_ = static_cast<Integer>($self->rows());
            if (i >= 0 && i < $self->rows())
                return MatrixRow((*$self)[i], $self->columns());
            else if (i < 0 && -i <= rows_)
                return MatrixRow((*$self)[rows_+i], $self->columns());
            else
                throw std::out_of_range("matrix indexes out of range");
        }
        #elif defined(SWIGR)
        Real ref(Size i, Size j) {
            if (i < $self->rows() && j < $self->columns())
                return (*self)[i][j];
            else
                throw std::out_of_range("matrix indexes out of range");
        }
        void setitem(Size i, Size j, Real x) {
            if (i < $self->rows() && j < $self->columns())
                (*self)[i][j] = x;
            else
                throw std::out_of_range("matrix indexes out of range");
        }
        #elif defined(SWIGCSHARP) || defined(SWIGJAVA)
        Real get(Size i, Size j) {
            if (i < $self->rows() && j < $self->columns())
                return (*self)[i][j];
            else
                throw std::out_of_range("matrix indexes out of range");
        }
        void set(Size i, Size j, Real x) {
            if (i < $self->rows() && j < $self->columns())
                (*self)[i][j] = x;
            else
                throw std::out_of_range("matrix indexes out of range");
        }
        #endif
        #if defined(SWIGR)
        Array dataVector() {
            Size nrows = self->rows();
            Size ncols = self->columns();
            Size nelems = nrows * ncols;
            Array a(nelems);
            for (int i=0; i < nrows; i++)
                for (int j=0; j < ncols; j++)
                    a[j*nrows+i] = (*self)[i][j];
            return a;
        }
        #endif
        #if defined(SWIGPYTHON)
        Matrix __rmul__(Real x) {
            return x*(*self);
        }
        Array __rmul__(const Array& x) {
            return x*(*self);
        }
        Matrix __rmul__(const Matrix& x) {
            return x*(*self);
        }
        #endif
    }
};


// functions

%{
using QuantLib::inverse;
using QuantLib::pseudoSqrt;
using QuantLib::SalvagingAlgorithm;
%}

struct SalvagingAlgorithm {
    #if defined(SWIGPYTHON)
    %rename(NoAlgorithm) None;
    #endif
    enum Type { None, Spectral };
};

Matrix inverse(const Matrix& m);
Matrix transpose(const Matrix& m);
Matrix outerProduct(const Array& v1, const Array& v2);
Matrix pseudoSqrt(const Matrix& m, SalvagingAlgorithm::Type a);

class SVD {
  public:
    SVD(const Matrix&);
    const Matrix& U() const;
    const Matrix& V() const;
    Matrix S() const;
    const Array& singularValues() const;
};

%{
using QuantLib::BiCGstab;
using QuantLib::GMRES;
%}

#if defined(SWIGPYTHON)
%{
Array extractArray(PyObject* source, const std::string& methodName) {
    QL_ENSURE(source != NULL,
              "failed to call " + methodName + " on Python object");

    QL_ENSURE(source != Py_None, methodName + " returned None");
        
    Array* ptr;            
    const int err = SWIG_ConvertPtr(
        source, (void **) &ptr, SWIGTYPE_p_Array, SWIG_POINTER_EXCEPTION);

    if (err != 0) {
        Py_XDECREF(source);
        QL_FAIL("return type must be of type QuantLib Array in " 
            + methodName);
    }
    
    Array tmp(*ptr);          
    Py_XDECREF(source);
     
    return tmp;
}

class MatrixMultiplicationProxy {
  public:
    MatrixMultiplicationProxy(PyObject* matrixMult)
    : matrixMult_(matrixMult) {
        Py_XINCREF(matrixMult_);    
    }
    
    MatrixMultiplicationProxy(const MatrixMultiplicationProxy& p) 
    : matrixMult_(p.matrixMult_) {
        Py_XINCREF(matrixMult_);
    }
        
    MatrixMultiplicationProxy& operator=(const MatrixMultiplicationProxy& f) {
        if ((this != &f) && (matrixMult_ != f.matrixMult_)) {
            Py_XDECREF(matrixMult_);
            matrixMult_ = f.matrixMult_;
            Py_XINCREF(matrixMult_);
        }
        return *this;
    }
        
    ~MatrixMultiplicationProxy() {
        Py_XDECREF(matrixMult_);    
    }
    
    Array operator()(const Array& x) const {
        PyObject* pyArray = SWIG_NewPointerObj(
            SWIG_as_voidptr(&x), SWIGTYPE_p_Array, 0);
            
        PyObject* pyResult 
            = PyObject_CallFunction(matrixMult_, "O", pyArray);
        
        Py_XDECREF(pyArray);
        
        return extractArray(pyResult, "matrix multiplication");         
    }
    
  private:
    PyObject* matrixMult_;      
};
%}

class MatrixMultiplicationProxy {
  public:
    MatrixMultiplicationProxy(PyObject* matrixMult);
    
    %extend {
	    Array operator()(const Array& x) const {
	    	Array retVal = self->operator()(x);
	    	return retVal;
	    }
	}    
};

#elif defined(SWIGJAVA) || defined(SWIGCSHARP)

%{
class MatrixMultiplicationDelegate {
  public:
    virtual ~MatrixMultiplicationDelegate() {}
      
    virtual Array apply(const Array& x) const {
        QL_FAIL("implementation of MatrixMultiplicationDelegate.apply is missing");        
    }
};

class MatrixMultiplicationProxy {
  public:
    MatrixMultiplicationProxy(MatrixMultiplicationDelegate* delegate)
    : delegate_(delegate) {}
    
    Array operator()(const Array& x) const {
        Array retVal = delegate_->apply(x);        
        return retVal;
    }
               
  private:
      MatrixMultiplicationDelegate* const delegate_; 
};
%}

class MatrixMultiplicationDelegate {
  public:
    virtual ~MatrixMultiplicationDelegate();      
    virtual Array apply(const Array& x) const;
};

#endif

#if defined(SWIGPYTHON) || defined(SWIGJAVA) || defined(SWIGCSHARP)

%shared_ptr(BiCGstab)
class BiCGstab  {
  public:
    %extend {
        Array solve(const Array& b, const Array& x0 = Array()) const {
                return self->solve(b, x0).x; 
        }
#if defined(SWIGPYTHON)
        BiCGstab(const MatrixMultiplicationProxy& proxy, Size maxIter, Real relTol) {              
            return new BiCGstab(BiCGstab::MatrixMult(proxy), maxIter, relTol);                       
        }
        
        BiCGstab(const MatrixMultiplicationProxy& proxy, Size maxIter, Real relTol,
                 const MatrixMultiplicationProxy& preconditioner) {              
            return new BiCGstab(
                BiCGstab::MatrixMult(proxy), maxIter, relTol,
                BiCGstab::MatrixMult(preconditioner));                       
        }
#else
        BiCGstab(MatrixMultiplicationDelegate* delegate, Size maxIter, Real relTol) {
        	MatrixMultiplicationProxy proxy(delegate);          
            return new BiCGstab(BiCGstab::MatrixMult(proxy), maxIter, relTol);                       
        }
        
        BiCGstab(MatrixMultiplicationDelegate* delegate, Size maxIter, Real relTol,
                 MatrixMultiplicationDelegate* preconditioner) {              
        	MatrixMultiplicationProxy p1(delegate); 
        	MatrixMultiplicationProxy p2(preconditioner);
            return new BiCGstab(
                BiCGstab::MatrixMult(p1), maxIter, relTol, BiCGstab::MatrixMult(p2));                       
        }
#endif        
    }
};


%shared_ptr(GMRES)
class GMRES  {
  public:
    %extend {
        Array solve(const Array& b, const Array& x0 = Array()) const {
            return self->solve(b, x0).x;
        }
        Array solveWithRestart(
            Size restart, const Array& b, const Array& x0 = Array()) const {
            return self->solveWithRestart(restart, b, x0).x;
        }

#if defined(SWIGPYTHON)
        GMRES(const MatrixMultiplicationProxy& proxy, Size maxIter, Real relTol) {              
            return new GMRES(GMRES::MatrixMult(proxy), maxIter, relTol);                       
        }
        
        GMRES(const MatrixMultiplicationProxy& proxy, Size maxIter, Real relTol,
              const MatrixMultiplicationProxy& preconditioner) {              
            return new GMRES(
                GMRES::MatrixMult(proxy), maxIter, relTol,
                GMRES::MatrixMult(preconditioner));                       
        }
#else
        GMRES(MatrixMultiplicationDelegate* delegate, Size maxIter, Real relTol) {
        	MatrixMultiplicationProxy proxy(delegate);              
            return new GMRES(GMRES::MatrixMult(proxy), maxIter, relTol);                       
        }
        
        GMRES(MatrixMultiplicationDelegate* delegate, Size maxIter, Real relTol,
              const MatrixMultiplicationProxy& preconditioner) {
        	MatrixMultiplicationProxy p1(delegate); 
        	MatrixMultiplicationProxy p2(preconditioner);                                      
            return new GMRES(
                GMRES::MatrixMult(p1), maxIter, relTol, GMRES::MatrixMult(p2));                       
        }
#endif        
    }    
};

#endif 

%{
using QuantLib::close;
using QuantLib::close_enough;
%}

bool close(Real x, Real y);
bool close(Real x, Real y, Size n);

bool close_enough(Real x, Real y);
bool close_enough(Real x, Real y, Size n);

#endif
