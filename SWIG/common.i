
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl

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

#ifndef quantlib_common_i
#define quantlib_common_i

%{
    namespace QuantLib { namespace ext {} }
    namespace ext = QuantLib::ext;
%}
#define SWIG_SHARED_PTR_NAMESPACE ext

%include stl.i
%include exception.i
%include boost_shared_ptr.i

%define QL_TYPECHECK_BOOL       7210    %enddef

%{
// This is necessary to avoid compile failures on 
// GCC 4
// see http://svn.boost.org/trac/boost/ticket/1793

#if defined(NDEBUG)
#define BOOST_DISABLE_ASSERTS 1
#endif

#include <boost/algorithm/string/case_conv.hpp>
%}

#if defined(SWIGPYTHON)
%typemap(in) ext::optional<bool> %{
	if($input == Py_None)
		$1 = ext::nullopt;
	else if ($input == Py_True)
		$1 = true;
	else
		$1 = false;
%}
%typecheck (QL_TYPECHECK_BOOL) ext::optional<bool> {
if (PyBool_Check($input) || Py_None == $input) 
	$1 = 1;
else
	$1 = 0;
}
#else
#if defined(SWIGCSHARP)
%typemap(cscode) ext::optional<bool> %{
    public static implicit operator OptionalBool(bool b) => new OptionalBool(b);
%}
#endif
namespace ext {
    template<class T>
    class optional {
      public:
        optional(T t);
    };
}
%template(OptionalBool) ext::optional<bool>;
#endif

%{
// generally useful classes
using QuantLib::Error;
using QuantLib::Handle;
using QuantLib::RelinkableHandle;
%}

namespace ext {

    %extend shared_ptr {
        T* operator->() {
            return (*self).operator->();
        }
        #if defined(SWIGPYTHON)
        bool __nonzero__() {
            return !!(*self);
        }
        bool __bool__() {
            return !!(*self);
        }
        #else
        bool isNull() {
            return !(*self);
        }
        #endif
    }
}


template <class T>
class Handle {
  public:
  Handle(const ext::shared_ptr<T>& = ext::shared_ptr<T>());
    ext::shared_ptr<T> operator->();
    ext::shared_ptr<T> currentLink();
    #if defined(SWIGPYTHON)
    %extend {
        bool __nonzero__() {
            return !self->empty();
        }
        bool __bool__() {
            return !self->empty();
        }
    }
    #else
    bool empty();
    #endif
};

template <class T>
class RelinkableHandle : public Handle<T> {
  public:
    RelinkableHandle(const ext::shared_ptr<T>& = ext::shared_ptr<T>());
    void linkTo(const ext::shared_ptr<T>&);
    %extend {
        // could be defined in C++ class, added here in the meantime
        void reset() {
            self->linkTo(ext::shared_ptr<T>());
        }
    }
};

%define swigr_list_converter(ContainerRType,
                             ContainerCType, ElemCType)
#if defined(SWIGR)
%Rruntime %{
setMethod('print', 'ContainerCType',
function(x) print(as(x, "ElemCType")))

setAs("ContainerCType", "ElemCType",
function(from) {if (from$size()) from[1:from$size()] else NULL} )

setAs("ElemCType", "ContainerCType",
function(from) { a <- ContainerRType(length(from));
sapply(1:length(from), function(n) {
a[n] <- from[n] } )
a
})
%}
#endif
%enddef


%define deprecate_feature(OldName, NewName)
#if defined(SWIGPYTHON)
%pythoncode %{
def OldName(*args, **kwargs):
    from warnings import warn
    warn('%s is deprecated; use %s' % (OldName.__name__, NewName.__name__))
    return NewName(*args, **kwargs)
%}
#endif
%enddef


#endif
