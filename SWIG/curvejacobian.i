/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 Kyrylo Protsenko

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
*/

#ifndef quantlib_curve_jacobian_i
#define quantlib_curve_jacobian_i

%include common.i
%include types.i
%include vectors.i
%include linearalgebra.i
%include termstructures.i

namespace std {
    %template(YieldTermStructureVector) vector<ext::shared_ptr<YieldTermStructure> >;
    %template(ArrayVector) vector<Array>;
}

%{
#include <ql/experimental/termstructures/jacobian/curvejacobiangraph.hpp>

// The graph needs concrete curve types, while wrappers expose base pointers.
// Each supported type registers a downcasting adder.
class CurveJacobianGraphProxy;

typedef bool (*CurveJacobianAdder)(CurveJacobianGraphProxy&,
                                   const ext::shared_ptr<YieldTermStructure>&);

inline std::vector<CurveJacobianAdder>& curveJacobianAdders() {
    static std::vector<CurveJacobianAdder> adders;
    return adders;
}

/* Keeps curves alive and records their order for vector results. */
class CurveJacobianGraphProxy {
  public:
    explicit CurveJacobianGraphProxy(bool errorOnIncomplete = false)
    : graph_(errorOnIncomplete) {}

    bool isComplete() const {
        return graph_.isComplete();
    }

    void add(const ext::shared_ptr<YieldTermStructure>& curve) {
        QL_REQUIRE(curve, "null curve");
        for (auto adder : curveJacobianAdders())
            if (adder(*this, curve))
                return;
        QL_FAIL("the given curve type is not supported by "
                "CurveJacobianGraph.add(); supported curves must provide "
                "bootstrap Jacobians or expose their underlying curve to C++");
    }

    template <class Curve>
    void addCurve(const ext::shared_ptr<Curve>& curve) {
        graph_.add(curve);
        if constexpr (QuantLib::detail::supportsCurveJacobianNode<Curve>) {
            for (auto& c : curves_) {
                if (c.get() == curve.get())
                    return;  // graph_.add() already replaced it
            }
            curves_.push_back(curve);
        } else {
            for (auto& c : inspectedCurves_)
                if (c.get() == curve.get())
                    return;
            inspectedCurves_.push_back(curve);
        }
    }

    const std::vector<ext::shared_ptr<YieldTermStructure> >& curves() const {
        return curves_;
    }

    Matrix crossJacobian(const ext::shared_ptr<YieldTermStructure>& of,
                         const ext::shared_ptr<YieldTermStructure>& withRespectTo,
                         std::vector<bool>* analyticEquations = nullptr) const {
        QL_REQUIRE(of && withRespectTo, "null curve");
        return graph_.crossJacobian(*of, *withRespectTo, analyticEquations);
    }

    Matrix inverseJacobian(const ext::shared_ptr<YieldTermStructure>& of,
                           const ext::shared_ptr<YieldTermStructure>& withRespectTo,
                           std::vector<bool>* analyticEquations = nullptr) const {
        QL_REQUIRE(of && withRespectTo, "null curve");
        return graph_.inverseJacobian(*of, *withRespectTo, analyticEquations);
    }

    Matrix zeroNodeJacobian(const ext::shared_ptr<YieldTermStructure>& curve,
                            std::vector<bool>* analyticDerivatives = nullptr) const {
        QL_REQUIRE(curve, "null curve");
        return graph_.zeroNodeJacobian(*curve, analyticDerivatives);
    }

    Matrix nodeZeroJacobian(const ext::shared_ptr<YieldTermStructure>& curve,
                            std::vector<bool>* analyticDerivatives = nullptr) const {
        QL_REQUIRE(curve, "null curve");
        return graph_.nodeZeroJacobian(*curve, analyticDerivatives);
    }

    std::map<const YieldTermStructure*, Array> riskInput(
        const std::vector<ext::shared_ptr<YieldTermStructure> >& curves,
              const std::vector<Array>& nodeRisk) const {
        QL_REQUIRE(curves.size() == nodeRisk.size(),
                   "the number of curves (" << curves.size() <<
                   ") does not match the number of node-risk vectors (" <<
                   nodeRisk.size() << ")");
        std::map<const YieldTermStructure*, Array> input;
        for (Size i = 0; i < curves.size(); ++i) {
            QL_REQUIRE(curves[i], "null curve");
            auto [entry, inserted] = input.emplace(curves[i].get(), nodeRisk[i]);
            if (!inserted) {
                QL_REQUIRE(entry->second.size() == nodeRisk[i].size(),
                           "node-risk vectors for the same curve have different sizes");
                entry->second += nodeRisk[i];
            }
        }
        return input;
    }

    //! converts a risk map to registration order
    std::vector<Array> ordered(
        const std::map<const YieldTermStructure*, Array>& risk) const {
        std::vector<Array> result;
        result.reserve(curves_.size());
        for (const auto& c : curves_)
            result.push_back(risk.at(c.get()));
        return result;
    }

    std::map<const YieldTermStructure*, Array> parRiskMap(
        const std::vector<ext::shared_ptr<YieldTermStructure> >& curves,
               const std::vector<Array>& nodeRisk,
               std::vector<bool>* analyticEquations = nullptr) const {
        return graph_.parRisk(riskInput(curves, nodeRisk), analyticEquations);
    }

    std::map<const YieldTermStructure*, Array> zeroRiskMap(
        const std::vector<ext::shared_ptr<YieldTermStructure> >& curves,
                const std::vector<Array>& nodeRisk,
                std::vector<bool>* analyticEquations = nullptr) const {
        return graph_.zeroRisk(riskInput(curves, nodeRisk), analyticEquations);
    }

    void propagateMap(const std::vector<ext::shared_ptr<YieldTermStructure> >& curves,
                      const std::vector<Array>& nodeRisk,
                      std::map<const YieldTermStructure*, Array>* zeroRisk,
                      std::map<const YieldTermStructure*, Array>* parRisk,
                      std::vector<bool>* analyticEquations = nullptr) const {
        graph_.propagateNodeRisk(riskInput(curves, nodeRisk),
                                 zeroRisk, parRisk, analyticEquations);
    }

  private:
    QuantLib::CurveJacobianGraph graph_;
    std::vector<ext::shared_ptr<YieldTermStructure> > curves_;
    std::vector<ext::shared_ptr<YieldTermStructure> > inspectedCurves_;
};

template <class Curve>
bool addCurveToJacobianGraph(CurveJacobianGraphProxy& graph,
                             const ext::shared_ptr<YieldTermStructure>& curve) {
    if (auto c = ext::dynamic_pointer_cast<Curve>(curve)) {
        graph.addCurve(c);
        return true;
    }
    return false;
}

template <class Curve>
bool registerCurveJacobianAdder() {
    curveJacobianAdders().push_back(&addCurveToJacobianGraph<Curve>);
    return true;
}
%}

/* Registers a piecewise curve type when the module loads. */
%define export_curve_to_jacobian_graph(Name)
%{
bool Name ## _registeredWithCurveJacobianGraph = registerCurveJacobianAdder<Name>();
%}
%enddef

/* Registers an inspectable derived curve without a bootstrap block. */
%define export_derived_curve_to_jacobian_graph(Name,Curve)
%{
bool Name ## _registeredAsDerivedWithCurveJacobianGraph =
    registerCurveJacobianAdder<Curve>();
%}
%enddef

/* Adds Jacobian accessors to an exported piecewise curve. */
%define export_curve_jacobian_methods
%extend {
    /*! Jacobian of alive helper quotes with respect to free curve nodes.

        Rows follow the helpers and columns follow data()[1:]. Other curves'
        nodes are held fixed. analyticEquations identifies helper equations
        calculated analytically. The
        remaining rows use numerical differentiation.
    */
    Matrix jacobian() {
        return self->jacobian();
    }
    Matrix jacobian(std::vector<bool>& analyticEquations) {
        return self->jacobian(&analyticEquations);
    }

    /*! Jacobian of free curve nodes with respect to helper quotes.

        This is the inverse of jacobian() for a standalone curve. For a
        MultiCurve member, it includes feedback through the group. Columns
        contain this curve's quotes first, followed by the other members in
        registration order.
    */
    Matrix inverseJacobian() {
        return self->inverseJacobian();
    }
    Matrix inverseJacobian(std::vector<bool>& analyticEquations) {
        return self->inverseJacobian(&analyticEquations);
    }
}
%enddef

%rename(CurveJacobianGraph) CurveJacobianGraphProxy;
class CurveJacobianGraphProxy {
  public:
    //! cross-curve Jacobians for registered bootstrapped curves
    /*! Tracks curve dependencies and propagates node sensitivities to helper
        quotes. Rows follow alive helpers. Columns follow free curve nodes and
        exclude the reference-date value.
    */
    explicit CurveJacobianGraphProxy(bool errorOnIncomplete = false);

    //! whether all reported curve dependencies are registered
    bool isComplete() const;

    //! registers a curve and replaces any existing entry for it
    void add(const ext::shared_ptr<YieldTermStructure>& curve);

    //! the registered curves, in registration order
    std::vector<ext::shared_ptr<YieldTermStructure> > curves() const;

    %extend {
        void add(const Handle<YieldTermStructure>& curve) {
            QL_REQUIRE(!curve.empty(), "empty curve handle");
            self->add(curve.currentLink());
        }

        /*! Jacobian of the first curve's helper quotes with respect to the
            second curve's nodes. Other curve nodes are held fixed. Equal
            curves give the local Jacobian. analyticEquations identifies
            helper equations that did not require numerical differentiation.
        */
        Matrix crossJacobian(const ext::shared_ptr<YieldTermStructure>& of,
                             const ext::shared_ptr<YieldTermStructure>& withRespectTo) {
            return self->crossJacobian(of, withRespectTo);
        }
        Matrix crossJacobian(const ext::shared_ptr<YieldTermStructure>& of,
                             const ext::shared_ptr<YieldTermStructure>& withRespectTo,
                             std::vector<bool>& analyticEquations) {
            return self->crossJacobian(of, withRespectTo,
                                       &analyticEquations);
        }

        /*! Block of the inverse Jacobian for the registered curve system.
            Rows are the first curve's nodes and columns are the second
            curve's helper quotes. Dependency cycles are included.
        */
        Matrix inverseJacobian(const ext::shared_ptr<YieldTermStructure>& of,
                               const ext::shared_ptr<YieldTermStructure>& withRespectTo) {
            return self->inverseJacobian(of, withRespectTo);
        }
        Matrix inverseJacobian(const ext::shared_ptr<YieldTermStructure>& of,
                               const ext::shared_ptr<YieldTermStructure>& withRespectTo,
                               std::vector<bool>& analyticEquations) {
            return self->inverseJacobian(of, withRespectTo,
                                         &analyticEquations);
        }

        /*! Jacobian of continuously compounded zero rates at a curve's node
            dates with respect to its free stored nodes.
        */
        Matrix zeroNodeJacobian(
                const ext::shared_ptr<YieldTermStructure>& curve) {
            return self->zeroNodeJacobian(curve);
        }
        Matrix zeroNodeJacobian(
                const ext::shared_ptr<YieldTermStructure>& curve,
                std::vector<bool>& analyticDerivatives) {
            return self->zeroNodeJacobian(curve, &analyticDerivatives);
        }

        /*! Jacobian of a curve's free stored nodes with respect to
            continuously compounded zero rates at its node dates.
        */
        Matrix nodeZeroJacobian(
                const ext::shared_ptr<YieldTermStructure>& curve) {
            return self->nodeZeroJacobian(curve);
        }
        Matrix nodeZeroJacobian(
                const ext::shared_ptr<YieldTermStructure>& curve,
                std::vector<bool>& analyticDerivatives) {
            return self->nodeZeroJacobian(curve, &analyticDerivatives);
        }

        /*! Converts curve-node sensitivities to helper-quote sensitivities.
            Input curves must be registered. Each risk vector must match its
            curve's free nodes. Repeated curves are summed. Results follow
            curves().
        */
        std::vector<Array> parRisk(
            const std::vector<ext::shared_ptr<YieldTermStructure> >& curves,
                const std::vector<Array>& nodeRisk) {
            return self->ordered(self->parRiskMap(curves, nodeRisk));
        }

        std::vector<Array> parRisk(
            const std::vector<ext::shared_ptr<YieldTermStructure> >& curves,
                const std::vector<Array>& nodeRisk,
                std::vector<bool>& analyticEquations) {
            return self->ordered(
                self->parRiskMap(curves, nodeRisk, &analyticEquations));
        }

        /*! Propagates node sensitivities through curve dependencies. Zero risk
            fixes each curve while the others rebootstrap. The curve's bootstrap
            block converts zero risk to par risk only for acyclic dependencies.
            Cyclic components require a reduced solve.
            Inputs and results follow parRisk().
        */
        std::vector<Array> zeroRisk(
            const std::vector<ext::shared_ptr<YieldTermStructure> >& curves,
                 const std::vector<Array>& nodeRisk) {
            return self->ordered(self->zeroRiskMap(curves, nodeRisk));
        }

        std::vector<Array> zeroRisk(
            const std::vector<ext::shared_ptr<YieldTermStructure> >& curves,
                 const std::vector<Array>& nodeRisk,
                  std::vector<bool>& analyticEquations) {
            return self->ordered(
                self->zeroRiskMap(curves, nodeRisk, &analyticEquations));
        }

        //! zero risk on one registered curve
        Array zeroRisk(const std::vector<ext::shared_ptr<YieldTermStructure> >& curves,
                       const std::vector<Array>& nodeRisk,
                       const ext::shared_ptr<YieldTermStructure>& onCurve) {
            QL_REQUIRE(onCurve, "null curve");
            std::map<const YieldTermStructure*, Array> risk =
                self->zeroRiskMap(curves, nodeRisk);
            auto i = risk.find(onCurve.get());
            QL_REQUIRE(i != risk.end(),
                       "the given curve was not added to the graph");
            return i->second;
        }

        /*! Computes zero and par risk in one propagation. Results contain one
            entry per registered curve and follow curves().
        */
        void propagateRisk(
                const std::vector<ext::shared_ptr<YieldTermStructure> >& curves,
                const std::vector<Array>& nodeRisk,
                std::vector<Array>& zeroRisk,
                std::vector<Array>& parRisk) {
            std::map<const YieldTermStructure*, Array> z, p;
            self->propagateMap(curves, nodeRisk, &z, &p);
            zeroRisk = self->ordered(z);
            parRisk = self->ordered(p);
        }

        void propagateRisk(
                const std::vector<ext::shared_ptr<YieldTermStructure> >& curves,
                const std::vector<Array>& nodeRisk,
                std::vector<Array>& zeroRisk,
                std::vector<Array>& parRisk,
                std::vector<bool>& analyticEquations) {
            std::map<const YieldTermStructure*, Array> z, p;
            self->propagateMap(curves, nodeRisk, &z, &p,
                               &analyticEquations);
            zeroRisk = self->ordered(z);
            parRisk = self->ordered(p);
        }

        //! par risk on one registered curve
        Array parRisk(const std::vector<ext::shared_ptr<YieldTermStructure> >& curves,
                      const std::vector<Array>& nodeRisk,
                      const ext::shared_ptr<YieldTermStructure>& onCurve) {
            QL_REQUIRE(onCurve, "null curve");
            std::map<const YieldTermStructure*, Array> risk =
                self->parRiskMap(curves, nodeRisk);
            auto i = risk.find(onCurve.get());
            QL_REQUIRE(i != risk.end(),
                       "the given curve was not added to the graph");
            return i->second;
        }
    }
};

#endif
