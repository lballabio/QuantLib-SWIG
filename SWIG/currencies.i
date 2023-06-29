
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006, 2008 StatPro Italia srl
 Copyright (C) 2005 Johan Witters
 Copyright (C) 2023 Skandinaviska Enskilda Banken AB (publ)

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

#ifndef quantlib_currencies_i
#define quantlib_currencies_i

%include common.i
%include types.i
%include rounding.i

// currency objects

%{
using QuantLib::Currency;
using QuantLib::Money;
%}

class Currency {
  public:
    Currency() = default;
    Currency(const std::string& name,
             const std::string& code,
             Integer numericCode,
             const std::string& symbol,
             const std::string& fractionSymbol,
             Integer fractionsPerUnit,
             const Rounding& rounding,
             const std::string& formatString,
             const Currency& triangulationCurrency = Currency());
    
    const std::string& name() const;
    const std::string& code() const;
    Integer numericCode() const;
    const std::string& symbol() const;
    const std::string& fractionSymbol() const;
    Integer fractionsPerUnit() const;
    const Rounding& rounding() const;
    std::string format() const;
    bool empty() const;
    const Currency& triangulationCurrency() const;
    %extend {
        std::string __str__() {
            return self->name();
        }
        #if defined(SWIGPYTHON) || defined(SWIGJAVA)
        bool __eq__(const Currency& other) {
            return (*self) == other;
        }
        bool __ne__(const Currency& other) {
            return (*self) != other;
        }
        #endif
        #if defined(SWIGPYTHON)
        Money operator*(Decimal x) {
            return *self*x;
        }
        Money __rmul__(Decimal x) {
            return *self*x;
        }
        bool __nonzero__() {
            return !self->empty();
        }
        bool __bool__() {
            return !self->empty();
        }
        #endif
    }
    #if defined(SWIGPYTHON)
    %pythoncode %{
    def __hash__(self):
        return hash(self.name())
    %}
    #endif
};


namespace QuantLib {
class AEDCurrency : public Currency {};
class AOACurrency : public Currency {};
class ARSCurrency : public Currency {};
class ATSCurrency : public Currency {};
class AUDCurrency : public Currency {};
class BDTCurrency : public Currency {};
class BEFCurrency : public Currency {};
class BHDCurrency : public Currency {};
class BGLCurrency : public Currency {};
class BGNCurrency : public Currency {};
class BRLCurrency : public Currency {};
class BWPCurrency : public Currency {};
class BYRCurrency : public Currency {};
class CADCurrency : public Currency {};
class CHFCurrency : public Currency {};
class CLFCurrency : public Currency {};
class CLPCurrency : public Currency {};
class CNHCurrency : public Currency {};
class CNYCurrency : public Currency {};
class COPCurrency : public Currency {};
class COUCurrency : public Currency {};
class CYPCurrency : public Currency {};
class CZKCurrency : public Currency {};
class DEMCurrency : public Currency {};
class DKKCurrency : public Currency {};
class EEKCurrency : public Currency {};
class EGPCurrency : public Currency {};
class ESPCurrency : public Currency {};
class ETBCurrency : public Currency {};
class EURCurrency : public Currency {};
class FIMCurrency : public Currency {};
class FRFCurrency : public Currency {};
class GELCurrency : public Currency {};
class GBPCurrency : public Currency {};
class GHSCurrency : public Currency {};
class GRDCurrency : public Currency {};
class HKDCurrency : public Currency {};
class HRKCurrency : public Currency {};
class HUFCurrency : public Currency {};
class IDRCurrency : public Currency {};
class IEPCurrency : public Currency {};
class ILSCurrency : public Currency {};
class INRCurrency : public Currency {};
class IQDCurrency : public Currency {};
class IRRCurrency : public Currency {};
class ISKCurrency : public Currency {};
class ITLCurrency : public Currency {};
class JODCurrency : public Currency {};
class JPYCurrency : public Currency {};
class KESCurrency : public Currency {};
class KRWCurrency : public Currency {};
class KWDCurrency : public Currency {};
class KZTCurrency : public Currency {};
class LKRCurrency : public Currency {};
class LTLCurrency : public Currency {};
class LUFCurrency : public Currency {};
class LVLCurrency : public Currency {};
class MADCurrency : public Currency {};
class MTLCurrency : public Currency {};
class MURCurrency : public Currency {};
class MXNCurrency : public Currency {};
class MXVCurrency : public Currency {};
class MYRCurrency : public Currency {};
class NGNCurrency : public Currency {};
class NLGCurrency : public Currency {};
class NOKCurrency : public Currency {};
class NPRCurrency : public Currency {};
class NZDCurrency : public Currency {};
class OMRCurrency : public Currency {};
class PEHCurrency : public Currency {};
class PEICurrency : public Currency {};
class PENCurrency : public Currency {};
class PHPCurrency : public Currency {};
class PKRCurrency : public Currency {};
class PLNCurrency : public Currency {};
class PTECurrency : public Currency {};
class QARCurrency : public Currency {};
class ROLCurrency : public Currency {};
class RONCurrency : public Currency {};
class RSDCurrency : public Currency {};
class RUBCurrency : public Currency {};
class SARCurrency : public Currency {};
class SEKCurrency : public Currency {};
class SGDCurrency : public Currency {};
class SITCurrency : public Currency {};
class SKKCurrency : public Currency {};
class THBCurrency : public Currency {};
class TNDCurrency : public Currency {};
class TRLCurrency : public Currency {};
class TRYCurrency : public Currency {};
class TTDCurrency : public Currency {};
class TWDCurrency : public Currency {};
class UAHCurrency : public Currency {};
class UGXCurrency : public Currency {};
class USDCurrency : public Currency {};
class UYUCurrency : public Currency {};
class VEBCurrency : public Currency {};
class VNDCurrency : public Currency {};
class XOFCurrency : public Currency {};
class ZARCurrency : public Currency {};
class ZMWCurrency : public Currency {};

// Crypto
class BCHCurrency : public Currency {};
class BTCCurrency : public Currency {};
class DASHCurrency : public Currency {};
class ETCCurrency : public Currency {};
class ETHCurrency : public Currency {};
class LTCCurrency : public Currency {};
class XRPCurrency : public Currency {};
class ZECCurrency : public Currency {};
}


#endif
