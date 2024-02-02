
/*
 Copyright (C) 2004, 2005, 2009 StatPro Italia srl

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

#ifndef quantlib_settings_i
#define quantlib_settings_i

%include common.i
%include date.i

%{
using QuantLib::Settings;

#if defined(SWIGPYTHON)
Date* Settings_evaluationDate_get(Settings* self) {
    return new Date(self->evaluationDate());
}
void Settings_evaluationDate_set(Settings* self, const Date* d) {
    self->evaluationDate() = *d;
}
bool Settings_enforcesTodaysHistoricFixings_get(Settings* self) {
    return self->enforcesTodaysHistoricFixings();
}
void Settings_enforcesTodaysHistoricFixings_set(Settings* self, bool b) {
    self->enforcesTodaysHistoricFixings() = b;
}
bool Settings_includeReferenceDateEvents_get(Settings* self) {
    return self->includeReferenceDateEvents();
}
void Settings_includeReferenceDateEvents_set(Settings* self, bool b) {
    self->includeReferenceDateEvents() = b;
}
ext::optional<bool> Settings_includeTodaysCashFlows_get(Settings* self) {
    return self->includeTodaysCashFlows();
}
void Settings_includeTodaysCashFlows_set(Settings* self, ext::optional<bool> b) {
    self->includeTodaysCashFlows() = b;
}
#endif
%}

class Settings {
  private:
    Settings();
  public:
    static Settings& instance();
    void anchorEvaluationDate();
    void resetEvaluationDate();
    %extend {
        Date getEvaluationDate() {
            #if defined(SWIGPYTHON)
            cpp_deprecate_feature(getEvaluationDate, evaluationDate);
            #endif
            return self->evaluationDate();
        }
        void setEvaluationDate(const Date& d) {
            #if defined(SWIGPYTHON)
            cpp_deprecate_feature(setEvaluationDate, evaluationDate);
            #endif
            self->evaluationDate() = d;
        }
        bool getEnforcesTodaysHistoricFixings() {
            #if defined(SWIGPYTHON)
            cpp_deprecate_feature(getEnforcesTodaysHistoricFixings, enforcesTodaysHistoricFixings);
            #endif
            return self->enforcesTodaysHistoricFixings();
        }
        void setEnforcesTodaysHistoricFixings(bool b) {
            #if defined(SWIGPYTHON)
            cpp_deprecate_feature(setEnforcesTodaysHistoricFixings, enforcesTodaysHistoricFixings);
            #endif
            self->enforcesTodaysHistoricFixings() = b;
        }
        #if defined(SWIGPYTHON)
        %newobject evaluationDate;
        Date evaluationDate;
        bool enforcesTodaysHistoricFixings;
        bool includeReferenceDateEvents;
        ext::optional<bool> includeTodaysCashFlows;
        #else
        void includeReferenceDateEvents(bool b) {
            self->includeReferenceDateEvents() = b;
        }
        void includeTodaysCashFlows(bool b) {
            self->includeTodaysCashFlows() = b;
        }
        #endif
    }
};

#if defined(SWIGPYTHON)
%rename(SavedSettings) _SavedSettings;
%inline %{
class _SavedSettings {
    ext::optional<QuantLib::SavedSettings> saved_;
  public:
    void __enter__() {
        saved_.emplace();
    }
    void __exit__(PyObject*, PyObject*, PyObject*) {
        saved_.reset();
    }
};
%}
#endif

#endif
