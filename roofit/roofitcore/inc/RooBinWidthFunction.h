// Author Stephan Hageboeck, CERN, 6/2020
/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2020, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROOFIT_ROOFITCORE_INC_BINWIDTHFUNCTION_H_
#define ROOFIT_ROOFITCORE_INC_BINWIDTHFUNCTION_H_

#include "RooAbsReal.h"
#include "RooTemplateProxy.h"
#include "RooHistFunc.h"

#include "RooDataHist.h"

class RooBinWidthFunction : public RooAbsReal {
public:
//  RooBinWidthFunction() :
//    _histFunc("HistFuncForBinWidth", "Handle to a RooHistFunc, whose bin volumes should be returned.", this,
//        /*valueServer=*/true, /*shapeServer=*/true) { }

  RooBinWidthFunction(const char* name, const char* title, const RooHistFunc& histFunc) :
    RooAbsReal(name, title),
    _histFunc("HistFuncForBinWidth", "Handle to a RooHistFunc, whose bin volumes should be returned.", this, histFunc,
      /*valueServer=*/true, /*shapeServer=*/true) { }
  RooBinWidthFunction(const RooBinWidthFunction& other, const char* newname = nullptr) :
    RooAbsReal(other, newname),
    _histFunc("HistFuncForBinWidth", this, other._histFunc) { }

  virtual TObject* clone(const char* newname = nullptr) const override {
    return new RooBinWidthFunction(*this, newname);
  }

  double evaluate() const override {
    return 1;

    const RooDataHist& dataHist = _histFunc->dataHist();
    auto idx = _histFunc->getBin();

//    return 1./dataHist.binVolumes()[idx];
//    return dataHist.binVolumes()[idx] / _histFunc->totVolume();
    return dataHist.binVolumes()[idx];
  }

private:
  RooTemplateProxy<const RooHistFunc> _histFunc;

  ClassDefOverride(RooBinWidthFunction, 1)
};

#endif
