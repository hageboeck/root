// Authors: Stephan Hageboeck, CERN; Andrea Sciandra, SCIPP-UCSC/Atlas; Nov 2020

/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2018, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_BIN_SAMPLING__PDF
#define ROO_BIN_SAMPLING__PDF

#include "RooAbsReal.h"
#include "RooTemplateProxy.h"
#include "RooAbsPdf.h"

#include "Math/Integrator.h"

#include <memory>

class RooBinSamplingPdf : public RooAbsPdf {
public:

  RooBinSamplingPdf() { };
  RooBinSamplingPdf(const char *name, const char *title, RooAbsRealLValue& observable, RooAbsPdf& inputPdf,
      double epsilon = 1.E-5);
  virtual ~RooBinSamplingPdf() {};

  RooBinSamplingPdf(const RooBinSamplingPdf& other, const char* name = 0) :
    RooAbsPdf(other, name),
    _pdf("inputPdf", this, other._pdf) { }

  virtual TObject* clone(const char* newname) const override {
    return new RooBinSamplingPdf(*this, newname);
  }

  // Analytical Integration handling
  Bool_t forceAnalyticalInt(const RooAbsArg& dep) const override {
    return _pdf->forceAnalyticalInt(dep);
  }

  /// Since contained PDF is already normalised, this always returns true.
  bool selfNormalized() const override { return true; }


  // Internal toy generation. Since our _func is not a PDF (if it is, it doesn't make sense to use this wrapper),
  // we cannot do anything.
  /// Get specialised generator. Since the underlying function is not a PDF, this will always return zero.
  Int_t getGenerator(const RooArgSet& directVars, RooArgSet& generateVars, bool staticInitOK = true) const override {
    return _pdf->getGenerator(directVars, generateVars, staticInitOK);
  }
  void initGenerator(Int_t code) override { _pdf->initGenerator(code); }
  void generateEvent(Int_t code) override { _pdf->generateEvent(code); }
  Bool_t isDirectGenSafe(const RooAbsArg& arg) const override { return _pdf->isDirectGenSafe(arg); }


  // Hints for optimized brute-force sampling
  Int_t getMaxVal(const RooArgSet& vars) const override { return _pdf->getMaxVal(vars); }
  Double_t maxVal(Int_t code) const override { return _pdf->maxVal(code); }
  Int_t minTrialSamples(const RooArgSet& arGenObs) const override { return _pdf->minTrialSamples(arGenObs); }

  // Plotting and binning hints
  Bool_t isBinnedDistribution(const RooArgSet& obs) const override { return _pdf->isBinnedDistribution(obs); }
  std::list<Double_t>* binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const override {
    return _pdf->binBoundaries(obs, xlo, xhi);
  }
  std::list<Double_t>* plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const override {
    return _pdf->plotSamplingHint(obs, xlo, xhi);
  }

  ROOT::Math::IntegratorOneDim& integrator();


private:
  double evaluate() const override;
  RooSpan<double> evaluateSpan(BatchHelpers::RunContext& evalData, const RooArgSet* normSet) const override;
  RooSpan<const double> binBoundaries() const;

  struct BinSamplingInfo {
  //  std::map<std::pair<std::size_t, std::size_t>, std::vector<double>> xValues; //cache for x values for high-resolution sampling
    std::vector<double> binBoundaries;
  };

  RooTemplateProxy<RooAbsPdf> _pdf;
  RooTemplateProxy<RooAbsRealLValue> _observable;
  std::unique_ptr<ROOT::Math::IntegratorOneDim> _integrator; /// Integrator used to sample bins.
  mutable std::unique_ptr<BinSamplingInfo> _binSamplingInfo; //! Workspace to store data for bin sampling
  mutable RooArgSet* _normSetForIntegrator{nullptr}; //!

  ClassDefOverride(RooBinSamplingPdf,0)
};

#endif
