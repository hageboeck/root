// Authors: Stephan Hageboeck, CERN; Andrea Sciandra, SCIPP-UCSC/Atlas; Nov 2020

/*****************************************************************************
 * RooFit
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2019, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/**
 * \class RooBinSamplingPdf
 * The RooBinSamplingPdf TODO
 */


#include "RooBinSamplingPdf.h"

#include "RooHelpers.h"
#include "RooRealBinding.h"
#include "BatchHelpers.h"
#include "RunContext.h"

#include "Math/Integrator.h"

#include <algorithm>

////////////////////////////////////////////////////////////////////////////////
/// Construct a new RooBinSamplingPdf.
/// \param[in] name A name to identify this object.
/// \param[in] title Title (for e.g. plotting)
/// \param[in] observable Observable to integrate over (the one that is binned).
/// \param[in] inputPdf A PDF whose bins should be sampled with higher precision.
/// \param[in] epsilon Relative epsilon for the integrator, which is used to sample the bins.
RooBinSamplingPdf::RooBinSamplingPdf(const char *name, const char *title, RooAbsRealLValue& observable,
    RooAbsPdf& inputPdf, double epsilon) :
      RooAbsPdf(name, title),
      _pdf("inputPdf", "Function to be converted into a PDF", this, inputPdf),
      _observable("observable", "Observable to integrator over", this, observable, true, true),
      _integrator(new ROOT::Math::IntegratorOneDim(ROOT::Math::IntegrationOneDim::kADAPTIVE, -1., epsilon)) {
  if (!_pdf->dependsOn(*_observable)) {
    throw std::invalid_argument(std::string("RooBinSamplingPDF(") + GetName()
        + "): The PDF " + _pdf->GetName() + " needs to depend on the observable "
        + _observable->GetName());
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Integrate the PDF over the current bin.
double RooBinSamplingPdf::evaluate() const {
  const unsigned int bin = _observable->getBin();
  const double low = _observable->getBinning().binLow(bin);
  const double high = _observable->getBinning().binHigh(bin);

  const double oldX = _observable->getVal();
  // Important: When the integrator samples x, caching of sub-tree values needs to be off.
  RooHelpers::DisableCachingRAII disableCaching(inhibitDirty());

  auto evalFunction = [this](double x){ _observable->setVal(x); return _pdf->getVal(_normSet); };
  _integrator->SetFunction(evalFunction);

  const double result = _integrator->Integral(low, high);

  _observable->setVal(oldX);

  return result;
}


////////////////////////////////////////////////////////////////////////////////
/// Integrate the PDF over all its bins, and return a batch with those values.
/// \param[in/out] evalData Struct with evaluation data.
/// \param[in] normSet Normalisation set that's passed on to the PDF.
RooSpan<double> RooBinSamplingPdf::evaluateSpan(BatchHelpers::RunContext& evalData, const RooArgSet* normSet) const {
  // Retrieve binning, which we need to compute the probabilities
  auto boundaries = binBoundaries();
  auto xValues = _observable->getValues(evalData, normSet);
  auto results = evalData.makeBatch(this, xValues.size());

  // Important: When the integrator samples x, caching of sub-tree values needs to be off.
  RooHelpers::DisableCachingRAII disableCaching(inhibitDirty());

  auto evalFunction = [this,normSet](double x){ _observable->setVal(x); return _pdf->getVal(normSet); };
  _integrator->SetFunction(evalFunction);

  // Now approximate integral of PDF over bin by summing trapezoids. Then compute logs:
  for (unsigned int i=0; i < xValues.size(); ++i) {
    const double x = xValues[i];
    auto upperIt = std::upper_bound(boundaries.begin(), boundaries.end(), x);
    unsigned int bin = std::distance(boundaries.begin(), upperIt) - 1;
    assert(bin < boundaries.size());

    results[bin] = _integrator->Integral(boundaries[bin], boundaries[bin+1]);
  }

  return results;
}


////////////////////////////////////////////////////////////////////////////////
/// Get the bin boundaries for the observable.
/// These will be recomputed whenever the shape of this object is dirty.
RooSpan<const double> RooBinSamplingPdf::binBoundaries() const {
  if (isShapeDirty() || !_binSamplingInfo) {
    _binSamplingInfo.reset(new BinSamplingInfo);
    const RooAbsBinning& binning = _observable->getBinning(nullptr);
    const double* boundaries = binning.array();

    for (int i=0; i < binning.numBoundaries(); ++i) {
      _binSamplingInfo->binBoundaries.push_back(boundaries[i]);
    }

    assert(std::is_sorted(_binSamplingInfo->binBoundaries.begin(), _binSamplingInfo->binBoundaries.end()));

    clearShapeDirty();
  }

  return {_binSamplingInfo->binBoundaries};
}


////////////////////////////////////////////////////////////////////////////////
/// Return reference to the integrator that's used to sample the bins.
/// This can be used to alter the integration method or sampling accuracy.
ROOT::Math::IntegratorOneDim& RooBinSamplingPdf::integrator() {
  return *_integrator;
}
