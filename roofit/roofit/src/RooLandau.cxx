/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/** \class RooLandau
    \ingroup Roofit

Landau distribution p.d.f
\image html RF_Landau.png "PDF of the Landau distribution."
**/

#include "RooLandau.h"
#include "TMath.h"
#include "RooFit.h"

#include "RooRandom.h"
#include "BatchHelpers.h"

#include "TError.h"

#include "LandauBatchEvaluate.h" // for compute()

using namespace std;
using namespace BatchHelpers;
using namespace LandauBatchEvaluate;

ClassImp(RooLandau);

////////////////////////////////////////////////////////////////////////////////

RooLandau::RooLandau(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _mean, RooAbsReal& _sigma) :
  RooAbsPdf(name,title),
  x("x","Dependent",this,_x),
  mean("mean","Mean",this,_mean),
  sigma("sigma","Width",this,_sigma)
{
}

////////////////////////////////////////////////////////////////////////////////

RooLandau::RooLandau(const RooLandau& other, const char* name) :
  RooAbsPdf(other,name),
  x("x",this,other.x),
  mean("mean",this,other.mean),
  sigma("sigma",this,other.sigma)
{
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooLandau::evaluate() const
{
  return TMath::Landau(x, mean, sigma);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute \f$ Landau(x,mean,sigma) \f$ in batches.
/// The local proxies {x, mean, sigma} will be searched for batch input data,
/// and if found, the computation will be batched over their
/// values. If batch data are not found for one of the proxies, the proxies value is assumed to
/// be constant over the batch.
/// \param[in] batchIndex Index of the batch to be computed.
/// \param[in] batchSize Size of each batch. The last batch may be smaller.
/// \return A span with the computed values.

RooSpan<double> RooLandau::evaluateBatch(std::size_t begin, std::size_t batchSize) const {
  auto output = _batchData.makeWritableBatchUnInit(begin, batchSize);

  auto xData = x.getValBatch(begin, batchSize);
  auto meanData = mean.getValBatch(begin, batchSize);
  auto sigmaData = sigma.getValBatch(begin, batchSize);
  const bool batchX = !xData.empty();
  const bool batchMean = !meanData.empty();
  const bool batchSigma = !sigmaData.empty();
  if (batchX && !batchMean && !batchSigma ) {
    compute(output, xData, BracketAdapter<RooRealProxy>(mean), BracketAdapter<RooRealProxy>(sigma));
  }
  else if (!batchX && batchMean && !batchSigma ) {
    compute(output, BracketAdapter<RooRealProxy>(x), meanData, BracketAdapter<RooRealProxy>(sigma));
  }
  else if (batchX && batchMean && !batchSigma ) {
    compute(output, xData, meanData, BracketAdapter<RooRealProxy>(sigma));
  }
  else if (!batchX && !batchMean && batchSigma ) {
    compute(output, BracketAdapter<RooRealProxy>(x), BracketAdapter<RooRealProxy>(mean), sigmaData);
  }
  else if (batchX && !batchMean && batchSigma ) {
    compute(output, xData, BracketAdapter<RooRealProxy>(mean), sigmaData);
  }
  else if (!batchX && batchMean && batchSigma ) {
    compute(output, BracketAdapter<RooRealProxy>(x), meanData, sigmaData);
  }
  else if (batchX && batchMean && batchSigma ) {
    compute(output, xData, meanData, sigmaData);
  }
  else{
    throw std::logic_error("Requested a batch computation, but no batch data available.");
  }

  return output;
}



////////////////////////////////////////////////////////////////////////////////

Int_t RooLandau::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t /*staticInitOK*/) const
{
  if (matchArgs(directVars,generateVars,x)) return 1 ;
  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

void RooLandau::generateEvent(Int_t code)
{
  R__ASSERT(code==1) ;
  Double_t xgen ;
  while(1) {
    xgen = RooRandom::randomGenerator()->Landau(mean,sigma);
    if (xgen<x.max() && xgen>x.min()) {
      x = xgen ;
      break;
    }
  }
  return;
}
