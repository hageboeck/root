/*
 * testNLLVar.cxx
 *
 *  Created on: 8 Oct 2020
 *      Author: Stephan Hageboeck
 */

#include <RooRealVar.h>
#include <RooGenericPdf.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooFitResult.h>

#include <gtest/gtest.h>

#include <memory>

TEST(RooNLLVar, HighResolutionSampling) {
  RooRealVar x("x", "x", 0.1, 5.);
  x.setBins(10);

  RooRealVar a("a", "a", -0.3, -5., 5.);
  RooArgSet targetValues;
  RooArgSet(a).snapshot(targetValues);

  RooGenericPdf pdf("gaussian", "std::pow(x, a)", RooArgSet(x, a));
  std::unique_ptr<RooDataHist> dataH(pdf.generateBinned(x,  10000));
  RooRealVar w("w", "weight", 0., 0., 10000.);
  RooDataSet data("data", "data", RooArgSet(x, w), RooFit::WeightVar(w));
  for (int i=0; i < dataH->numEntries(); ++i) {
    data.add(*dataH->get(i), dataH->weight());
  }

  std::unique_ptr<RooPlot> frame( x.frame() );
  dataH->plotOn(frame.get(), RooFit::MarkerColor(kRed));
  data.plotOn(frame.get(), RooFit::Name("data"));


  a.setVal(3.);
  std::unique_ptr<RooFitResult> fit1( pdf.fitTo(data, RooFit::Save(), RooFit::PrintLevel(-1)) );
  pdf.plotOn(frame.get(), RooFit::LineColor(kRed), RooFit::Name("standard"));

  a.setVal(3.);
  std::unique_ptr<RooFitResult> fit2( pdf.fitTo(data, RooFit::Save(), RooFit::PrintLevel(-1),
      RooFit::BatchMode(true), RooFit::HighResolutionSampling(50)) );
  pdf.plotOn(frame.get(), RooFit::LineColor(kBlue), RooFit::Name("highRes"));


  auto getVal = [](const char* name, const RooArgSet& set) {
    return dynamic_cast<const RooRealVar&>(set[name]).getVal();
  };
  auto getErr = [](const char* name, const RooArgSet& set) {
    return dynamic_cast<const RooRealVar&>(set[name]).getError();
  };

  EXPECT_GT(fabs(getVal("a", targetValues) - getVal("a", fit1->floatParsFinal())), 1. * getErr("a", fit1->floatParsFinal()))
      << "Expecting a bias when sampling PDF in bin centre.";

  EXPECT_NEAR(getVal("a", targetValues), getVal("a", fit2->floatParsFinal()), 1. * getErr("a", fit2->floatParsFinal()))
      << "Expect reduced bias with high-resolution sampling.";

  EXPECT_GT(frame->chiSquare("standard", "data", 1) * 0.9, frame->chiSquare("highRes",  "data", 1))
      << "Expect chi2/ndf at least 10% better.";

//  fit1->Print();
//  fit2->Print();
}
