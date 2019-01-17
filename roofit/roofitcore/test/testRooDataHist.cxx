// Tests for the RooWorkspace
// Authors: Stephan Hageboeck, CERN  01/2019

#include "RooDataHist.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooHelpers.h"
#include "TH1D.h"

#include "gtest/gtest.h"

/// ROOT-8163
/// The RooDataHist warns that it has to adjust the binning of x to the next bin boundary
/// although the boundaries match perfectly.
TEST(RooDataHist, BinningRangeCheck_8163)
{
  RooHelpers::HijackMessageStream hijack(RooFit::INFO, RooFit::DataHandling, "dataHist");

  RooRealVar x("x", "x", 0., 1.);
  TH1D hist("hist", "", 10, 0., 1.);

  RooDataHist dataHist("dataHist", "", RooArgList(x), &hist);
  EXPECT_TRUE(hijack.str().empty()) << "Messages issued were: " << hijack.str();
}

TEST(RooDataHist, assign_8522) {
   RooDataHist dhist;
   TH1D hist("hist", "Test histogram", 100, -5.0, 5.0);
   hist.FillRandom("gaus", 1000);

   // Define a dummy RooRealVar:
   RooRealVar x("x", "x", -5.0, 5.0);
   dhist = RooDataHist("test", "TestDHist", x, &hist);
   EXPECT_EQ(dhist.numEntries(), 1000);
}
