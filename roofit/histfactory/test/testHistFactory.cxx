// Tests for the HistFactory
// Authors: Stephan Hageboeck, CERN  01/2019

#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "RooStats/HistFactory/Sample.h"
#include "RooStats/ModelConfig.h"

#include "RooWorkspace.h"
#include "RooArgSet.h"
#include "RooSimultaneous.h"
#include "RooRealSumPdf.h"
#include "RooRealVar.h"
#include "RooHelpers.h"
#include "RooFitResult.h"
#include "RooPlot.h"



#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "gtest/gtest.h"

using namespace RooStats;
using namespace RooStats::HistFactory;

TEST(Sample, CopyAssignment)
{
  Sample s("s");
  {
    Sample s1("s1");
    auto hist1 = new TH1D("hist1", "hist1", 10, 0, 10);
    s1.SetHisto(hist1);
    s = s1;
    //Now go out of scope. Should delete hist1, that's owned by s1.
  }
  
  auto hist = s.GetHisto();
  ASSERT_EQ(hist->GetNbinsX(), 10);
}


TEST(HistFactory, Read_ROOT6_16_Model) {
  std::string filename = "./ref_6.16_example_UsingC_channel1_meas_model.root";
  std::unique_ptr<TFile> file(TFile::Open(filename.c_str()));
  if (!file || !file->IsOpen()) {
    filename = TROOT::GetRootSys() + "/roofit/histfactory/test/" + filename;
    file.reset(TFile::Open(filename.c_str()));
  }

  ASSERT_TRUE(file && file->IsOpen());
  RooWorkspace* ws;
  file->GetObject("channel1", ws);
  ASSERT_NE(ws, nullptr);

  auto mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
  ASSERT_NE(mc, nullptr);

  RooAbsPdf* pdf = mc->GetPdf();
  ASSERT_NE(pdf, nullptr);

  const RooArgSet* obs = mc->GetObservables();
  ASSERT_NE(obs, nullptr);

  EXPECT_NEAR(pdf->getVal(), 0.17488817, 1.E-8);
  EXPECT_NEAR(pdf->getVal(*obs), 0.95652174, 1.E-8);
}


TEST(HistFactory, Read_ROOT6_16_Combined_Model) {
  std::string filename = "./ref_6.16_example_UsingC_combined_meas_model.root";
  std::unique_ptr<TFile> file(TFile::Open(filename.c_str()));
  if (!file || !file->IsOpen()) {
    filename = TROOT::GetRootSys() + "/roofit/histfactory/test/" + filename;
    file.reset(TFile::Open(filename.c_str()));
  }

  ASSERT_TRUE(file && file->IsOpen());
  RooWorkspace* ws;
  file->GetObject("combined", ws);
  ASSERT_NE(ws, nullptr);

  auto mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
  ASSERT_NE(mc, nullptr);

  RooAbsPdf* pdf = mc->GetPdf();
  ASSERT_NE(pdf, nullptr);

  const RooArgSet* obs = mc->GetObservables();
  ASSERT_NE(obs, nullptr);

  EXPECT_NEAR(pdf->getVal(), 0.17488817, 1.E-8);
  EXPECT_NEAR(pdf->getVal(*obs), 0.95652174, 1.E-8);
}


enum MakeModelMode {kEquidistantBinning, kCustomBinning};
class MakeModel : public testing::TestWithParam<MakeModelMode> {
public:
  void SetUp() {
    _inputFile = "TestMakeModel.root";

    TFile example(_inputFile.c_str(), "RECREATE");
    TH1F *data, *signal, *bkg1, *bkg2, *statUnc = nullptr;
    if (GetParam() == kEquidistantBinning) {
      data = new TH1F("data","data", 2,1,2);
      signal = new TH1F("signal","signal histogram (pb)", 2,1,2);
      bkg1 = new TH1F("background1","background 1 histogram (pb)", 2,1,2);
      bkg2 = new TH1F("background2","background 2 histogram (pb)", 2,1,2);
      statUnc = new TH1F("background1_statUncert", "statUncert", 2,1,2);
    } else if (GetParam() == kCustomBinning) {
      double binning[3] = {0., 1.8, 2.};
      data = new TH1F("data","data", 3, binning);
      signal = new TH1F("signal","signal histogram (pb)", 2, binning);
      bkg1 = new TH1F("background1","background 1 histogram (pb)", 2, binning);
      bkg2 = new TH1F("background2","background 2 histogram (pb)", 2, binning);
      statUnc = new TH1F("background1_statUncert", "statUncert", 2, binning);
    } else {
      throw std::logic_error("Test suite set up wrongly");
    }

    data->SetBinContent(1, 140);
    data->SetBinContent(2, 120);

    signal->SetBinContent(1,20);
    signal->SetBinContent(2,10);

    bkg1->SetBinContent(1,100);
    bkg2->SetBinContent(2,100);

    // A small statistical uncertainty
    statUnc->SetBinContent(1, .05);  // 5% uncertainty
    statUnc->SetBinContent(2, .05);  // 5% uncertainty

    for (auto hist : {data, signal, bkg1, bkg2, statUnc}) {
      example.WriteTObject(hist);
    }
  }

  void TearDown() {

  }

  std::string _inputFile;
};

TEST_P(MakeModel, MakingModels) {
  // Create the measurement
  Measurement meas("meas", "meas");

  meas.SetOutputFilePrefix( "example_variableBins" );
  meas.SetPOI( "SigXsecOverSM" );
  meas.AddConstantParam("alpha_syst1");
  meas.AddConstantParam("Lumi");
  meas.SetExportOnly(true);

  meas.SetLumi( 1.0 );
  meas.SetLumiRelErr( 0.10 );


  // Create a channel
  Channel chan( "channel1" );
  chan.SetData( "data", _inputFile);
  chan.SetStatErrorConfig( 0.05, "Poisson" );


  // Now, create some samples

  // Create the signal sample
  Sample signal( "signal", "signal", _inputFile);
  signal.AddOverallSys( "syst1",  0.95, 1.05 );
  signal.AddNormFactor( "SigXsecOverSM", 1, 0, 3 );
  chan.AddSample( signal );

  // Background 1
  Sample background1( "background1", "background1", _inputFile);
  background1.ActivateStatError( "background1_statUncert", _inputFile);
  background1.AddOverallSys( "syst2", 0.95, 1.05  );
  chan.AddSample( background1 );

  // Background 1
  Sample background2( "background2", "background2", _inputFile);
  background2.ActivateStatError();
  background2.AddOverallSys( "syst3", 0.95, 1.05  );
  chan.AddSample( background2 );


  // Done with this channel
  // Add it to the measurement:
  meas.AddChannel( chan );

  RooHelpers::HijackMessageStream hijackW(RooFit::WARNING, RooFit::HistFactory);
  RooHelpers::HijackMessageStream hijack(RooFit::INFO, RooFit::HistFactory);

  // Collect the histograms from their files,
  meas.CollectHistograms();

  // Now, create the measurement
  RooWorkspace* ws = MakeModelAndMeasurementFast( meas );

  EXPECT_TRUE(hijackW.str().empty()) << "Warning stream for HistFactory is " << hijackW.str();

  auto simPdf = dynamic_cast<RooSimultaneous*>(ws->pdf("simPdf"));
  EXPECT_NE(simPdf, nullptr);

  auto channelPdf = dynamic_cast<RooRealSumPdf*>(ws->pdf("channel1_model"));
  EXPECT_NE(channelPdf, nullptr);
//  channelPdf->Print("T");
  channelPdf->graphVizTree("/tmp/graphVizTree.dot");

  auto obs = dynamic_cast<RooRealVar*>(ws->var("obs_x_channel1"));
  ASSERT_NE(obs, nullptr);
  if (GetParam() == kEquidistantBinning) {
    EXPECT_DOUBLE_EQ(obs->getBinWidth(0), 0.5);
    EXPECT_DOUBLE_EQ(obs->getBinWidth(1), 0.5);
    EXPECT_EQ(obs->numBins(), 2);
  } else {
    EXPECT_DOUBLE_EQ(obs->getBinWidth(0), 0.8);
    EXPECT_DOUBLE_EQ(obs->getBinWidth(1), 0.2);
    EXPECT_EQ(obs->numBins(), 2);
  }

  RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
  ASSERT_NE(mc, nullptr);

  for (const auto& systName : std::initializer_list<std::string>{"alpha_syst1", "alpha_syst2", "alpha_syst3"}) {
    auto& var = *ws->var(systName.c_str());

    EXPECT_TRUE(channelPdf->dependsOnValue(var)) << "Expect channel pdf to depend on " << systName;
    if (!var.isConstant())
      EXPECT_NE(mc->GetNuisanceParameters()->find(systName.c_str()), nullptr) << systName << " should be in list of nuisance parameters.";
  }

  EXPECT_EQ(*mc->GetParametersOfInterest()->begin(), ws->var("SigXsecOverSM"));

  obs->setBin(0);
  EXPECT_NEAR(channelPdf->getVal(/*mc->GetObservables()*/), 120, 1.E-4);
  obs->setBin(1);
  EXPECT_NEAR(channelPdf->getVal(mc->GetObservables()), 110, 1.E-4);

  RooAbsData* data = dynamic_cast<RooAbsData*>(ws->data("obsData"));
  ASSERT_NE(data, nullptr);
  data->Print("V");
  std::cout << data->IsA()->GetName() << std::endl;

  auto fitResult = simPdf->fitTo(*data, RooFit::GlobalObservables(*mc->GetGlobalObservables()), RooFit::Save(), RooFit::PrintLevel(0));
  ASSERT_NE(fitResult, nullptr);
  EXPECT_EQ(fitResult->status(), 0);

  // Model is set up such that both backgrounds should be close to 1, and signal = 2
  EXPECT_NEAR(dynamic_cast<RooRealVar*>(fitResult->floatParsFinal().find("SigXsecOverSM"))->getVal(), 2., 0.001);
  EXPECT_NEAR(dynamic_cast<RooRealVar*>(fitResult->floatParsFinal().find("alpha_syst2"))->getVal(), 0., 0.001);
  EXPECT_NEAR(dynamic_cast<RooRealVar*>(fitResult->floatParsFinal().find("alpha_syst3"))->getVal(), 0., 0.001);

  auto frame = obs->frame();
  data->plotOn(frame);
  channelPdf->plotOn(frame);
  channelPdf->plotOn(frame, RooFit::Components("signal_channel1_shapes"), RooFit::LineColor(kRed));
  TCanvas canv;
  frame->Draw();
  canv.Draw();
  canv.SaveAs("/tmp/HFTest.png");
}


INSTANTIATE_TEST_SUITE_P(HistFactory,
    MakeModel,
    testing::Values(kEquidistantBinning, kCustomBinning));

