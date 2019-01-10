// Tests for the RooWorkspace
// Authors: Stephan Hageboeck, CERN  01/2019

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooGlobalFunc.h"
#include "RooGaussian.h"
#include "RooStats/ModelConfig.h"

#include "TFile.h"
#include "TSystem.h"

#include "gtest/gtest.h"

class TestRooWorkspaceWithGaussian : public ::testing::Test {
protected:
  TestRooWorkspaceWithGaussian() :
  Test()
  {
    RooRealVar x("x", "x", 1, 0, 10);
    RooRealVar mu("mu", "mu", 1, 0, 10);
    RooRealVar sigma("sigma", "sigma", 1, 0, 10);

    RooGaussian pdf("Gauss", "Gauss", x, mu, sigma);

    TFile outfile(_filename, "RECREATE");

    // now create the model config for this problem
    RooWorkspace w("ws");
    RooStats::ModelConfig modelConfig("ModelConfig", &w);
    modelConfig.SetPdf(pdf);
    modelConfig.SetParametersOfInterest(RooArgSet(sigma));
    modelConfig.SetGlobalObservables(RooArgSet(mu));
    w.import(modelConfig);

    outfile.WriteObject(&w, "ws");
  }

  ~TestRooWorkspaceWithGaussian() {
    gSystem->Unlink(_filename);
  }

  const char* _filename = "ROOT-9777.root";
};


/// ROOT-9777, cloning a RooWorkspace. The ModelConfig did not get updated
/// when a workspace was cloned, and was hence pointing to a non-existing workspace.
///
TEST_F(TestRooWorkspaceWithGaussian, CloneModelConfig_ROOT_9777)
{
  std::unique_ptr<RooWorkspace> w2;
  {
    TFile infile(_filename, "READ");
    RooWorkspace *w;
    infile.GetObject("ws", w);
    ASSERT_TRUE(w) << "Workspace not read from file.";

    w2 = std::make_unique<RooWorkspace>(*w);
    delete w;
  }

  RooStats::ModelConfig *mc = dynamic_cast<RooStats::ModelConfig*>(w2->genobj("ModelConfig"));

  EXPECT_TRUE(mc) << "ModelConfig not retrieved.";
  EXPECT_TRUE(mc->GetGlobalObservables()) << "GlobalObservables in ModelConfig broken.";
  EXPECT_TRUE(mc->GetParametersOfInterest()) << "ParametersOfInterest in ModelConfig broken.";
}

