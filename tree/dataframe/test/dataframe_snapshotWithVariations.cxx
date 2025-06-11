#include <ROOT/RDataFrame.hxx>

#include <TDataType.h>
#include <TFile.h>
#include <TInterpreter.h>
#include <TSystem.h>

#include <gtest/gtest.h>

using ROOT::RDF::Experimental::VariationsFor;

void printTree(TTree & tree) {
   tree.Print();
   // tree.SetScanField(40);
   tree.Scan("*", "", "colsize=15", 10);
}

template<typename X_t, typename Y_t, typename F>
void checkOutput(TTree & tree, std::vector<std::string> const & systematics, F&& activeCuts) {
   X_t x;
   Y_t y;

   for (auto const & sysName : systematics) {
      ASSERT_EQ(TTree::kMatch, tree.SetBranchAddress(("x" + sysName).c_str(), &x)) << ("x" + sysName);
      ASSERT_EQ(TTree::kMatch, tree.SetBranchAddress(("y" + sysName).c_str(), &y)) << ("y" + sysName);

      for (unsigned int i = 0; i < tree.GetEntries(); ++i) {
         ASSERT_GT(tree.GetEntry(i), 0);

         EXPECT_EQ(x, -1 * y);
         if (!activeCuts(x, y)) {
            EXPECT_EQ(x, X_t{});
            EXPECT_EQ(y, Y_t{});
         }
      }
   }
}

TEST(RDFVarySnapshot, SimpleRDFWithFilters)
{
   constexpr auto filename = "VarySnapshot.root";
   constexpr unsigned int N = 10;
   ROOT::RDF::RSnapshotOptions options;
   options.fLazy = true;
   options.fOverwriteIfExists = true;

   auto cuts = [](float x, double y) { return x < 50 || y < -70; };

   auto h = ROOT::RDataFrame(N)
               .Define("x", [](ULong64_t e) -> float { return 10.f * e; }, {"rdfentry_"})
               .Vary(
                  "x",
                  [](float x) {
                     return ROOT::RVecF{x - 0.5f, x + 0.5f};
                  },
                  {"x"}, 2, "xVar")
               .Define("y", [](float x) -> double { return -1. * x; }, {"x"})
               .Filter(cuts, {"x", "y"})
               .Snapshot<float, double>("t", filename, {"x", "y"}, options);
   auto variation = VariationsFor(h);
   variation.GetPtr();

   {
      std::unique_ptr<TFile> file{TFile::Open(filename)};
      auto tree = file->Get<TTree>("t");
      printTree(*tree);

      EXPECT_EQ(tree->GetEntries(), N);
      for (const auto branchName : {"x", "y", "x__xVar_0", "x__xVar_1", "y__xVar_0", "y__xVar_0"})
         EXPECT_NE(tree->FindBranch(branchName), nullptr);

      checkOutput<float, double>(*tree, std::vector<std::string>{"__xVar_0", "__xVar_1"}, cuts);

      if (HasFailure()) {
         tree->Print();
         tree->Scan();
      }
   }

   if (!HasFailure()) std::remove(filename);
}

TEST(RDFVarySnapshot, RDFFromTTree)
{
   constexpr auto inFile = "VarySnapshot_1.root";
   constexpr unsigned int N = 10;
   auto in = ROOT::RDataFrame(N)
      .Define("x", [](ULong64_t e) -> int { return 10 * int(e); }, {"rdfentry_"})
      .Define("y", [](int x) ->float { return -1*x; }, {"x"})
      .Snapshot<int, float>("t", inFile, {"x", "y"});
   auto nextRDF = in.GetValue();
   {
      std::unique_ptr<TFile> fileIn{TFile::Open(inFile, "READ")};
      TTree* tree = fileIn->Get<TTree>("t");
      TBranch* branch = nullptr;

      ASSERT_NE(branch = tree->GetBranch("x"), nullptr);
      EXPECT_STREQ(static_cast<TLeaf*>(branch->GetListOfLeaves()->At(0))->GetTypeName(), "Int_t");

      ASSERT_NE(branch = tree->GetBranch("y"), nullptr);
      EXPECT_STREQ(static_cast<TLeaf*>(branch->GetListOfLeaves()->At(0))->GetTypeName(), "Float_t");
      tree->Scan("*", "", "", 5);
   }

   constexpr auto filename = "VarySnapshot_2.root";
   ROOT::RDF::RSnapshotOptions options;
   options.fLazy = true;
   auto filter = [](int x, Long64_t y){ return (20 <= x && x < 70) || y > 0; };

   auto h = nextRDF
      .Vary(
         "x",
         [](int x) {
            return ROOT::RVecI{x - 1, x + 1};
         },
         {"x"}, 2, "xVariation")
      .Redefine("y", [](int x) -> Long64_t { return -1 * x; }, {"x"})
      .Filter(filter, {"x", "y"})
      .Snapshot<int, Long64_t>("t", filename, {"x", "y"}, options);
   auto variation = VariationsFor(h);
   auto thirdRDF = variation.GetValue();

   {
      std::unique_ptr<TFile> file{TFile::Open(filename)};
      auto tree = file->Get<TTree>("t");
      tree->Print();
      tree->Scan("*", "", "colsize=20", 10);

      EXPECT_EQ(tree->GetEntries(), N);
      for (const auto [branchName, branchType] : std::initializer_list<std::pair<char const*,char const*>>{{"x", "Int_t"},
            {"y", "Long64_t"},
            {"x__xVariation_0", "Int_t"},
            {"x__xVariation_1", "Int_t"},
            {"y__xVariation_0", "Long64_t"},
            {"y__xVariation_1", "Long64_t"}}) {
         TBranch *branch;
         ASSERT_NE(branch = tree->GetBranch(branchName), nullptr);
         EXPECT_STREQ(static_cast<TLeaf*>(branch->GetListOfLeaves()->At(0))->GetTypeName(), branchType);
      }

      checkOutput<int, Long64_t>(*tree, std::vector<std::string>{"__xVariation_0", "__xVariation_1"}, filter);

      if (HasFailure()) {
         tree->Print();
         tree->Scan();
      }
   }

   if (!HasFailure()) {
      std::remove(inFile);
      std::remove(filename);
   }
}

TEST(RDFVarySnapshot, WritingOfBitmask)
{
   constexpr auto filename = "VarySnapshot_bitmask.root";
   std::string const treename = "testTree";
   constexpr unsigned int N = 15;
   constexpr unsigned int NSystematics = 130; // Will use three bitmask branches
   ROOT::RDF::RSnapshotOptions options;
   options.fLazy = true;
   options.fOverwriteIfExists = true;

   auto cuts = [](int x, int y) { return x % 2 == 0 && y % 3 == 0; };

   auto h = ROOT::RDataFrame(N)
               .Define("x", [](ULong64_t e) -> int { return e; }, {"rdfentry_"})
               .Vary(
                  "x",
                  [](int x) {
                     std::vector<int> systOffsets(NSystematics);
                     std::iota(systOffsets.begin(), systOffsets.end(), 0);
                     std::transform(systOffsets.begin(), systOffsets.end(), systOffsets.begin(), [x](int offset){ return x + offset; });
                     return ROOT::RVecI{systOffsets};
                  },
                  {"x"}, NSystematics, "xVar")
               .Define("y", [](int x) -> int { return -1 * x; }, {"x"})
               .Filter(cuts, {"x", "y"})
               .Snapshot<int, int>(treename, filename, {"x", "y"}, options);
   auto variation = VariationsFor(h);
   variation.GetPtr();

   {
      TFile file(filename, "READ");
      std::unique_ptr<TTree> tree{file.Get<TTree>(treename.data())};
      ASSERT_NE(tree, nullptr);
      TBranch *branch = tree->GetBranch(("R_rdf_mask_" + treename + "_0").c_str());
      ASSERT_NE(branch, nullptr);

      auto *branchToIndexMap = file.Get<std::unordered_map<std::string, std::pair<std::string, unsigned int>>>(("R_rdf_branchToBitmaskMapping_" + treename).c_str());
      ASSERT_NE(branchToIndexMap, nullptr);
      for (const auto branchName : {"x", "y", "x__xVar_0", "x__xVar_1", "y__xVar_0", "y__xVar_0"}) {
         ASSERT_NE(branchToIndexMap->find(branchName), branchToIndexMap->end());
      }
      EXPECT_EQ((*branchToIndexMap)["x"], (*branchToIndexMap)["y"]);

      for (unsigned int systematic = 0; systematic < NSystematics; ++systematic) {
         tree->ResetBranchAddresses();
         int x,y;
         uint64_t bitmask;
         std::string systematicName = std::string{"__xVar_"} + std::to_string(systematic);
         ASSERT_EQ(tree->SetBranchAddress(("x" + systematicName).c_str(), &x), TTree::kMatch);
         ASSERT_EQ(tree->SetBranchAddress(("y" + systematicName).c_str(), &y), TTree::kMatch);
         const auto it = branchToIndexMap->find("x" + systematicName);
         ASSERT_NE(it, branchToIndexMap->end());
         ASSERT_EQ(tree->SetBranchAddress(it->second.first.c_str(), &bitmask), TTree::kMatch) << it->second.first;
         for (unsigned int i = 0; i < N; ++i) {
            const int xOrig = i + systematic;
            const int yOrig = -1 * xOrig;
            tree->GetEntry(i);
            if (cuts(xOrig, yOrig)) {
               ASSERT_EQ(x, xOrig) << "event=" << i << " systematic=" << systematic;
               ASSERT_EQ(y, yOrig) << "event=" << i << " systematic=" << systematic;
            }

            const std::bitset<64> bs{bitmask};
            const auto bitIndex = it->second.second;
            EXPECT_EQ(cuts(xOrig, yOrig), bs[bitIndex]) << "event=" << i << " syst=" << systematic << " x=" << xOrig << " y=" << yOrig << " bitset: " << bs << " bitIndex: " << bitIndex;
         }
      }

      if (HasFailure()) {
         tree->Print();
         tree->Scan();
      }
   }

   if (!HasFailure()) std::remove(filename);
}

TEST(DISABLED_RDFVarySnapshot, SnapshotCollections)
{
   constexpr auto filename = "VarySnapshot.root";
   constexpr unsigned int N = 10;
   ROOT::RDF::RSnapshotOptions options;
   options.fLazy = true;

   gInterpreter->GenerateDictionary("std::vector<ROOT::VecOps::RVec<int> >", "vector;ROOT/RVec.hxx");

   auto cuts = [](int x, ROOT::RVecI const & y) { return x % 2 == 0 && y.size() < 18; };

   auto h = ROOT::RDataFrame(N)
            .Define("x", [](ULong64_t e) -> int { return int(e); }, {"rdfentry_"})
            .Vary(
               "x",
               [](int x) {
                  return ROOT::RVecI{x - 1, x + 1};
               },
               {"x"}, 2, "xVariation")
            .Define("y", [](int x) { return ROOT::RVecI{ x, x+1, x+2 }; }, {"x"})
            .Filter(cuts, {"x", "y"})
            .Snapshot<int, ROOT::RVecI>("t", filename, {"x", "y"}, options);
   auto variation = VariationsFor(h);
   variation.GetPtr();

   std::unique_ptr<TFile> file{TFile::Open(filename)};
   auto tree = file->Get<TTree>("t");
   printTree(*tree);

   EXPECT_EQ(tree->GetEntries(), N);
   for (const auto branchName : {"x", "xVariation:0:x", "xVariation:1:x", "xVariation:0:y", "xVariation:1:y"})
      EXPECT_NE(tree->GetBranch(branchName), nullptr) << branchName;

   for (std::string sysPrefix : {"", "xVariation:0:", "xVariation:1:"}) {
      std::vector<int> x;
      std::vector<ROOT::RVecI> y;
      TBranch *xBranch, *yBranch;
      ASSERT_NE(xBranch = tree->Branch((sysPrefix + "x").c_str(), &x), nullptr) << (sysPrefix + "x");
      ASSERT_NE(yBranch = tree->Branch((sysPrefix + "y").c_str(), &y), nullptr) << (sysPrefix + "y");
      // checkOutput<int, ROOT::RVecI>(*tree, std::vector<std::string>{"xVariation:0:", "xVariation:1:"}, cuts);
   }

   file.reset();
   gSystem->Unlink(filename);
}