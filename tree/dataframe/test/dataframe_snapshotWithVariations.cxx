#include <ROOT/RDataFrame.hxx>

#include <TDataType.h>
#include <TFile.h>
#include <TInterpreter.h>
#include <TSystem.h>

#include <gtest/gtest.h>

using ROOT::RDF::Experimental::VariationsFor;

void printTree(TTree & tree) {
   tree.Print();
   tree.Scan("*", "", "colsize=20", 10);
}

template<typename X_t, typename Y_t, typename F>
void checkOutput(TTree & tree, std::vector<std::string> const & systematics, F&& activeCuts) {
   std::vector<X_t> *x = nullptr;
   std::vector<Y_t> y;

   for (auto const & sysPrefix : systematics) {
      TBranch *xBranch, *yBranch;
      ASSERT_NE(xBranch = tree.Branch((sysPrefix + "x").c_str(), &x), nullptr) << (sysPrefix + "x");
      ASSERT_NE(yBranch = tree.Branch((sysPrefix + "y").c_str(), &y), nullptr) << (sysPrefix + "y");

      for (unsigned int i = 0; i < tree.GetEntries(); ++i) {
         ASSERT_GT(tree.GetEntry(i), 0);
         if (!x->empty() && !y.empty()) {
            EXPECT_EQ(y.front(), 2.f * x->front());
            EXPECT_TRUE(activeCuts(x->front(), y.front())) << x->front() << " " << y.front();
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

   auto cuts = [](int x, float y) { return x % 2 == 0 && y < 18; };

   auto h = ROOT::RDataFrame(N)
               .Define("x", [](ULong64_t e) -> int { return int(e); }, {"rdfentry_"})
               .Vary(
                  "x",
                  [](int x) {
                     return ROOT::RVecI{x - 1, x + 1};
                  },
                  {"x"}, 2, "xVariation")
               .Define("y", [](int x) -> float { return 2.f*x; }, {"x"})
               .Filter(cuts, {"x", "y"})
               .Snapshot<int, float>("t", filename, {"x", "y"}, options);
   auto variation = VariationsFor(h);
   variation.GetPtr();

   std::unique_ptr<TFile> file{TFile::Open(filename)};
   auto tree = file->Get<TTree>("t");
   printTree(*tree);

   EXPECT_EQ(tree->GetEntries(), N);
   for (const auto branchName : {"x", "xVariation:0:x", "xVariation:1:x", "xVariation:0:y", "xVariation:1:y"})
      EXPECT_NE(tree->GetBranch(branchName), nullptr);

   checkOutput<int, float>(*tree, std::vector<std::string>{"xVariation:0:", "xVariation:1:"}, cuts);

   file.reset();
   gSystem->Unlink(filename);
}

TEST(RDFVarySnapshot, RDFFromTTree)
{
   constexpr auto inFile = "VarySnapshot_in.root";
   constexpr unsigned int N = 10;
   auto in = ROOT::RDataFrame(N)
      .Define("x", [](ULong64_t e) -> int { return int(e); }, {"rdfentry_"})
      .Define("y", [](int x) ->float { return 2.f*x; }, {"x"})
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
   }


   constexpr auto filename = "VarySnapshot.root";
   ROOT::RDF::RSnapshotOptions options;
   options.fLazy = true;

   auto h = nextRDF
      .Vary(
         "x",
         [](int x) {
            return ROOT::RVecI{x - 1, x + 1};
         },
         {"x"}, 2, "xVariation")
      .Redefine("y", [](int x){ return 2.f * x; }, {"x"})
      .Snapshot<int, float>("t", filename, {"x", "y"}, options);
   auto variation = VariationsFor(h);
   variation.GetValue();

   std::unique_ptr<TFile> file{TFile::Open(filename)};
   auto tree = file->Get<TTree>("t");
   tree->Print();
   tree->Scan("*", "", "colsize=20", 10);

   EXPECT_EQ(tree->GetEntries(), N);
   for (const auto branchName : {"x", "xVariation:0:x", "xVariation:1:x", "xVariation:0:y", "xVariation:1:y"})
      EXPECT_NE(tree->GetBranch(branchName), nullptr);
   checkOutput<int, float>(*tree, std::vector<std::string>{"xVariation:0:", "xVariation:1:"}, [](int, float){ return true; });

   file.reset();
   gSystem->Unlink(filename);
   gSystem->Unlink(inFile);
}

TEST(RDFVarySnapshot, SnapshotCollections)
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
               .Define("y", [](int x) { return ROOT::RVecI{x, x+1, x+2 }; }, {"x"})
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