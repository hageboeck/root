/**
 \file ROOT/RDF/Snapshot.hxx
 \ingroup dataframe
 \author Stephan Hageboeck, CERN
 \date 2025-02
*/

/*************************************************************************
 * Copyright (C) 1995-2025, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RDFSNAPSHOT
#define ROOT_RDFSNAPSHOT

#include "ROOT/RDF/RAction.hxx"
#include "ROOT/RDF/RActionImpl.hxx"
#include "ROOT/RDF/Utils.hxx"
#include "ROOT/RSnapshotOptions.hxx"

#include "ROOT/RDF/ActionHelpers.hxx" // TODO: Get rid of this include?

#include "TError.h"
#include "TFile.h"
#include "TClass.h"

#include <algorithm>
#include <memory>

namespace ROOT::Internal::RDF {
   using namespace ROOT::TypeTraits;
   using namespace ROOT::RDF;

/// Store file and tree in one common place to share them between instances.
struct FileHandle {
   std::unique_ptr<TFile> fFile;
   std::unique_ptr<TTree> fTree;
   std::string fDirectoryName;
   RLoopManager * fOutputLoopManager;

   FileHandle(TFile * file, TTree * tree) : fFile{file}, fTree{tree} { }
   ~FileHandle() {
      // use AutoSave to flush TTree contents because TTree::Write writes in gDirectory, not in fDirectory
      if (fTree) {
         fTree->AutoSave("flushbaskets");

         // Now connect the data source to the loop manager so it can be used for further processing
         std::string tree = fTree->GetName();
         if (!fDirectoryName.empty()) tree = fDirectoryName + '/' + tree;
         std::string file = fFile->GetName();

         fTree.reset();
         fFile.reset();

         if (fOutputLoopManager)
            fOutputLoopManager->SetDataSource(std::make_unique<ROOT::Internal::RDF::RTTreeDS>(tree, file));
      }
   }
};

struct BranchData {
   TBranch * fOutputBranch = nullptr;
   std::function<void(void*)> fDeleterOfEmptyInstance;
   void * fEmptyInstance;
   TBranch * fInputBranch = nullptr;

   BranchData(TBranch * branch, void* emptyInstance, std::function<void(void*)> deleter) :
      fOutputBranch{branch},
      fDeleterOfEmptyInstance{std::move(deleter)},
      fEmptyInstance{emptyInstance}
   {

   }
   BranchData(const BranchData&) = delete;
   BranchData(BranchData && other)
   {
      *this = std::move(other);
   }
   ~BranchData() {
      if (fEmptyInstance && fDeleterOfEmptyInstance)
         fDeleterOfEmptyInstance(fEmptyInstance);
   }
   BranchData& operator=(const BranchData&) = delete;
   BranchData& operator=(BranchData && other) {
      fOutputBranch = other.fOutputBranch;
      fDeleterOfEmptyInstance = std::move(other.fDeleterOfEmptyInstance);
      fEmptyInstance = other.fEmptyInstance;
      fInputBranch = other.fInputBranch;
      other.fEmptyInstance = nullptr;
      return *this;
   }

   /// Point the branch address to an empty instance of the type represented by this branch.
   /// This is used in case of variations, when certain defines/actions don't execute. We
   /// nevertheless need to write something, so we point the branch to an empty instance.
   void ResetBranchAddressToEmtpyInstance() {
      assert(fEmptyInstance);
      fOutputBranch->SetAddress(fEmptyInstance);
   }
};

template <typename T>
void SetBranchesHelper(TTree *inputTree, TTree &outputTree, const std::string &inName, const std::string &name,
      T *address, BranchData & bd, bool /*isDefine*/, int basketSize)
{
   static const TClassRef TBOClRef("TBranchObject");

   if (!bd.fOutputBranch) {
      if (!bd.fInputBranch && inputTree) {
         bd.fInputBranch = inputTree->GetBranch(inName.c_str());
         if (!bd.fInputBranch) // try harder
            bd.fInputBranch = inputTree->FindBranch(inName.c_str());
      }

      if (bd.fInputBranch) {
         // Respect the original bufsize and splitlevel arguments
         // In particular, by keeping splitlevel equal to 0 if this was the case for `inputBranch`, we avoid
         // writing garbage when unsplit objects cannot be written as split objects (e.g. in case of a polymorphic
         // TObject branch, see https://bit.ly/2EjLMId ).
         // A user-provided basket size value takes precedence.
         const auto bufSize = (basketSize > 0) ? basketSize : bd.fInputBranch->GetBasketSize();
         const auto splitLevel = bd.fInputBranch->GetSplitLevel();

         if (bd.fInputBranch->IsA() == TBOClRef) {
            // Need to pass a pointer to pointer
            bd.fOutputBranch =
               outputTree.Branch(name.c_str(), reinterpret_cast<T **>(bd.fInputBranch->GetAddress()), bufSize, splitLevel);
         } else {
            bd.fOutputBranch = outputTree.Branch(name.c_str(), address, bufSize, splitLevel);
         }
      } else {
         // Set Custom basket size for new branches.
         const auto buffSize = (basketSize > 0) ? basketSize : 32000;
         bd.fOutputBranch = outputTree.Branch(name.c_str(), address, buffSize);
      }

      // Create an empty instance of this type. This will be written to the tree if a systematic
      // uncertainty didn't pass the cuts, but another did.
      auto const dict = TDictionary::GetDictionary(typeid(T));
      if (auto dataType = dynamic_cast<TDataType const*>(dict); dataType) {
         bd.fEmptyInstance = new T{}; // replace by TDataType if the template goes away
         bd.fDeleterOfEmptyInstance = [](void * ptr){ std::default_delete<T>{}(static_cast<T*>(ptr)); };
      } else if (auto tClass = dynamic_cast<TClass const*>(dict); tClass) {
         bd.fEmptyInstance = tClass->New();
         bd.fDeleterOfEmptyInstance = [tClass](void* ptr) {
            tClass->GetDestructor()(ptr);
            tClass->GetDelete()(ptr);
         };
      }
   }

   // the output branch was already created, we just need to (re)set its address
   if (bd.fInputBranch && bd.fInputBranch->IsA() == TBOClRef) {
      bd.fOutputBranch->SetAddress(reinterpret_cast<T **>(bd.fInputBranch->GetAddress()));
   } else if (bd.fOutputBranch->IsA() != TBranch::Class()) {
      bd.fOutputBranch->SetAddress(&address);
   } else {
      bd.fOutputBranch->SetAddress(address);
   }
}

/// Helper object for a single-thread Snapshot action
template <typename... ColTypes>
class R__CLING_PTRCHECK(off) SnapshotHelperWithVariations : public ROOT::Detail::RDF::RActionImpl<SnapshotHelperWithVariations<ColTypes...>> {
   RSnapshotOptions fOptions;
   std::shared_ptr<FileHandle> fOutputHandle;
   // ColumnNames_t fInputBranchNames; // This contains the resolved aliases
   ColumnNames_t fOutputBranchNames;
   std::unique_ptr<std::vector<BranchData>> fBranchData;
   std::vector<std::vector<BranchData>*> fBranchDataToClear;
   ROOT::Detail::RDF::RLoopManager *fInputLoopManager = nullptr;
   ROOT::Detail::RDF::RLoopManager *fOutputLoopManager = nullptr;


   template<typename T>
   void CreateOutputBranch(std::string const & name)
   {
      // TODO: Review
      std::string name_noColon{name};
      std::replace(name_noColon.begin(), name_noColon.end(), ':', '_');

      fBranchData->emplace_back(nullptr, nullptr, nullptr);
      // TODO: Review the dummy arguments
      SetBranchesHelper(nullptr, *fOutputHandle->fTree, "", name_noColon, static_cast<T*>(nullptr), fBranchData->back(), /*isDefine*/false, fOptions.fBasketSize);
   }

   template <std::size_t... S>
   void CreateOutputBranches(std::vector<std::string> const & branchNames, std::index_sequence<S...>)
   {
      // TODO: Check basket size of equivalent input branches?
      (CreateOutputBranch<ColTypes>(branchNames[S]), ...);
   }

   template <std::size_t... S>
   void SetBranches(ColTypes &... values, std::index_sequence<S...> /*dummy*/)
   {
      // TODO: Review empty arguments
      (SetBranchesHelper(nullptr,
            *fOutputHandle->fTree, "", fOutputBranchNames[S], &values, (*fBranchData)[S], /*fIsDefine[S]*/false, fOptions.fBasketSize),
       ...);
      // TODO: Review:
      // fOutputBranches.AssertNoNullBranchAddresses();
   }

   SnapshotHelperWithVariations(SnapshotHelperWithVariations & other, std::string const & variationPrefix)
   : fOptions{other.fOptions}, fOutputHandle{other.fOutputHandle}, fBranchData{new std::vector<BranchData>{}}
   {
      fOutputBranchNames.resize(other.fOutputBranchNames.size());
      std::transform(other.fOutputBranchNames.begin(), other.fOutputBranchNames.end(), fOutputBranchNames.begin(),
      [&variationPrefix](std::string const & originalName){ return originalName + variationPrefix;});

      fOutputBranchNames = ReplaceDotWithUnderscore(fOutputBranchNames);
      using ind_t = std::index_sequence_for<ColTypes...>;
      CreateOutputBranches<>(fOutputBranchNames, ind_t{});

      other.fBranchDataToClear.push_back(fBranchData.get());
   }

public:
   using ColumnTypes_t = TypeList<ColTypes...>;
   SnapshotHelperWithVariations(std::string_view filename, std::string_view dirname, std::string_view treename,
                  const ColumnNames_t &/*vbnames*/, const ColumnNames_t &bnames, const RSnapshotOptions &options,
                  std::vector<bool> &&/*isDefine*/, ROOT::Detail::RDF::RLoopManager * outputLoopMgr, ROOT::Detail::RDF::RLoopManager * inputLoopMgr)
      : fOptions(options),
      fOutputBranchNames{ReplaceDotWithUnderscore(bnames)},
      fBranchData{new std::vector<BranchData>{}},
      fBranchDataToClear{fBranchData.get()},
      fInputLoopManager{inputLoopMgr},
      fOutputLoopManager{outputLoopMgr}
   {
      EnsureValidSnapshotTTreeOutput(fOptions, std::string(treename), std::string(filename));

      TFile::TContext fileCtxt;
      fOutputHandle = std::make_shared<FileHandle>(TFile::Open(filename.data(), fOptions.fMode.c_str(), /*ftitle=*/"",
      ROOT::CompressionSettings(fOptions.fCompressionAlgorithm, fOptions.fCompressionLevel)), nullptr);
      if(!fOutputHandle->fFile)
         throw std::runtime_error(std::string{"Snapshot: could not create output file "} + std::string{filename});

      TDirectory *outputDir = fOutputHandle->fFile.get();
      if (!dirname.empty()) {
         fOutputHandle->fDirectoryName = dirname;
         TString checkupdate = fOptions.fMode;
         checkupdate.ToLower();
         if (checkupdate == "update")
            outputDir = fOutputHandle->fFile->mkdir(std::string{dirname}.c_str(), "", true);  // do not overwrite existing directory
         else
            outputDir = fOutputHandle->fFile->mkdir(std::string{dirname}.c_str());
      }

      fOutputHandle->fTree =
         std::make_unique<TTree>(std::string{treename}.c_str(), std::string{treename}.c_str(), fOptions.fSplitLevel, /*dir=*/outputDir);
      fOutputHandle->fOutputLoopManager = fOutputLoopManager;
      if (fOptions.fAutoFlush)
         fOutputHandle->fTree->SetAutoFlush(fOptions.fAutoFlush);

      using ind_t = std::index_sequence_for<ColTypes...>;
      CreateOutputBranches<>(fOutputBranchNames, ind_t{});
   }

   SnapshotHelperWithVariations(SnapshotHelperWithVariations &&) = default;
   ~SnapshotHelperWithVariations()
   {
      // if (!fTreeName.empty() /*not moved from*/ && !fOutputFile /* did not run */ && fOptions.fLazy) {
      //    const auto fileOpenMode = [&]() {
      //       TString checkupdate = fOptions.fMode;
      //       checkupdate.ToLower();
      //       return checkupdate == "update" ? "updated" : "created";
      //    }();
      //    Warning("Snapshot",
      //            "A lazy Snapshot action was booked but never triggered. The tree '%s' in output file '%s' was not %s. "
      //            "In case it was desired instead, remember to trigger the Snapshot operation, by storing "
      //            "its result in a variable and for example calling the GetValue() method on it.",
      //            fTreeName.c_str(), fFileName.c_str(), fileOpenMode);
      // }
   }

   void InitTask(TTreeReader * /*r*/, unsigned int /* slot */)
   {

   }

   void Exec(unsigned int /* slot */, ColTypes &... values)
   {
      using ind_t = std::index_sequence_for<ColTypes...>;
      SetBranches(values..., ind_t{});
   }

   /// Call the Fill of the output tree, and reset all braches to empty values.
   /// This function must be called from exactly one snapshot action.
   /// It triggers the fill of the shared tree at the end of each event.
   TTree& PartialUpdate(unsigned int /*slot*/) {
      if (!fOutputHandle->fTree)
         throw std::runtime_error("The TTree associated to the Snapshot action doesn't exist, any more.");

      fOutputHandle->fTree->Fill();

      for (auto vecPtr : fBranchDataToClear) {
         for (auto & branchData : *vecPtr) branchData.ResetBranchAddressToEmtpyInstance();
      }
      return *fOutputHandle->fTree;
   }

   void Initialize()
   {

   }

   void Finalize()
   {
      fOutputHandle.reset();
   }

   std::string GetActionName() { return "Snapshot"; }

   ROOT::RDF::SampleCallback_t GetSampleCallback() final
   {
      // TODO: Needed?
      return [this](unsigned int, const RSampleInfo &) mutable { ; };
   }

   SnapshotHelperWithVariations MakeNew(void * /*newName*/, std::string_view variation = "nominal")
   {
      return SnapshotHelperWithVariations{*this, "__" + std::string{variation}};
   }
};

}

#endif