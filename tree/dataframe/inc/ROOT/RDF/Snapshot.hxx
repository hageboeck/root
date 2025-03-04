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

#include "TFile.h"
#include "TError.h"

#include <algorithm>
#include <memory>

namespace ROOT::Internal::RDF {
   using namespace ROOT::TypeTraits;
   using namespace ROOT::RDF;

/// Store file and tree in one common place to share them between instances.
struct FileHandle {
   std::unique_ptr<TFile> fFile;
   std::unique_ptr<TTree> fTree;

   FileHandle(TFile * file, TTree * tree) : fFile{file}, fTree{tree} { }
   ~FileHandle() {
      // use AutoSave to flush TTree contents because TTree::Write writes in gDirectory, not in fDirectory
      if (fTree) fTree->AutoSave("flushbaskets");
   }
};

/// Helper object for a single-thread Snapshot action
template <typename... ColTypes>
class R__CLING_PTRCHECK(off) SnapshotHelperWithVariations : public ROOT::Detail::RDF::RActionImpl<SnapshotHelperWithVariations<ColTypes...>> {
   RSnapshotOptions fOptions;
   std::shared_ptr<FileHandle> fOutputHandle;
   // ColumnNames_t fInputBranchNames; // This contains the resolved aliases
   ColumnNames_t fOutputBranchNames;
   // TODO we might be able to unify fBranches, fBranchAddresses and fOutputBranches
   // std::vector<TBranch *> fBranches; // Addresses of branches in output, non-null only for the ones holding C arrays
   // std::vector<void *> fBranchAddresses; // Addresses of objects associated to output branches
   // RBranchSet fOutputBranches;
   // std::vector<bool> fIsDefine;

   // Local buffers for writing nominal or systematics.
   // These need to be stable in memory because TBranch only has the address, so storing them by pointer.
   using BranchBufferTuple_t = std::tuple<ColTypes...>;
   using BranchBufferTupleWithVariations_t = std::tuple<std::vector<ColTypes>...>;
   std::unique_ptr<BranchBufferTuple_t> fBranchBuffers;
   std::unique_ptr<BranchBufferTupleWithVariations_t> fBranchBuffersWithVariations;
   std::vector<BranchBufferTupleWithVariations_t*> fBranchBuffersToClear; /// Branch buffers that must be cleared at the end of the event.

   template<size_t N>
   void CreateOutputBranch(std::string const & name, int bufSize, int splitLevel)
   {
      // TODO: Deal with TBOClRef in original snapshot implementation
      if (fOutputHandle->fTree->GetBranch(name.c_str()))
         throw std::logic_error(std::string{"Snapshot: Output branch "} + name + " already present in tree.");
      if (fBranchBuffers)
         fOutputHandle->fTree->Branch(name.c_str(), &std::get<N>(*fBranchBuffers), bufSize, splitLevel);
      else
         fOutputHandle->fTree->Branch(name.c_str(), &std::get<N>(*fBranchBuffersWithVariations), bufSize, splitLevel);
   }

   template <std::size_t... S>
   void CreateOutputBranches(std::index_sequence<S...>)
   {
      (CreateOutputBranch<S>(fOutputBranchNames[S], fOptions.fBasketSize.value_or(32000), fOptions.fSplitLevel), ...);
   }

   template <std::size_t... S>
   void WriteValuesToBuffers(ColTypes &... values, std::index_sequence<S...>)
   {
      assert((fBranchBuffers || fBranchBuffersWithVariations) && !(fBranchBuffers && fBranchBuffersWithVariations));

      if (fBranchBuffers)
         ((std::get<S>(*fBranchBuffers) = values),...);
      else
         (std::get<S>(*fBranchBuffersWithVariations).push_back(values),...);
   }
   template <std::size_t... S>
   void ClearBranchBuffers(std::index_sequence<S...>)
   {
      for (auto branchBuffer : fBranchBuffersToClear) {
         (std::get<S>(*branchBuffer).clear(),...);
      }
   }

   SnapshotHelperWithVariations(SnapshotHelperWithVariations & other, std::string const & variationPrefix)
   : fOptions{other.fOptions}, fOutputHandle{other.fOutputHandle},
   // fInputBranchNames{other.fInputBranchNames},
   fOutputBranchNames{},
   // fBranches{}, fBranchAddresses{}, fOutputBranches{}, fIsDefine{other.fIsDefine},
   fBranchBuffersWithVariations{new BranchBufferTupleWithVariations_t}
   {
      fOutputBranchNames.resize(other.fOutputBranchNames.size());
      std::transform(other.fOutputBranchNames.begin(), other.fOutputBranchNames.end(), fOutputBranchNames.begin(),
      [&variationPrefix](std::string const & originalName){ return variationPrefix + originalName;});
      other.fBranchBuffersToClear.push_back(fBranchBuffersWithVariations.get());
   }

public:
   using ColumnTypes_t = TypeList<ColTypes...>;
   SnapshotHelperWithVariations(std::string_view filename, std::string_view dirname, std::string_view treename,
                  const ColumnNames_t &/*vbnames*/, const ColumnNames_t &bnames, const RSnapshotOptions &options,
                  std::vector<bool> &&/*isDefine*/)
      : fOptions(options),
      // fInputBranchNames(vbnames),
        fOutputBranchNames(ReplaceDotWithUnderscore(bnames)),
      //   fBranches(vbnames.size(), nullptr),
      //   fBranchAddresses(vbnames.size(), nullptr), fIsDefine(std::move(isDefine)),
        fBranchBuffers{new BranchBufferTuple_t}
   {
      ValidateSnapshotOutput(fOptions, std::string(treename), std::string(filename));

      fOutputHandle = std::make_shared<FileHandle>(TFile::Open(filename.data(), fOptions.fMode.c_str(), /*ftitle=*/"",
      ROOT::CompressionSettings(fOptions.fCompressionAlgorithm, fOptions.fCompressionLevel)), nullptr);
      if(!fOutputHandle->fFile)
         throw std::runtime_error(std::string{"Snapshot: could not create output file "} + std::string{filename});

      TDirectory *outputDir = fOutputHandle->fFile.get();
      if (!dirname.empty()) {
         TString checkupdate = fOptions.fMode;
         checkupdate.ToLower();
         if (checkupdate == "update")
            outputDir = fOutputHandle->fFile->mkdir(std::string{dirname}.c_str(), "", true);  // do not overwrite existing directory
         else
            outputDir = fOutputHandle->fFile->mkdir(std::string{dirname}.c_str());
      }

      fOutputHandle->fTree =
         std::make_unique<TTree>(std::string{treename}.c_str(), std::string{treename}.c_str(), fOptions.fSplitLevel, /*dir=*/outputDir);

      if (fOptions.fAutoFlush)
         fOutputHandle->fTree->SetAutoFlush(fOptions.fAutoFlush);
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
      WriteValuesToBuffers(values..., ind_t{});
   }

   /// Call the Fill of the output tree.
   /// This function must be called from one callback defined for a single snapshot action.
   /// It triggers the fill of the shared tree at the end of each event.
   TTree& PartialUpdate(unsigned int /*slot*/) {
      if (!fOutputHandle->fTree)
         throw std::runtime_error("The TTree associated to the Snapshot action doesn't exist, any more.");

      fOutputHandle->fTree->Fill();
      ClearBranchBuffers(std::index_sequence_for<ColTypes...>{});
      return *fOutputHandle->fTree;
   }

   void Initialize()
   {
      using ind_t = std::index_sequence_for<ColTypes...>;
      CreateOutputBranches(ind_t{});
   }

   void Finalize()
   {
      // TODO: Reset here or in destructor?
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
      // We enter systematics mode, so now write vectors instead of bare values
      if (fBranchBuffers) fBranchBuffers.reset();
      if (!fBranchBuffersWithVariations) {
         fBranchBuffersWithVariations.reset(new BranchBufferTupleWithVariations_t{});
         fBranchBuffersToClear.push_back(fBranchBuffersWithVariations.get());
      }

      return SnapshotHelperWithVariations{*this, std::string{variation} + ":"};
   }
};

}

#endif