// Author: Enrico Guiraud, Danilo Piparo CERN  9/2017

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RROOTTDS
#define ROOT_RROOTTDS

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDataSource.hxx"
#include <TChain.h>

#include <memory>

namespace ROOT {

namespace Internal {

namespace RDF {

/// This class is unused and it has only been implemented as a proof of concept.
/// It shows how to implement the RDataSource API for a complex kind of source such as TTrees.
class RRootDS final : public ROOT::RDF::RDataSource {
private:
   std::string fTreeName;
   std::string fFileNameGlob;
   mutable TChain fModelChain; // Mutable needed for getting the column type name
   std::vector<double *> fAddressesToFree;
   std::vector<std::string> fListOfBranches;
   std::vector<std::pair<ULong64_t, ULong64_t>> fEntryRanges;
   std::vector<std::vector<void *>> fBranchAddresses; // first container-> slot, second -> column;
   std::vector<std::unique_ptr<TChain>> fChains;

   std::vector<void *> GetColumnReadersImpl(std::string_view, const std::type_info &) final;

protected:
   std::string AsString() final { return "ROOT data source"; };

public:
   RRootDS(std::string_view treeName, std::string_view fileNameGlob);
   // Rule of five
   RRootDS(const RRootDS &) = delete;
   RRootDS &operator=(const RRootDS &) = delete;
   RRootDS(RRootDS &&) = delete;
   RRootDS &operator=(RRootDS &&) = delete;
   ~RRootDS() final;

   std::size_t GetNFiles() const final;
   std::string GetTypeName(std::string_view colName) const final;
   const std::vector<std::string> &GetColumnNames() const final;
   bool HasColumn(std::string_view colName) const final;
   void InitSlot(unsigned int slot, ULong64_t firstEntry) final;
   void FinalizeSlot(unsigned int slot) final;
   std::vector<std::pair<ULong64_t, ULong64_t>> GetEntryRanges() final;
   bool SetEntry(unsigned int slot, ULong64_t entry) final;
   void SetNSlots(unsigned int nSlots) final;
   void Initialize() final;
   std::string GetLabel() final;
};

} // ns RDF

} // ns Internal

} // ns ROOT

#endif
