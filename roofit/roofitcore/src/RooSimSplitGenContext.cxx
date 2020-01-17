/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
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

/**
\file RooSimSplitGenContext.cxx
\class RooSimSplitGenContext
\ingroup Roofitcore

RooSimSplitGenContext is an efficient implementation of the generator context
specific for RooSimultaneous PDFs when generating more than one of the
component pdfs.
**/

#include "RooFit.h"
#include "Riostream.h"

#include "RooSimSplitGenContext.h"
#include "RooSimultaneous.h"
#include "RooRealProxy.h"
#include "RooDataSet.h"
#include "Roo1DTable.h"
#include "RooCategory.h"
#include "RooMsgService.h"
#include "RooRandom.h"
#include "RooGlobalFunc.h"

using namespace RooFit ;

#include <string>

using namespace std;

ClassImp(RooSimSplitGenContext);
;
  

////////////////////////////////////////////////////////////////////////////////
/// Constructor of specialized generator context for RooSimultaneous p.d.f.s. This
/// context creates a dedicated context for each component p.d.f.s and delegates
/// generation of events to the appropriate component generator context

RooSimSplitGenContext::RooSimSplitGenContext(const RooSimultaneous &model, const RooArgSet &vars, Bool_t verbose, Bool_t autoBinned, const char* binnedTag) :
  RooAbsGenContext(model,vars,0,0,verbose), _pdf(&model)
{
  // Determine if we are requested to generate the index category
  RooAbsCategory *idxCat = (RooAbsCategory*) model._indexCat.absArg() ;
  RooArgSet pdfVars(vars) ;

  RooArgSet allPdfVars(pdfVars) ;

  if (!idxCat->isDerived()) {
    pdfVars.remove(*idxCat,kTRUE,kTRUE) ;
    Bool_t doGenIdx = allPdfVars.find(idxCat->GetName())?kTRUE:kFALSE ;

    if (!doGenIdx) {
      oocoutE(_pdf,Generation) << "RooSimSplitGenContext::ctor(" << GetName() << ") ERROR: This context must"
			       << " generate the index category" << endl ;
      _isValid = kFALSE ;
      // coverity[UNINIT_CTOR]
      return ;
    }
  } else {
    TIterator* sIter = idxCat->serverIterator() ;
    RooAbsArg* server ;
    Bool_t anyServer(kFALSE), allServers(kTRUE) ;
    while((server=(RooAbsArg*)sIter->Next())) {
      if (vars.find(server->GetName())) {
	anyServer=kTRUE ;
	pdfVars.remove(*server,kTRUE,kTRUE) ;
      } else {
	allServers=kFALSE ;
      }
    }
    delete sIter ;    

    if (anyServer && !allServers) {
      oocoutE(_pdf,Generation) << "RooSimSplitGenContext::ctor(" << GetName() << ") ERROR: This context must"
			       << " generate all components of a derived index category" << endl ;
      _isValid = kFALSE ;
      // coverity[UNINIT_CTOR]
      return ;
    }
  }

  // We must extended likelihood to determine the relative fractions of the components
  _idxCatName = idxCat->GetName() ;
  if (!model.canBeExtended()) {
    oocoutE(_pdf,Generation) << "RooSimSplitGenContext::RooSimSplitGenContext(" << GetName() << "): All components of the simultaneous PDF "
			     << "must be extended PDFs. Otherwise, it is impossible to calculate the number of events to be generated per component." << endl ;
    _isValid = kFALSE ;
    // coverity[UNINIT_CTOR]
    return ;
  }

  // Initialize fraction threshold array (used only in extended mode)
  _fracThresh = new Double_t[_pdf->_pdfProxyMap.size()+1] ;
  _fracThresh[0] = 0 ;
  
  // Generate index category and all registered PDFS
  _allVarsPdf.add(allPdfVars) ;
  Int_t i(1) ;
  for (const auto& proxyIter : _pdf->_pdfProxyMap) {
    const auto& pdf = proxyIter.second.arg();

    // Create generator context for this PDF
    RooArgSet* compVars = pdf.getObservables(pdfVars) ;
    RooAbsGenContext* cx = pdf.autoGenContext(*compVars,0,0,verbose,autoBinned,binnedTag) ;
    delete compVars ;

    const auto state = idxCat->lookupIndex(proxyIter.second.name());

    cx->SetName(proxyIter.second.name()) ;
    _gcList.push_back(cx) ;
    _gcIndex.push_back(state);
    
    // Fill fraction threshold array
    _fracThresh[i] = _fracThresh[i-1] + pdf.expectedEvents(&allPdfVars) ;
    i++ ;
  }   

  for (unsigned int j=0 ; j < _pdf->_pdfProxyMap.size() ; j++) {
    _fracThresh[j] /= _fracThresh[_pdf->_pdfProxyMap.size()] ;
  }
    
  // Clone the index category
  _idxCatSet = (RooArgSet*) RooArgSet(model._indexCat.arg()).snapshot(kTRUE) ;
  if (!_idxCatSet) {
    oocoutE(_pdf,Generation) << "RooSimSplitGenContext::RooSimSplitGenContext(" << GetName() << ") Couldn't deep-clone index category, abort," << endl ;
    throw std::string("RooSimSplitGenContext::RooSimSplitGenContext() Couldn't deep-clone index category, abort") ;
  }
  
  _idxCat = (RooAbsCategoryLValue*) _idxCatSet->find(model._indexCat.arg().GetName()) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor. Delete all owned subgenerator contexts

RooSimSplitGenContext::~RooSimSplitGenContext()
{
  delete[] _fracThresh ;
  delete _idxCatSet ;
  for (vector<RooAbsGenContext*>::iterator iter = _gcList.begin() ; iter!=_gcList.end() ; ++iter) {
    delete (*iter) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Attach the index category clone to the given event buffer

void RooSimSplitGenContext::attach(const RooArgSet& args) 
{
  if (_idxCat->isDerived()) {
    _idxCat->recursiveRedirectServers(args,kTRUE) ;
  }

  // Forward initGenerator call to all components
  for (vector<RooAbsGenContext*>::iterator iter = _gcList.begin() ; iter!=_gcList.end() ; ++iter) {
    (*iter)->attach(args) ;
  }
  
}


////////////////////////////////////////////////////////////////////////////////
/// Perform one-time initialization of generator context

void RooSimSplitGenContext::initGenerator(const RooArgSet &theEvent)
{
  // Attach the index category clone to the event
  if (_idxCat->isDerived()) {
    _idxCat->recursiveRedirectServers(theEvent,kTRUE) ;
  } else {
    _idxCat = (RooAbsCategoryLValue*) theEvent.find(_idxCat->GetName()) ;
  }
  
  // Forward initGenerator call to all components
  for (vector<RooAbsGenContext*>::iterator iter = _gcList.begin() ; iter!=_gcList.end() ; ++iter) {
    (*iter)->initGenerator(theEvent) ;
  }

}



////////////////////////////////////////////////////////////////////////////////

RooDataSet* RooSimSplitGenContext::generate(Double_t nEvents, Bool_t skipInit, Bool_t extendedMode)
{
  if(!isValid()) {
    coutE(Generation) << ClassName() << "::" << GetName() << ": context is not valid" << endl;
    return 0;
  }


  // Calculate the expected number of events if necessary
  if(nEvents <= 0) {
    nEvents= _expectedEvents;
  }
  coutI(Generation) << ClassName() << "::" << GetName() << ":generate: will generate "
		    << nEvents << " events" << endl;
    
  if (_verbose) Print("v") ;

  // Perform any subclass implementation-specific initialization
  // Can be skipped if this is a rerun with an identical configuration
  if (!skipInit) {
    initGenerator(*_theEvent);
  }

  // Generate lookup table from expected event counts
  vector<Double_t> nGen(_pdf->_pdfProxyMap.size()) ;
  if (extendedMode ) {
    Int_t i(0) ;
    for (const auto& proxyIter : _pdf->_pdfProxyMap) {
      const auto& pdf = proxyIter.second.arg();
      //nGen[i] = Int_t(pdf->expectedEvents(&_allVarsPdf)+0.5) ;
      nGen[i] = pdf.expectedEvents(&_allVarsPdf) ;
      i++ ;
    }     
    
  } else {
    Int_t i(1) ;
    _fracThresh[0] = 0 ;
    for (const auto& proxyIter : _pdf->_pdfProxyMap) {
      const auto& pdf = proxyIter.second.arg();
      _fracThresh[i] = _fracThresh[i-1] + pdf.expectedEvents(&_allVarsPdf) ;
      i++ ;
    }

    for (unsigned int j=0 ; j < _pdf->_pdfProxyMap.size() ; j++) {
      _fracThresh[j] /= _fracThresh[_pdf->_pdfProxyMap.size()] ;
    }
    
    // Determine from that total number of events to be generated for each component
    Double_t nGenSoFar(0) ;
    while (nGenSoFar<nEvents) {
      Double_t rand = RooRandom::uniform() ;
      for (unsigned int j=0; j < _pdf->_pdfProxyMap.size(); j++) {
        if (rand>_fracThresh[j] && rand<_fracThresh[j+1]) {
          nGen[j]++ ;
          nGenSoFar++ ;
          break ;
        }
      }
    }
  }

    

  // Now loop over states
  map<string,RooAbsData*> dataMap ;
  Int_t icomp(0) ;
  for (const auto& proxyIter : _pdf->_pdfProxyMap) {
    // Calculate number of events to generate for this state    
    if (_gcList[icomp]) {
      dataMap[proxyIter.second.GetName()] = _gcList[icomp]->generate(nGen[icomp],skipInit,extendedMode) ;
    }
    
    icomp++ ;
  }
  
  // Put all datasets together in a composite-store RooDataSet that links and owns the component datasets
  RooDataSet* hmaster = new RooDataSet("hmaster","hmaster",_allVarsPdf,RooFit::Index((RooCategory&)*_idxCat),RooFit::Link(dataMap),RooFit::OwnLinked()) ;    
  return hmaster ;
}



////////////////////////////////////////////////////////////////////////////////
/// Forward to components

void RooSimSplitGenContext::setExpectedData(Bool_t flag) 
{
  for (vector<RooAbsGenContext*>::iterator iter=_gcList.begin() ; iter!=_gcList.end() ; ++iter) {
    (*iter)->setExpectedData(flag) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// this method is empty because it is not used by this context

RooDataSet* RooSimSplitGenContext::createDataSet(const char* , const char* , const RooArgSet& )
{
  return 0 ;
}



////////////////////////////////////////////////////////////////////////////////
/// this method is empty because it is not used in this type of context

void RooSimSplitGenContext::generateEvent(RooArgSet &, Int_t )
{
  assert(0) ;
}




////////////////////////////////////////////////////////////////////////////////
/// this method is empty because proto datasets are not supported by this context

void RooSimSplitGenContext::setProtoDataOrder(Int_t* )
{
  assert(0) ;
}


////////////////////////////////////////////////////////////////////////////////
/// Detailed printing interface

void RooSimSplitGenContext::printMultiline(ostream &os, Int_t content, Bool_t verbose, TString indent) const 
{
  RooAbsGenContext::printMultiline(os,content,verbose,indent) ;
  os << indent << "--- RooSimSplitGenContext ---" << endl ;
  os << indent << "Using PDF ";
  _pdf->printStream(os,kName|kArgs|kClassName,kSingleLine,indent);
}
