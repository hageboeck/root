// Author: Stephan Hageboeck, CERN  01/2019

/*****************************************************************************
 * RooFit
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2019, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/


#ifndef ROOFIT_ROOFITCORE_INC_ROOWORKSPACEHANDLE_H_
#define ROOFIT_ROOFITCORE_INC_ROOWORKSPACEHANDLE_H_

#include "RooWorkspace.h"

///////////////////////////////////////////////////////////////////////////////////////////////
/// An interface to set and retrieve a workspace.
/// This is needed for all generic objects that can be saved in a workspace, which itself depend
/// on the workspace (e.g. the RooStats::ModelConfig).
/// Because of a circular dependency, a workspace with a ModelConfig cannot be (deep) cloned.
/// The handle hides this dependency.
class RooWorkspaceHandle {
public:
   virtual ~RooWorkspaceHandle() = default;

   ///Set the workspace. If it exists, it is up to the implementing class to decide how to proceed.
   virtual void SetWS(RooWorkspace &ws) = 0;

   ///Set the workspace irrespective of what the previous workspace is.
   virtual void ReplaceWS(RooWorkspace *ws) = 0;

   ///Retrieve the workspace.
   virtual RooWorkspace *GetWS() const = 0;

   ClassDef(RooWorkspaceHandle, 0)
};

#endif /* ROOFIT_ROOFITCORE_INC_ROOWORKSPACEHANDLE_H_ */
