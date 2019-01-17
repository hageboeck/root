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
\file RooNumber.cxx
\class RooNumber
\ingroup Roofitcore

Class RooNumber implements numeric constants used by RooFit
**/

#include "RooFit.h"
#include "RooNumber.h"

#include <limits>

using namespace std;

ClassImp(RooNumber);
;

Double_t RooNumber::_Infinity = numeric_limits<Double_t>::infinity();


////////////////////////////////////////////////////////////////////////////////
/// Return internal infinity representation

Double_t RooNumber::infinity() 
{
  return _Infinity ;
}


////////////////////////////////////////////////////////////////////////////////
/// Return true if x is infinite by RooNumBer internal specification

Int_t RooNumber::isInfinite(Double_t x) 
{
  return (x >= +_Infinity) ? +1 : ((x <= -_Infinity) ? -1 : 0);
}

