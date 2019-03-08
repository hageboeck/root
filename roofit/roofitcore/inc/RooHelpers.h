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

#ifndef ROOFIT_ROOFITCORE_INC_ROOHELPERS_H_
#define ROOFIT_ROOFITCORE_INC_ROOHELPERS_H_

#include "RooMsgService.h"
#include "RooAbsArg.h"

#include <sstream>

namespace RooHelpers {

/// Switches the message service to verbose while the instance alive.
class MakeVerbose {
  public:
    MakeVerbose() {
      auto& msg = RooMsgService::instance();
      fOldKillBelow = msg.globalKillBelow();
      msg.setGlobalKillBelow(RooFit::DEBUG);
      fOldConf = msg.getStream(0);
      msg.getStream(0).minLevel= RooFit::DEBUG;
      msg.setStreamStatus(0, true);
    }

    ~MakeVerbose() {
      auto& msg = RooMsgService::instance();
      msg.setGlobalKillBelow(fOldKillBelow);
      msg.getStream(0) = fOldConf;
      msg.setStreamStatus(0, true);
    }

  private:
    RooFit::MsgLevel fOldKillBelow;
    RooMsgService::StreamConfig fOldConf;
};


/// Hijacks all messages with given level and topic (and optionally object name) while alive.
/// Use like ostringstream afterwards. Useful for unit tests and debugging.
class HijackMessageStream : public std::ostringstream {
  public:
    HijackMessageStream(RooFit::MsgLevel level, RooFit::MsgTopic topics, const char* objectName = nullptr);

    virtual ~HijackMessageStream();

  private:
    RooFit::MsgLevel _oldKillBelow;
    std::vector<RooMsgService::StreamConfig> _oldConf;
    Int_t _thisStream;
};


std::vector<std::string> tokenise(const std::string &str, const std::string &delims);



class CachingError : public std::logic_error {
  public:
    CachingError(const std::string& indent, const std::string& str) :
    std::logic_error(str),
    _indent{indent} { }


    std::string indent() const {
      return _indent + " ";
    }

  private:
    std::string _indent;
};

class FormatPdfTree {
  public:
    FormatPdfTree& operator<<(const CachingError& arg) {
      _stream << arg.what() << "\n" << arg.indent();
      return *this;
    }

    template <class T,
    typename std::enable_if<std::is_base_of<RooAbsArg, T>::value>::type* = nullptr >
    FormatPdfTree& operator<<(const T& arg) {
      _stream << arg.ClassName() << "::" << arg.GetName() << " " << &arg << " ";
      arg.printArgs(_stream);
      return *this;
    }

    template <class T,
    typename std::enable_if< ! std::is_base_of<RooAbsArg, T>::value>::type* = nullptr >
    FormatPdfTree& operator<<(const T& arg) {
      _stream << arg;
      return *this;
    }

    operator std::string() const {
      return _stream.str();
    }

  private:
    std::ostringstream _stream;
};

}

#endif /* ROOFIT_ROOFITCORE_INC_ROOHELPERS_H_ */
