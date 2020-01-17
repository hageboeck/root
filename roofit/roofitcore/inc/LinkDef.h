#include "roofit/roofitcore/inc/LinkDef1.h"
#include "roofit/roofitcore/inc/LinkDef2.h"
#include "roofit/roofitcore/inc/LinkDef3.h"
#include "roofit/roofitcore/inc/LinkDef4.h"

#pragma link C++ class RooSTLRefCountList<RooAbsArg>+;
#pragma read sourceClass="RooDataHist" targetClass="RooDataHist" version="[3-4]" \
source="int _arrSize; double* _wgt;" \
target="_wgtVec" \
include="TVirtualStreamerInfo.h" \
code="{ TClass::GetClass(\"RooDataHist@@4\")->GetStreamerInfo()->ls(); \
        _wgtVec.assign(onfile._wgt, onfile._wgt + onfile._arrSize); }"
