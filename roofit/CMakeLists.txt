# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

add_subdirectory(roofitcore)
add_subdirectory(roofit)
if(mathmore)
  add_subdirectory(roofitmore)
endif()
add_subdirectory(roostats)
if(xml)
  add_subdirectory(histfactory)
endif()

# For recent clang, this can facilitate auto-vectorisation.
# In RooFit, the errno side effect is not needed, anyway:
target_compile_options(RooFitCore PUBLIC -fno-math-errno)
target_compile_options(RooFit     PUBLIC -fno-math-errno)
if(TARGET RooFitMore)
  target_compile_options(RooFitMore PUBLIC -fno-math-errno)
endif()
if(TARGET HistFactory)
  target_compile_options(HistFactory PUBLIC -fno-math-errno)
endif()
