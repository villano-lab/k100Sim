# --------------------------------------------------------------
# GNUmakefile for physics list user.  
# JPW. Fri Jul 25 10:39:58 CEST 2003
# --------------------------------------------------------------

name := k100
G4TARGET := $(name)
G4EXLIB := true 

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

#ifdef G4EXPAT_PATH
#  EXTRALIBS += -L$(G4EXPAT_PATH)
#endif

include $(G4INSTALL)/config/architecture.gmk

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
#CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
#CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))

