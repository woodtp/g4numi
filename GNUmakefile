# $Id: GNUmakefile,v 1.2 2005/04/05 20:29:13 zarko Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := g4numi
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

  CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)

  ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --glibs) -lMinuit -lHtml
  ROOTLIBS      := $(filter-out -lNew,$(ROOTLIBS))
  ROOTLIBS      := $(filter-out -lThread,$(ROOTLIBS))
  ROOTLIBS      := $(filter-out -lpthread,$(ROOTLIBS))
  INTYLIBS      += $(ROOTLIBS)

# Root (exlude libNew and libpthread from library list)

#ROOTINC       = -I$(ROOTSYS)/include

#ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --glibs) -lMinuit -lHtml
#ROOTLIBS      := $(filter-out -lNew,$(ROOTLIBS))
#ROOTLIBS      := $(filter-out -lpthread,$(ROOTLIBS))

# Extra flags for G4
#CPPFLAGS += $(ROOTINC)
#LDLIBS   += $(ROOTLIBS)

include $(G4INSTALL)/config/architecture.gmk

include hadronic_lists.gmk

#G4ANALYSIS_USE:= true

ifdef G4ANALYSIS_USE
  CPPFLAGS += `aida-config --include`
  LDFLAGS += `aida-config --lib`
endif


.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

