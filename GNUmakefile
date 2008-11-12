# --------------------------------------------------------------
# GNUmakefile for physics list user.  
# JPW. Fri Jul 25 10:39:58 CEST 2003
# --------------------------------------------------------------

name := g4numi

G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

#
# define G4LISTS_BASE, if you have your own physics lists area installed
# point G4LISTS_BASE to the directory, that contains the subdirectory 'lists'.
#
ifndef G4LISTS_BASE
  EXTRALIBS += -L$(G4LIB)/.lists_build/$(G4SYSTEM)
  G4LISTS_BASE = $(G4INSTALL)/hadronic_lists
  EXTRALIBS += -L$(G4LIB)/plists/$(G4SYSTEM)	
else
  EXTRALIBS += -L$(G4LISTS_BASE)/hadronic/plists/lib/$(G4SYSTEM)
  EXTRALIBS += -L$(G4LIB)/plists/$(G4SYSTEM)
endif

#for root
#CPPFLAGS   +=  -I$(ROOTSYS)/include
  CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)

  ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --glibs) -lMinuit -lHtml
  ROOTLIBS      := $(filter-out -lNew,$(ROOTLIBS))
  ROOTLIBS      := $(filter-out -lThread,$(ROOTLIBS))
  ROOTLIBS      := $(filter-out -lpthread,$(ROOTLIBS))
  INTYLIBS      += $(ROOTLIBS)

#for debuging
   CPPFLAGS += -g	

EXTRALIBS += -lQGSP
EXTRALIBS += -lPackaging

CPPFLAGS += -I$(G4INSTALL)/physics_lists/hadronic/QGSP/include
CPPFLAGS += -I$(G4INSTALL)/physics_lists/hadronic/Packaging/include
CPPFLAGS += -I$(G4INCLUDE)


.PHONY: all
all: g4numiCint lib  libg4numiDict bin 
	cp $(G4WORKDIR)/bin/Linux-g++/$(G4TARGET) .

include $(G4INSTALL)/config/binmake.gmk

CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))

g4numiCint: include/data_t.hh include/hadmmtuple_t.hh  Linkdef.h
	rootcint -f ./src/g4numiCint.cc -c -I./include ../include/data_t.hh ../include/hadmmtuple_t.hh ../Linkdef.h

libg4numiDict:  $(G4WORKDIR)/tmp/$(G4SYSTEM)/g4numi/data_t.o   $(G4WORKDIR)/tmp/$(G4SYSTEM)/g4numi/hadmmtuple_t.o $(G4WORKDIR)/tmp/$(G4SYSTEM)/g4numi/g4numiCint.o
	gcc -g -shared -o libg4numiDict.so    $(G4WORKDIR)/tmp/$(G4SYSTEM)/g4numi/data_t.o   $(G4WORKDIR)/tmp/$(G4SYSTEM)/g4numi/hadmmtuple_t.o $(G4WORKDIR)/tmp/$(G4SYSTEM)/g4numi/g4numiCint.o
