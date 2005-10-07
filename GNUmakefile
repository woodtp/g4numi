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
else
  EXTRALIBS += -L$(G4LISTS_BASE)/hadronic/plists/lib/$(G4SYSTEM)
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
#
# Select your physics lists to link against.
#
# EXTRALIBS += -lFTFC
# EXTRALIBS += -lFTFP
# EXTRALIBS += -lLBE
# EXTRALIBS += -lLHEP
# EXTRALIBS += -lLHEP_GN
# EXTRALIBS += -lLHEP_HP
# EXTRALIBS += -lLHEP_LEAD
# EXTRALIBS += -lLHEP_BERT_HP
# EXTRALIBS += -lLHEP_BIC_HP
# EXTRALIBS += -lLHEP_LEAD_HP
# EXTRALIBS += -lLHEP_PRECO
# EXTRALIBS += -lQGSP_BERT
# EXTRALIBS += -lLHEP_PRECO_HP
# EXTRALIBS += -lQGSC
# EXTRALIBS += -lQGSC_LEAD
# EXTRALIBS += -lQGSC_LEAD_HP
EXTRALIBS += -lQGSP
# EXTRALIBS += -lQGSP_GN
# EXTRALIBS += -lQGSP_HP
# EXTRALIBS += -lLHEP_BERT
# EXTRALIBS += -lLHEP_BIC
# EXTRALIBS += -lQGSP_BIC
EXTRALIBS += -lPackaging

EXTRALIBS += -lG4hadronic_proc
EXTRALIBS += -lG4hadronic_HE
EXTRALIBS += -lG4hadronic_LE
EXTRALIBS += -lG4hadronic_iso
EXTRALIBS += -lG4had_neu_hp
EXTRALIBS += -lG4hadronic_coherent_elastic
EXTRALIBS += -lG4hadronic_hetcpp_evaporation
EXTRALIBS += -lG4hadronic_hetcpp_utils
EXTRALIBS += -lG4hadronic_bert_cascade
EXTRALIBS += -lG4hadronic_interface_ci
EXTRALIBS += -lG4hadronic_body_ci
EXTRALIBS += -lG4hadronic_leading_particle
EXTRALIBS += -lG4hadronic_stop
EXTRALIBS += -lG4hadronic_radioactivedecay

EXTRALIBS += -lG4had_theo_max 
EXTRALIBS += -lG4hadronic_qgstring
EXTRALIBS += -lG4had_string_diff
EXTRALIBS += -lG4had_string_frag
EXTRALIBS += -lG4had_string_man

EXTRALIBS += -lG4hadronic_binary 
EXTRALIBS += -lG4had_im_r_matrix 
EXTRALIBS += -lG4had_preequ_exciton 
EXTRALIBS += -lG4hadronic_deex_evaporation 
EXTRALIBS += -lG4hadronic_deex_fermi_breakup
EXTRALIBS += -lG4hadronic_deex_fission
EXTRALIBS += -lG4hadronic_deex_gem_evaporation
EXTRALIBS += -lG4hadronic_deex_handler
EXTRALIBS += -lG4hadronic_deex_management
EXTRALIBS += -lG4hadronic_deex_multifragmentation
EXTRALIBS += -lG4hadronic_deex_photon_evaporation
EXTRALIBS += -lG4hadronic_deex_util

EXTRALIBS += -lG4had_mod_man 
EXTRALIBS += -lG4had_mod_util

EXTRALIBS += -lG4hadronic_mgt
EXTRALIBS += -lG4hadronic_xsect
EXTRALIBS += -lG4hadronic_util

EXTRALIBS += -lG4shortlived

CPPFLAGS += -I$(G4LISTS_BASE)/lists/FTFC/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/FTFP/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/LBE/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/LHEP/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/LHEP_BERT/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/LHEP_BIC/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/LHEP_GN/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/LHEP_HP/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/LHEP_BERT_HP/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/LHEP_BIC_HP/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/LHEP_LEAD/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/LHEP_LEAD_HP/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/LHEP_PRECO/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/LHEP_PRECO_HP/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/Packaging/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/QGSC/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/QGSC_LEAD/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/QGSC_LEAD_HP/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/QGSP/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/QGSP_BERT/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/QGSP_BIC/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/QGSP_GN/include
CPPFLAGS += -I$(G4LISTS_BASE)/lists/QGSP_HP/include


LISTS_INCLUDE = $(G4LISTS_BASE)
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/FTFC/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/FTFP/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LBE/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_BERT/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_BIC/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_GN/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_HP/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_BERT_HP/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_BIC_HP/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_LEAD/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_LEAD_HP/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_PRECO/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LHEP_PRECO_HP/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/Packaging/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/QGSC/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/QGSC_LEAD/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/QGSC_LEAD_HP/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/QGSP/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/QGSP_BERT/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/QGSP_BIC/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/QGSP_GN/include
CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/QGSP_HP/include

G4BASE = $(G4INSTALL)/source
CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/geometry/navigation/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/interface/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/event/include \
            -I$(G4BASE)/particles/shortlived/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/run/include \
            -I$(G4BASE)/tracking/include \
            -I$(G4BASE)/digits_hits/hits/include \
            -I$(G4BASE)/digits_hits/digits/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/decay/include \
            -I$(G4BASE)/processes/cuts/include \
            -I$(G4BASE)/processes/optical/include \
            -I$(G4BASE)/processes/transportation/include \
            -I$(G4BASE)/processes/electromagnetic/lowenergy/include \
            -I$(G4BASE)/processes/electromagnetic/standard/include \
            -I$(G4BASE)/processes/electromagnetic/muons/include \
            -I$(G4BASE)/processes/electromagnetic/utils/include \
            -I$(G4BASE)/processes/electromagnetic/xrays/include \
            -I$(G4BASE)/processes/hadronic/cross_sections/include \
            -I$(G4BASE)/processes/hadronic/stopping/include \
            -I$(G4BASE)/processes/hadronic/management/include \
            -I$(G4BASE)/processes/hadronic/processes/include \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(LISTS_INCLUDE)/hadronic/Packaging/include \
            -I$(G4BASE)/processes/hadronic/models/management/include \
            -I$(G4BASE)/processes/hadronic/models/util/include \
            -I$(G4BASE)/processes/hadronic/models/binary_cascade/include \
            -I$(G4BASE)/processes/hadronic/models/cascade/cascade/include \
            -I$(G4BASE)/processes/hadronic/models/cascade/evaporation/include \
            -I$(G4BASE)/processes/hadronic/models/cascade/utils/include \
            -I$(G4BASE)/processes/hadronic/models/chiral_inv_phase_space/body/include \
            -I$(G4BASE)/processes/hadronic/models/chiral_inv_phase_space/interface/include \
            -I$(G4BASE)/processes/hadronic/models/coherent_elastic/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/evaporation/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/fermi_breakup/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/fission/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/gem_evaporation/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/handler/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/management/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/multifragmentation/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/photon_evaporation/include \
            -I$(G4BASE)/processes/hadronic/models/de_excitation/util/include \
            -I$(G4BASE)/processes/hadronic/models/high_energy/include \
            -I$(G4BASE)/processes/hadronic/models/im_r_matrix/include \
            -I$(G4BASE)/processes/hadronic/models/isotope_production/include \
            -I$(G4BASE)/processes/hadronic/models/leading_particle/include \
            -I$(G4BASE)/processes/hadronic/models/low_energy/include \
            -I$(G4BASE)/processes/hadronic/models/neutron_hp/include \
            -I$(G4BASE)/processes/hadronic/models/parton_string/diffraction/include \
            -I$(G4BASE)/processes/hadronic/models/parton_string/hadronization/include \
            -I$(G4BASE)/processes/hadronic/models/parton_string/management/include \
            -I$(G4BASE)/processes/hadronic/models/parton_string/qgsm/include \
            -I$(G4BASE)/processes/hadronic/models/pre_equilibrium/exciton_model/include \
            -I$(G4BASE)/processes/hadronic/models/radioactive_decay/include \
            -I$(G4BASE)/processes/hadronic/models/theo_high_energy/include \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/materials/include


.PHONY: all
all: g4numiCint lib  libg4numiDict bin 

include $(G4INSTALL)/config/binmake.gmk

CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))

g4numiCint: include/data_t.hh include/hadmmtuple_t.hh  Linkdef.h
	rootcint -f ./src/g4numiCint.cc -c -I./include ../include/data_t.hh ../include/hadmmtuple_t.hh ../Linkdef.h

libg4numiDict:  $(G4WORKDIR)/tmp/$(G4SYSTEM)/g4numi/data_t.o   $(G4WORKDIR)/tmp/$(G4SYSTEM)/g4numi/hadmmtuple_t.o $(G4WORKDIR)/tmp/$(G4SYSTEM)/g4numi/g4numiCint.o
	gcc -g -shared -o libg4numiDict.so    $(G4WORKDIR)/tmp/$(G4SYSTEM)/g4numi/data_t.o   $(G4WORKDIR)/tmp/$(G4SYSTEM)/g4numi/hadmmtuple_t.o $(G4WORKDIR)/tmp/$(G4SYSTEM)/g4numi/g4numiCint.o

