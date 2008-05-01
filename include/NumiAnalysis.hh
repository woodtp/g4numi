//
// NumiAnalysis.hh
//
// Modified Jul 2005 by A. MArino to make data_t and hadmmtuple_t classes

#ifndef NUMIANALYSIS_HH
#define NUMIANALYSIS_HH

#include "globals.hh"
#include "NumiTrajectory.hh"
#include "NumiDataInput.hh"

//root
#include "TSystem.h"

//G4
#include "G4ios.hh"
#include "G4TrajectoryContainer.hh"

class TFile;
class TTree;
class G4Track;
class data_t;
class hadmmtuple_t;

class NumiAnalysis
{
public:

  NumiAnalysis();
  ~NumiAnalysis();

  void book();
  void finish();
  void FillNeutrinoNtuple(const G4Track& track);
  void FillHadmmNtuple(const G4Track& track,Int_t mm_num,Int_t cellNum);
  void FillHadmmNtuple();
  void WriteHadmmNtuple();
  NumiTrajectory* GetParentTrajectory(G4int parentID);
  static NumiAnalysis* getInstance();

  void SetCount(G4int count);
  G4int GetCount();
  void SetEntry(G4int entry);
  G4int GetEntry();

private:
  std::string GetOFileName(std::string ifilename);

private:
  static NumiAnalysis* instance;

  G4double x;
  G4double y;
  G4double z;

  G4double noProtons;
  char asciiFileName[50], nuNtupleFileName[50], hadmmNtupleFileName[50];

  NumiDataInput* NumiData;

  TFile* hadmmNtuple;
  TFile* nuNtuple;
 
  TTree* tree;
  TTree* hadmmtree;

  data_t *g4data;
  hadmmtuple_t *g4hmmdata;

  G4int fcount;
  G4int fentry;

};

#endif 
