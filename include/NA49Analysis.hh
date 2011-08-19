
#ifndef NA49ANALYSIS_HH
#define NA49ANALYSIS_HH

#include "globals.hh"

//root
#include "TSystem.h"

//G4
#include "G4ios.hh"
#include "G4TrajectoryContainer.hh"
#include "G4ParticleDefinition.hh"
#include "G4Element.hh"

#include <map.h>
#include <map>
#include <vector>

class G4ParticleDefinition;
class G4Step;
class TFile;
class TTree;
class G4Track;
class G4VTrajectory;

class ProdTuple_t;

class NA49Analysis
{
public:

   NA49Analysis();
   ~NA49Analysis();

   void book();
   void FillNtuple(const G4Track& track);
   void WriteNtuple();
   static NA49Analysis* getInstance();
   void finish();
   
private:
   
private:
   static NA49Analysis* instance;
   
  char NtupleFileName[50];
   TFile* FileNtuple;  
   TTree* ProdTree;
   ProdTuple_t* g4Proddata; 

};

#endif 
