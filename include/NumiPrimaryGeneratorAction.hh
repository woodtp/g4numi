
#ifndef NumiPrimaryGeneratorAction_h
#define NumiPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class NumiDataInput;
class TFile;
class TTree;
class NumiRunManager;

class NumiPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
  NumiPrimaryGeneratorAction();
  ~NumiPrimaryGeneratorAction();

 public:
  void GeneratePrimaries(G4Event* anEvent);
  void SetProtonBeam();
  void SetMuonBeam();
  G4bool OpenNtuple(G4String ntupleName);
  void CloseNtuple();
  
  // Primary proton information
  G4ThreeVector GetProtonOrigin(){
    return fProtonOrigin;
  }
  G4ThreeVector GetProtonMomentum(){
    return fProtonMomentum;
  }
  G4ThreeVector GetProtonIntVertex(){
    return fProtonIntVertex;
  }
  // Info about a particle leaving the target (when using external ntuple)
  G4ThreeVector GetParticlePosition(){
    return fParticlePosition;
  }
  G4ThreeVector GetParticleMomentum(){
    return fParticleMomentum;
  }
  G4double GetWeight(){
    return fWeight;
  }
  G4int GetTgen(){
    return fTgen;
  }
  G4int GetParticleType(){
    return fType;
  }
  // *********************************************************************
  G4int GetNoOfPrimaries(){
    return fNoOfPrimaries;
  }    
 

 private:
  G4ThreeVector fParticleMomentum,fParticlePosition;    
  G4int fNoOfPrimaries,fTgen,fType; 
  G4double fWeight;
  G4bool fIsFirst;
  G4ThreeVector tunnelPos;
    
  G4ThreeVector fProtonOrigin;
  G4ThreeVector fProtonMomentum;
  G4ThreeVector fProtonIntVertex;

  TFile *fRootFile;
  TTree *fPrimaryNtuple;
  G4int fCurrentPrimaryNo;
  NumiDataInput* fND;
  NumiRunManager* fRunManager;
  G4ParticleGun* fParticleGun;

};

#endif


