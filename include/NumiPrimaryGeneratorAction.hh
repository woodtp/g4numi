
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
class NumiAnalysis;
class NumiRunManager;

class NumiPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    NumiPrimaryGeneratorAction();
    ~NumiPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);
    void SetProtonBeam();
    G4bool OpenNtuple(G4String ntupleName);
    void CloseNtuple();
  
  G4ThreeVector GetProtonOrigin(){
    return protonOrigin;
  }
  G4ThreeVector GetProtonMomentum(){
    return protonMomentum;
  }
  G4ThreeVector GetProtonIntVertex(){
    return protonIntVertex;
  }
  
    G4int noOfPrimaries,tgen,type;
    G4double weight;
    G4ThreeVector ParticleMomentum,ParticlePosition;

  private:
    G4bool isFirst;
    G4ThreeVector protonOrigin;
    G4ThreeVector protonMomentum;
    G4ThreeVector protonIntVertex;

    TFile *rootFile;
    TTree *primaryNtuple;
    G4int currentPrimaryNo;
    NumiDataInput* ND;
    NumiAnalysis* numiAnalysis;
    NumiRunManager* pRunManager;
    G4ParticleGun* particleGun;

};

#endif


