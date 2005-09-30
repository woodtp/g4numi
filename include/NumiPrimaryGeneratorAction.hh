
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
  
  G4double GetProtonX(){
    return protonX;
  }
  G4double GetProtonY(){
    return protonY;
  }
  G4ThreeVector GetProtonDirection(){
    return protonDirection;
  }
  
    G4int noOfPrimaries,tgen,type;
    G4double weight;
    G4ThreeVector ParticleMomentum,ParticlePosition;

  private:
    G4bool isFirst;
    G4double protonX,protonY;
    G4ThreeVector protonDirection;   
    TFile *rootFile;
    TTree *primaryNtuple;
    G4int currentPrimaryNo;
    NumiDataInput* ND;
    NumiAnalysis* numiAnalysis;
    NumiRunManager* pRunManager;
    G4ParticleGun* particleGun;
    G4double meanx;
    G4double meany;
    G4double meanz;
};

#endif


