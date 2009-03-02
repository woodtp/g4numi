
#ifndef NumiPrimaryGeneratorAction_h
#define NumiPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class NumiPrimaryMessenger;
class NumiDataInput;
class TFile;
class TTree;
class NumiRunManager;

#include "NtpMuon.hh"

class NumiPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
  NumiPrimaryGeneratorAction();
  ~NumiPrimaryGeneratorAction();

 public:
  void GeneratePrimaries(G4Event* anEvent);
  G4ParticleGun* GetParticleGun() {return fParticleGun;};
  void SetProtonBeam();
  void SetMuonBeam();
  G4bool OpenNtuple(G4String ntupleName);
  void CloseNtuple();

  // Parameters from Messenger that allow input to mirror
  // FLUKA input in macro file
  void PG_CosX(G4double _cx) { cosx = _cx; };
  void PG_CosY(G4double _cy) { cosy = _cy; };
  void PG_Z(G4double _z0) { z0 = _z0; };
  void PG_Momentum(G4double _p) { p = _p; };
  void PG_Spread(G4double _spread) { spread = _spread; };
  void PG_Divergence(G4double _div) { div = _div; };
  void PG_InnerR(G4double _innerR) { inR = _innerR; };
  void PG_OuterR(G4double _outerR) { outR = _outerR; };
  void PG_Particle(G4ParticleDefinition* _particle) { particle = _particle; };


  
   // Primary proton information
   G4ThreeVector GetProtonOrigin()         { return fProtonOrigin; }
   G4ThreeVector GetProtonMomentum()       { return fProtonMomentum; }
   G4ThreeVector GetProtonIntVertex()      { return fProtonIntVertex; }
   // Info about a particle leaving the target (when using external ntuple)
   G4ThreeVector GetParticlePosition()     { return fParticlePosition; }
   G4ThreeVector GetParticleMomentum()     { return fParticleMomentum; }
   G4double GetWeight()                    { return fWeight; }
   G4int GetTgen()                         { return fTgen; }
   G4int GetParticleType()                 { return fType; }
   G4double GetMuWeight()                  { return fmuweight; }
   G4ThreeVector GetMuTParentMomentum()    { return fMuTParentMomentum; }
   G4ThreeVector GetMuParentMomentum()     { return fMuParentMomentum; }
   G4ThreeVector GetMuTParentPosition()    { return fMuTParentPosition; }
   G4ThreeVector GetMuParentPosition()     { return fMuParentPosition; }
   G4ThreeVector GetMuParentProdPosition() { return fMuParentProdPosition; }
   G4int GetTParentType()                  { return ftpptype; }
   G4double GetImpWeight()                 { return fnimpwt; }
   G4int GetMuParentType()                 { return fpptype; }
   G4int GetMuParentProdMedium()           { return fppmedium; }
   G4int GetMuParentGen()                  { return fpgen; }
   G4int GetEvtno()                        { return fevtno; }
   G4int GetNoOfPrimaries()                { return fNoOfPrimaries; }    
   

 private:
  double                  DoubleRand() {return 2*G4UniformRand()-1.;}

  G4ThreeVector fParticleMomentum,fParticlePosition, flocalParticleMomentum;    
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

  // Inputs (as in ToyDP.inp)
  NumiPrimaryMessenger*       beamMessenger;
  G4double cosx;
  G4double cosy;
  G4double z0;
  G4double p;
  G4double spread;
  G4double div;
  G4double outR;
  G4double inR;
  G4ParticleDefinition* particle;
  G4String part;

  // Muon Inputs
   G4double fmuweight;
   G4ThreeVector fMuTParentMomentum;
   G4ThreeVector fMuParentMomentum;
   G4ThreeVector fMuTParentPosition;
   G4ThreeVector fMuParentPosition;
   G4ThreeVector fMuParentProdPosition;
   G4int ftpptype;
   G4double fnimpwt;
   G4int fpptype;
   G4int fppmedium;
   G4int fpgen;
   G4int fevtno;
   
   NtpMuon* fMuon;

};

#endif


