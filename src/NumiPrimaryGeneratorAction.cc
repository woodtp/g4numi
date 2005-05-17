
#include "NumiPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "NumiDataInput.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh"

NumiPrimaryGeneratorAction::NumiPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  ND=NumiDataInput::GetNumiDataInput();
  particleGun = new G4ParticleGun(n_particle);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  particleGun->SetParticleDefinition(particleTable->FindParticle("proton"));
  particleGun->SetParticleEnergy(ND->protonKineticEnergy);
  particleGun->SetParticlePosition(ND->beamPosition);
  particleGun->SetParticleMomentumDirection(ND->beamDirection);

  G4cout << "Proton Momentum: "<<ND->protonMomentum/GeV << "GeV" <<G4endl;
  G4cout << "Proton Kinetic Energy: "<<ND->protonKineticEnergy/GeV<<"GeV" <<G4endl;
  G4cout << "Proton BeamDirection: "<<ND->beamDirection <<G4endl;  
}

NumiPrimaryGeneratorAction::~NumiPrimaryGeneratorAction()
{
  delete particleGun;
}

void NumiPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4UImanager* UI=G4UImanager::GetUIpointer();
  G4bool reGet=true;

  if (anEvent->GetEventID()==0){
    meanx = UI->GetCurrentDoubleValue("/gun/position",1,reGet)*cm;
    meany = UI->GetCurrentDoubleValue("/gun/position",2,reGet)*cm;
    meanz = UI->GetCurrentDoubleValue("/gun/position",3,reGet)*cm;
  }
  G4double x0;
  G4double y0; 
  G4double z0;
  G4double sigmax=ND->beamSigmaX;
  G4double sigmay=ND->beamSigmaY;
  //  G4cout<<"xpos = " <<meanx/m<< " ypos = "<<meany/m<<" zpos = "<< meanz/m<<G4endl;
  x0 = G4RandGauss::shoot(ND->beamPosition[0],sigmax);
  y0 = G4RandGauss::shoot(ND->beamPosition[1],sigmay);
  z0 = ND->beamPosition[2];

  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  particleGun->GeneratePrimaryVertex(anEvent);
}


