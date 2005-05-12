
#include "NumiPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"

NumiPrimaryGeneratorAction::NumiPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  particleGun->SetParticleDefinition(particleTable->FindParticle("proton"));
  particleGun->SetParticleEnergy(120.0*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, -3.0*m));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.0,0.0,1.0));
}

NumiPrimaryGeneratorAction::~NumiPrimaryGeneratorAction()
{
  delete particleGun;
}

void NumiPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4double x0;
  G4double y0;
  G4double z0;
  G4double meanx=0.; G4double sigmax=0.9*mm;
  G4double meany=0.; G4double sigmay=1.*mm;
  x0 = G4RandGauss::shoot(meanx,sigmax);
  y0 = G4RandGauss::shoot(meany,sigmay);
  z0 = -3.0*m;

  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  particleGun->GeneratePrimaryVertex(anEvent);
}


