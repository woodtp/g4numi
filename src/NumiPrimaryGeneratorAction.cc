
#include "NumiPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

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
  particleGun->GeneratePrimaryVertex(anEvent);
}


