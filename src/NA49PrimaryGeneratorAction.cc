
#include "NA49PrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "Randomize.hh"
#include "NA49Analysis.hh"
#include <time.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NA49PrimaryGeneratorAction::NA49PrimaryGeneratorAction()
{
  particleGun  = new G4ParticleGun(1);
  /*
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  G4long seed=time(0);
  CLHEP::HepRandom::setTheSeed(seed);
  //randFlat=new CLHEP::RandFlat(&ranecuEngine);
  //srand(time(0));

  G4double theta=CLHEP::RandFlat::shoot()*2*3.14159265;
  G4double pt=CLHEP::RandFlat::shoot()*2.0;//max pT is 2 GeV
  G4double px0=pt*cos(theta);
  G4double py0=pt*sin(theta);
  G4double pz0=sqrt(158.*158.-(px0*px0+py0*py0));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(px0,py0,pz0));
  */
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NA49PrimaryGeneratorAction::~NA49PrimaryGeneratorAction()
{
  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NA49PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4double zVertex = -50.*mm;
  particleGun->SetParticlePosition(G4ThreeVector(0.,0.,zVertex));
  particleGun->GeneratePrimaryVertex(anEvent);
  NA49Analysis* analysis = NA49Analysis::getInstance();
  Double_t ener = particleGun->GetParticleEnergy();
  G4ParticleDefinition* part  = particleGun->GetParticleDefinition();

  analysis->GetPrimGenInfo(ener,part);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
