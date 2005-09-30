#include "NumiPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "NumiDataInput.hh"
#include "G4UImanager.hh"
#include "NumiRunManager.hh"
#include "NumiAnalysis.hh"

#include "TROOT.h"
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TSystem.h>


NumiPrimaryGeneratorAction::NumiPrimaryGeneratorAction()
{
  ND=NumiDataInput::GetNumiDataInput();
  
  G4int n_particle = 1;
  isFirst=true;
  particleGun = new G4ParticleGun(n_particle);
  numiAnalysis=NumiAnalysis::getInstance();
  pRunManager=(NumiRunManager*)NumiRunManager::GetRunManager();
}

NumiPrimaryGeneratorAction::~NumiPrimaryGeneratorAction()
{
  delete particleGun;
}

void NumiPrimaryGeneratorAction::SetProtonBeam()
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
 
  particleGun->SetParticleDefinition(particleTable->FindParticle("proton"));
  particleGun->SetParticleEnergy(ND->protonKineticEnergy);
  particleGun->SetParticlePosition(ND->beamPosition);
  particleGun->SetParticleMomentumDirection(ND->beamDirection);

  currentPrimaryNo=0;
}

G4bool NumiPrimaryGeneratorAction::OpenNtuple(G4String ntupleName)
{
  G4bool isOpen=false;
  rootFile=new TFile(ntupleName,"READ");
  if (!rootFile->IsZombie())
    {
      primaryNtuple=(TTree*)(rootFile->Get("h1"));
      currentPrimaryNo=0;
      noOfPrimaries=primaryNtuple->GetEntries();
      isOpen=true;
    }
  else
    {
      G4cout<<"NumiPrimaryGeneratorAction: Input (fluka/mars) root file doesn't exist"
	    <<"   Aborting run"<<G4endl;
      isOpen=false;
    }
  currentPrimaryNo=0;
  return isOpen;
}

void NumiPrimaryGeneratorAction::CloseNtuple()
{
    rootFile->Close();
    currentPrimaryNo=0;
}

void NumiPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  
  //G4UImanager* UI=G4UImanager::GetUIpointer();
  
    G4int totNoPrim=pRunManager->GetNumberOfEvents();
    //  G4cout <<currentPrimaryNo<<G4endl;
    if (pRunManager==NULL) G4cout<<"RunManager==NULL??"<<G4endl;
    if (totNoPrim>20){
      if (currentPrimaryNo%(totNoPrim/20)==0) 
	G4cout<<"Processing particles #: "
	      <<currentPrimaryNo<<" to "<< currentPrimaryNo+(totNoPrim/20)<<G4endl;
    }

    if (!ND->useFlukaInput&&!ND->useMarsInput){
    G4double x0;
    G4double y0; 
    G4double z0;
    G4double sigmax=ND->beamSigmaX;
    G4double sigmay=ND->beamSigmaY;
  
    x0 = G4RandGauss::shoot(ND->beamPosition[0],sigmax);
    y0 = G4RandGauss::shoot(ND->beamPosition[1],sigmay);
    z0 = ND->beamPosition[2];
    
    particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
    particleGun->GeneratePrimaryVertex(anEvent);
  }
  if (ND->useFlukaInput||ND->useMarsInput) {
    G4double x0,y0,z0,px,py,pz;
    G4String particleName;
    primaryNtuple->GetEntry(currentPrimaryNo);
  
    x0 = primaryNtuple->GetLeaf("x")->GetValue()*cm+ND->TargetX0;
    y0 = primaryNtuple->GetLeaf("y")->GetValue()*cm+ND->TargetY0;
    z0 = primaryNtuple->GetLeaf("z")->GetValue()*cm+ND->TargetZ0+35*cm;
    px = primaryNtuple->GetLeaf("px")->GetValue()*GeV;
    py = primaryNtuple->GetLeaf("py")->GetValue()*GeV;
    pz = primaryNtuple->GetLeaf("pz")->GetValue()*GeV;

    ParticlePosition=G4ThreeVector(x0,y0,z0);
    ParticleMomentum=G4ThreeVector(px,py,pz);

    weight = primaryNtuple->GetLeaf("weight")->GetValue();
    type = G4int(primaryNtuple->GetLeaf("type")->GetValue());
    tgen = G4int(primaryNtuple->GetLeaf("gener")->GetValue());
    particleName=numiAnalysis->GetParticleName(type);
  
    if (ND->useMarsInput){
      protonX=primaryNtuple->GetLeaf("protx")->GetValue()*cm;
      protonY=primaryNtuple->GetLeaf("proty")->GetValue()*cm;
      //      protonDirection=G4ThreeVector(primaryNtuple->GetLeaf("protdxdz")->GetValue(),
      //			      primaryNtuple->GetLeaf("protdydz")->GetValue(),
      //			      1.);
    }
 
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    particleGun->SetParticleDefinition(particleTable->FindParticle(particleName));
    G4double mass=particleTable->FindParticle(particleName)->GetPDGMass();
    particleGun->SetParticleEnergy(sqrt(mass*mass+px*px+py*py+pz*pz)-mass);
    particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
    particleGun->SetParticleMomentum(G4ThreeVector(px,py,pz));
    particleGun->GeneratePrimaryVertex(anEvent);
    
  }
  currentPrimaryNo++;
}


