#include "NumiPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "NumiDataInput.hh"
#include "G4UImanager.hh"
#include "NumiRunManager.hh"
#include "NumiParticleCode.hh"

#include "TROOT.h"
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TSystem.h>


NumiPrimaryGeneratorAction::NumiPrimaryGeneratorAction()
{
  fND=NumiDataInput::GetNumiDataInput();
  G4int n_particle = 1;
  fIsFirst=true;
  fParticleGun = new G4ParticleGun(n_particle);
  fRunManager=(NumiRunManager*)NumiRunManager::GetRunManager();
  tunnelPos = G4ThreeVector(0,0,fND->TunnelLength/2.+fND->TunnelZ0);

}

NumiPrimaryGeneratorAction::~NumiPrimaryGeneratorAction()
{
  delete fParticleGun;
}

void NumiPrimaryGeneratorAction::SetProtonBeam()
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
 
  fParticleGun->SetParticleDefinition(particleTable->FindParticle("proton"));
  fParticleGun->SetParticleEnergy(fND->protonKineticEnergy);
  fParticleGun->SetParticlePosition(fND->beamPosition);
  fParticleGun->SetParticleMomentumDirection(fND->beamDirection);

  fCurrentPrimaryNo=0;
}

// Used as a quicker check of muon response
// past the end of the decay pipe. The values are hard-coded
// for the muon 'beam' as to make it easier to change.

void NumiPrimaryGeneratorAction::SetMuonBeam()
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4float xpos, ypos, zpos;

  xpos = 0;
  ypos = 0;
  zpos = 700*m;
 
  G4ThreeVector muonBeamPos = G4ThreeVector(xpos, ypos, zpos) - tunnelPos;
  G4ThreeVector muonBeamDirection = G4ThreeVector(0, 0, 1); // straight down the pipe

  fParticleGun->SetParticleDefinition(particleTable->FindParticle("mu-"));
  fParticleGun->SetParticleEnergy(80*GeV);// just a guess right now
  fParticleGun->SetParticlePosition(muonBeamPos);
  fParticleGun->SetParticleMomentumDirection(muonBeamDirection);

  fCurrentPrimaryNo=0;
}

G4bool NumiPrimaryGeneratorAction::OpenNtuple(G4String ntupleName)
{
  G4bool fIsOpen=false;
  fRootFile=new TFile(ntupleName,"READ");
  if (!fRootFile->IsZombie())
    {
      fPrimaryNtuple=(TTree*)(fRootFile->Get("h3"));
      fCurrentPrimaryNo=0;
      fNoOfPrimaries=fPrimaryNtuple->GetEntries();
      fIsOpen=true;
    }
  else
    {
      G4cout<<"NumiPrimaryGeneratorAction: Input (fluka/mars) root file doesn't exist"
	    <<"   Aborting run"<<G4endl;
      fIsOpen=false;
    }
  fCurrentPrimaryNo=0;
  return fIsOpen;
}

void NumiPrimaryGeneratorAction::CloseNtuple()
{
    fRootFile->Close();
    fCurrentPrimaryNo=0;
}

void NumiPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //  std::cout<<"*************** anEvent: "<<anEvent->get_eventID()<<" *******************"<<std::endl;
  //G4UImanager* UI=G4UImanager::GetUIpointer();

  G4int totNoPrim=fRunManager->GetNumberOfEvents();
  if (totNoPrim>20){
    if (fCurrentPrimaryNo%(totNoPrim/20)==0) 
      G4cout<<"Processing particles #: "
	    <<fCurrentPrimaryNo<<" to "<< fCurrentPrimaryNo+(totNoPrim/20)<<G4endl;
  }
  
  if (!fND->useFlukaInput && !fND->useMarsInput && !fND->useMuonBeam){
    G4double x0;
    G4double y0; 
    G4double z0;
    G4double sigmax=fND->beamSigmaX;
    G4double sigmay=fND->beamSigmaY;
    
    x0 = G4RandGauss::shoot(fND->beamPosition[0],sigmax);
    y0 = G4RandGauss::shoot(fND->beamPosition[1],sigmay);
    z0 = fND->beamPosition[2];
    
    fProtonOrigin   = G4ThreeVector(x0, y0, z0);
    fProtonMomentum = G4ThreeVector(0, 0, fND->protonMomentum);
    fProtonIntVertex = G4ThreeVector(-9999.,-9999.,-9999.);
    fWeight=1.; //for primary protons set weight and tgen
    fTgen=0;
    
    fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
  else if(fND->useMuonBeam){
  
    G4double x0 = fND->DecayPipeRadius*2.0;
    G4double y0 = fND->DecayPipeRadius*2.0; 
    G4double z0;

    while(sqrt(pow(x0,2)+pow(y0,2)) > fND->DecayPipeRadius){
      x0 = 2*(G4UniformRand()-0.5)*fND->DecayPipeRadius;
      y0 = 2*(G4UniformRand()-0.5)*fND->DecayPipeRadius;
    }
      z0 = 722*m;// just before the endcap

    fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }


  if (fND->useFlukaInput || fND->useMarsInput) {
    /*
      Fluka and Mars input variables:
      FLUKA                           MARS                   
      -----------------------------------------------------------------------------------------------
      TYPE                            TYPE                                      - type of particle (see NumiAnalysis::GetParticleName())
      X, Y, Z                         X,Y,Z                                     - particle coordinates
      PX, PY, PZ                      PX, PY, PZ                                - particle momentum
      WEIGHT                          WEIGHT                                    - particle weight
      GENER                           GENER                                     - particle generation
      PROTVX, PROTVY, PROTVZ          PVTXX, PVTXY, PVTXZ                       - proton interaction vertex 
      PROTX, PROTY, PROTZ             PROTX, PROTY, PROTZ                       - proton initial coordinates
      PROTPX, PROTPY, PROTPZ          PROTPX, PROTPY, PROTPZ                    - proton initial momentum
      MOMPX,MOMPY,MOMPZ               PPX, PPY, PPZ                             - ???
      MOMTYPE                         PTYPE                                     - ???
    */
    
    G4double x0,y0,z0,px,py,pz;
    G4String particleName;
    fPrimaryNtuple->GetEntry(fCurrentPrimaryNo);
    
    x0 = fPrimaryNtuple->GetLeaf("x")->GetValue()*cm;
    y0 = fPrimaryNtuple->GetLeaf("y")->GetValue()*cm;
    z0 = fPrimaryNtuple->GetLeaf("z")->GetValue()*cm+fND->TargetZ0+35*cm;
    px = fPrimaryNtuple->GetLeaf("px")->GetValue()*GeV;
    py = fPrimaryNtuple->GetLeaf("py")->GetValue()*GeV;
    pz = fPrimaryNtuple->GetLeaf("pz")->GetValue()*GeV;
    
    fParticlePosition=G4ThreeVector(x0,y0,z0);
    fParticleMomentum=G4ThreeVector(px,py,pz);
    
    fWeight = fPrimaryNtuple->GetLeaf("weight")->GetValue();
    fType = G4int(fPrimaryNtuple->GetLeaf("type")->GetValue());
    fTgen = G4int(fPrimaryNtuple->GetLeaf("gener")->GetValue());
    particleName=NumiParticleCode::AsString(NumiParticleCode::IntToEnum(fType));
    fProtonOrigin   = G4ThreeVector(fPrimaryNtuple->GetLeaf("protx")->GetValue()*cm,
				    fPrimaryNtuple->GetLeaf("proty")->GetValue()*cm,
				    fPrimaryNtuple->GetLeaf("protz")->GetValue()*cm);
    fProtonMomentum = G4ThreeVector(fPrimaryNtuple->GetLeaf("protpx")->GetValue()*cm,
				    fPrimaryNtuple->GetLeaf("protpy")->GetValue()*cm,
				    fPrimaryNtuple->GetLeaf("protpz")->GetValue()*cm);
    
    if (fND->useMarsInput){
      fProtonIntVertex = G4ThreeVector(fPrimaryNtuple->GetLeaf("pvtxx")->GetValue()*cm,
				       fPrimaryNtuple->GetLeaf("pvtxy")->GetValue()*cm,
				       fPrimaryNtuple->GetLeaf("pvtxz")->GetValue()*cm);
    }
    else if (fND->useFlukaInput){
      fProtonIntVertex = G4ThreeVector(fPrimaryNtuple->GetLeaf("protvx")->GetValue()*cm,
				       fPrimaryNtuple->GetLeaf("protvy")->GetValue()*cm,
				       fPrimaryNtuple->GetLeaf("protvz")->GetValue()*cm);
    }
    
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    fParticleGun->SetParticleDefinition(particleTable->FindParticle(particleName));
    G4double mass=particleTable->FindParticle(particleName)->GetPDGMass();
    
    fParticleGun->SetParticleEnergy(sqrt(mass*mass+px*px+py*py+pz*pz)-mass);
    fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
    fParticleGun->SetParticleMomentum(G4ThreeVector(px,py,pz));
    fParticleGun->GeneratePrimaryVertex(anEvent);
    
  }
  fCurrentPrimaryNo++;
}


