//----------------------------------------------------------------------
// $Id
//----------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <iomanip>

#include "NumiPrimaryGeneratorAction.hh"

#include <math.h>
#include <map>
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "NumiPrimaryMessenger.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "NumiDataInput.hh"
#include "G4UImanager.hh"
#include "NumiRunManager.hh"
#include "NumiParticleCode.hh"
#include "NumiAnalysis.hh"

#include "TROOT.h"
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TSystem.h>

using namespace std;

NumiPrimaryGeneratorAction::NumiPrimaryGeneratorAction()
{
  beamMessenger = new NumiPrimaryMessenger(this);
  fND=NumiDataInput::GetNumiDataInput();
  G4int n_particle = 1;
  fIsFirst=true;
  fParticleGun = new G4ParticleGun(n_particle);
  fRunManager=(NumiRunManager*)NumiRunManager::GetRunManager();
  tunnelPos = G4ThreeVector(0,0,fND->TunnelLength/2.+fND->TunnelZ0);
  cosx=0;
  cosy=0;

}

NumiPrimaryGeneratorAction::~NumiPrimaryGeneratorAction()
{
  delete beamMessenger;
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
  zpos = fND->DecayPipeLength + fND->TunnelZ0;;
 
  G4ThreeVector muonBeamPos = G4ThreeVector(xpos, ypos, zpos) - tunnelPos;
  G4ThreeVector muonBeamDirection = G4ThreeVector(0, 0, 1); // straight down the pipe

  fParticleGun->SetParticleDefinition(particleTable->FindParticle("mu-"));
  //  fParticleGun->SetParticleEnergy(80*GeV);// just a guess right now
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

      if(fND->useMuonInput && fND->useMuonBeam)
      {
	fPrimaryNtuple=(TTree*)(fRootFile->Get("muon"));
	fMuon = new NtpMuon();
	fMuon -> SetTree(fPrimaryNtuple);

      }
      else
      {
	fPrimaryNtuple=(TTree*)(fRootFile->Get("h3"));
      }
      if(!fPrimaryNtuple)
      {
	G4cout<<"NumiPrimaryGeneratorAction: Can't find tree "<< G4endl;
      }

      fCurrentPrimaryNo=0;
      fNoOfPrimaries=fPrimaryNtuple->GetEntries();

      //G4cout << "*****************entries = " <<fPrimaryNtuple->GetEntries() << G4endl;
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
  
  if(fND->useMuonBeam && !(fND->useMuonInput))
  {
 
    G4double x0 = fND->DecayPipeRadius*2.0;
    G4double y0 = fND->DecayPipeRadius*2.0; 
    G4double z0;
    
    // Uniformly distributed circular muon beam 
    while(sqrt(pow(x0,2)+pow(y0,2)) > fND->DecayPipeRadius){
      x0 = 2*(G4UniformRand()-0.5)*fND->DecayPipeRadius;
      y0 = 2*(G4UniformRand()-0.5)*fND->DecayPipeRadius;
    }
    
    //Square beam the fits just inside the decay pipe without
    // producing particles that would immediatly hit rock outside the
    // decay pipe.
    /*
     x0 = 2*(G4UniformRand()-0.5)*fND->DecayPipeRadius/sqrt(2.0);
     y0 = 2*(G4UniformRand()-0.5)*fND->DecayPipeRadius/sqrt(2.0);
     */
    
    z0 = fND->DecayPipeLength + fND->TunnelZ0 - 1;// just before the endcap
    
    fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
  else if (fND->useMuonInput && fND->useMuonBeam)
  {
     
    G4String particleName;

    NumiAnalysis* analysis2 = NumiAnalysis::getInstance();
    analysis2->SetCount(0);
    analysis2->SetEntry(fCurrentPrimaryNo);
    analysis2->SetAlcEdepFlag(false);
    
    fMuon -> Clear();
    fMuon -> GetEntry(fCurrentPrimaryNo);
    fevtno = fMuon->evtno;
    fmuweight = fMuon->muweight;

    fMuTParentMomentum = G4ThreeVector(fMuon->tpx*GeV,fMuon->tpy*GeV,fMuon->tpz*GeV);
    fMuTParentPosition = G4ThreeVector(fMuon->tvx*cm,fMuon->tvy*cm,fMuon->tvz*cm);
    fMuParentMomentum = G4ThreeVector(fMuon->pdpx*GeV,fMuon->pdpy*GeV,fMuon->pdpz*GeV);
    fMuParentPosition = G4ThreeVector(fMuon->pdvx*cm,fMuon->pdvy*cm,fMuon->pdvz*cm);
    fMuParentProdPosition = G4ThreeVector(fMuon->ppvx*cm,fMuon->ppvy*cm,fMuon->ppvz*cm);
    ftpptype = fMuon->tptype;
    fnimpwt = fMuon->nimpwt;
    fpptype = fMuon->ptype;
    fppmedium = fMuon->ppmedium;
    fpgen = fMuon->pgen;


    //
    //z0 = fND->DecayPipeLength + fND->TunnelZ0 - 1;// just before the endcap this is what Jason uses.
    //z0 =        676681        + 45698.5       - 1 = 722378.5 mm
    //which means the endcap is at 722379.5 mm
    //fMuon->muvz*cm = 722380, so subtract 1mm to get 722379 just to be safe
    //

    G4double z0 = fMuon->muvz*cm - 1; //this is - 1 mm.
    
    fParticlePosition=G4ThreeVector(fMuon->muvx*cm,fMuon->muvy*cm,z0);
    fParticleMomentum=G4ThreeVector(fMuon->mupx*GeV,fMuon->mupy*GeV,fMuon->mupz*GeV);
       
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    fParticleGun->SetParticleDefinition(particleTable->FindParticle("mu-"));
    G4double mass=particleTable->FindParticle("mu-")->GetPDGMass();

    NumiAnalysis* analysis = NumiAnalysis::getInstance();
    if(!analysis)
    {
      G4cout << "Can't get NumiAnalysis pointer" << G4endl;
    }
    analysis->FillHadmmNtuple();
    
    if(fMuon->mupz < 1.0)
    {
      flocalParticleMomentum = G4ThreeVector(0.0,0.0,0.00*GeV);
      fParticleGun->SetParticleEnergy(sqrt(mass*mass+0.0*0.0*GeV+0.0*0.0*GeV+0.00*GeV*0.00*GeV)-mass);

      analysis->WriteHadmmNtuple();    
      
    }
    else
    {
      flocalParticleMomentum = G4ThreeVector(fMuon->mupx*GeV,fMuon->mupy*GeV,fMuon->mupz*GeV);
      fParticleGun->SetParticleEnergy(sqrt(mass*mass+fMuon->mupx*GeV*fMuon->mupx*GeV+fMuon->mupy*GeV*fMuon->mupy*GeV+fMuon->mupz*GeV*fMuon->mupz*GeV)-mass);

    }
  
    fParticleGun->SetParticlePosition(fParticlePosition);
    fParticleGun->SetParticleMomentum(flocalParticleMomentum);
    fParticleGun->GeneratePrimaryVertex(anEvent);
    
  }
  else if (fND->useTestBeam) { 
    // Test Beam for comparison with FLUGG
    // Fluka-analagous paramaters defined in input file
    
    G4double x0, y0;
    G4double px, py, pz;
    
    if (outR > 0) {
      G4double r0   = G4UniformRand()*(outR - inR) + inR;
      G4double phi0 = G4UniformRand()*2*M_PI;
      
      x0 = cos(phi0) * r0;
      y0 = sin(phi0) * r0;
    }
    else {
      x0 = 0;
      y0 = 0;
    }
    
    if (spread > 0) {
      p += DoubleRand()*spread/2.;
    }
    
    if (cosx + cosy == 0) {
      G4double ptheta  = sqrt(G4UniformRand())*div/2.;
      G4double pphi    = G4UniformRand()*2*M_PI; 
      px = p * sin(ptheta)*cos(pphi);
      py = p * sin(ptheta)*sin(pphi);
      pz = p * cos(ptheta);
    }
    else {
      px = p * cosx;
      py = p * cosy;
      pz = p * sqrt(1 - cosx*cosx - cosy*cosy);
    }
    
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleEnergy(0.0); 
    fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
    fParticleGun->SetParticleMomentum(G4ThreeVector(px, py, pz));
    
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
  else if (fND->useMacro){
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
  else if (fND->useFlukaInput || fND->useMarsInput) {
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
  else {
    // If nothing else is set, use a proton beam
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
  
  fCurrentPrimaryNo++;
}


