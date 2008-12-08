//
// NumiRunAction.cc
//

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "NumiRunAction.hh"
#include "NumiRunActionMessenger.hh"
#include "NumiAnalysis.hh"
#include "NumiTrajectory.hh"
#include "NumiDataInput.hh"
#include "Randomize.hh"
#include "NumiRunManager.hh"

NumiRunAction::NumiRunAction()
{ 
  runMessenger = new NumiRunActionMessenger(this);
}


NumiRunAction::~NumiRunAction()
{ 
  delete runMessenger;
}


void NumiRunAction::BeginOfRunAction(const G4Run* aRun)
{
  NumiDataInput *ND = NumiDataInput::GetNumiDataInput();
  NumiRunManager *pRunManager = (NumiRunManager*)NumiRunManager::GetRunManager();

  G4cout << "Starting run " << aRun->GetRunID()<<G4endl;
  G4cout<<"Random seed used for this run "<<HepRandom::getTheSeed();
  G4String randomFile="rndm/beginOfRun_";
  char runN[4];
  sprintf(runN,"%04d",aRun->GetRunID());
  randomFile.append(runN);
  randomFile.append(".rndm");
  HepRandom::saveEngineStatus(randomFile);
  G4cout << "; Random engine status saved in "<<randomFile<<G4endl;

  if (ND->useFlukaInput){
    G4cout<<"Using Fluka input ntuple "<<ND->GetExtNtupleFileName()<<G4endl;
  }
  else if (ND->useMarsInput){
    G4cout<<"Using Mars input ntuple "<<ND->GetExtNtupleFileName()<<G4endl;
  }
  else if(ND->useMuonBeam  && !(ND->useMuonInput)){
    G4cout<<" *** Using Muon beam:"<<G4endl;
  }
  else if(ND->useMuonBeam && ND->useMuonInput){
    G4cout<<" *** Using muon input ntuple " <<ND->GetExtNtupleFileName()<<G4endl;
  }
  else if(ND->useTestBeam){
      G4cout<<" *** Using Fluka-like Gun:"<<G4endl;
  }
  else{
    G4cout<<"Proton beam:"<<G4endl;
      G4cout << " Momentum: "<<ND->protonMomentum/GeV << "GeV" <<G4endl;
      G4cout << " Kinetic Energy: "<<ND->protonKineticEnergy/GeV<<"GeV" <<G4endl;
      G4cout << " Position "<<ND->beamPosition/m<<" m"<<G4endl;
      G4cout << " SigmaX = "<<ND->beamSigmaX/mm<<" mm"<<G4endl;
      G4cout << " SigmaY = "<<ND->beamSigmaY/mm<<" mm"<<G4endl;
  }

  // Outputs whether the rock density and alcove walls have been
  // changed, and by what fraction of the sigma values.
  if ( ND->GetMaterialSigma() != 0 ){
    G4cout << " Material Sigma = " << ND->GetMaterialSigma() << G4endl;
    if ( ND->GetMaterialSigma() < 0 ) ND->SetGeometryTag("_Y");
    else if ( ND->GetMaterialSigma() > 0 ) ND->SetGeometryTag("_V");
  }

  G4cout << "Processing "<<pRunManager->GetNumberOfEvents()<<" particles"<<G4endl;
  //Book histograms and ntuples
  NumiAnalysis* analysis = NumiAnalysis::getInstance();
  analysis->book();
 
}


void NumiRunAction::EndOfRunAction(const G4Run* aRun)
{

  G4String randomFile="rndm/endOfRun_";
  char runN[4];
  sprintf(runN,"%04d",aRun->GetRunID());
  randomFile.append(runN);
  randomFile.append(".rndm");
  HepRandom::saveEngineStatus(randomFile);
  G4cout << "Random engine status at the end of the run saved in "<<randomFile<<G4endl;
  NumiAnalysis* analysis = NumiAnalysis::getInstance();
  analysis->finish();
}




















