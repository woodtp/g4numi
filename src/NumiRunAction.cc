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
   NumiDataInput *NumiData = NumiDataInput::GetNumiDataInput();
   if(NumiData->fPrintInfo > 0 || NumiData->IsDebugOn())
   {
      G4cout << "NumiRunAction Constructor Called." << G4endl;
   }
   
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

  G4cout << G4endl;
  G4cout << G4endl;
  G4cout << "********************************************************************" << G4endl;
  G4cout << "********************************************************************" << G4endl;
  G4cout << "NumiRunAction::BeginOfRunAction - Starting Run..." << G4endl;
  G4cout << "********************************************************************" << G4endl;

  G4cout << "  Starting run " << aRun->GetRunID() << G4endl;
  G4cout << "  Random seed used for this run " << CLHEP::HepRandom::getTheSeed();
  G4String randomFile="rndm/beginOfRun_";
  char runN[4];
  sprintf(runN,"%04d",aRun->GetRunID());
  randomFile.append(runN);
  randomFile.append(".rndm");
  CLHEP::HepRandom::saveEngineStatus(randomFile);
  G4cout << "; Random engine status saved in "<<randomFile<<G4endl;

  if (ND->useFlukaInput){
    G4cout<<"  Using Fluka input ntuple "<<ND->GetExtNtupleFileName()<<G4endl;
  }
  else if (ND->useMarsInput){
    G4cout<<"  Using Mars input ntuple "<<ND->GetExtNtupleFileName()<<G4endl;
  }
  else if(ND->useMuonBeam  && !(ND->useMuonInput)){
    G4cout<<"  Using Muon beam:"<<G4endl;
  }
  else if(ND->useMuonBeam && ND->useMuonInput){
    G4cout<<"  Using muon input ntuple " <<ND->GetExtNtupleFileName()<<G4endl;
  }
  else if(ND->useTestBeam){
      G4cout<<"  Using Fluka-like Gun:"<<G4endl;
  }
  else if(ND->useMacro){
    G4cout<<"  Following the Macro set Parameters"<<G4endl;
  }
  else if(ND->GetDetailedProtonBeam())
  {
     G4cout << "  Using the Detailed Proton Beam." << G4endl;  	 
  }
  else
  {
     G4cout << "  Using the Simple Proton beam:"<<G4endl;
     G4cout << "     Momentum: "<<ND->protonMomentum/GeV << "GeV" <<G4endl;
     G4cout << "     Kinetic Energy: "<<ND->protonKineticEnergy/GeV<<"GeV" <<G4endl;
     G4cout << "     Position "<<ND->beamPosition/m<<" m"<<G4endl;
     G4cout << "     SigmaX = "<<ND->beamSigmaX/mm<<" mm"<<G4endl;
     G4cout << "     SigmaY = "<<ND->beamSigmaY/mm<<" mm"<<G4endl;
  }


  G4cout << G4endl;
  G4cout << G4endl;
  ND -> Print();
  G4cout << G4endl;
  G4cout << G4endl;

  //---------------------------------------
  //Check that the user has set up the simulation the way he/she really wants it.
  //If he/she hasn't, don't let him/her run.
  //

  if((ND -> GetSimulation()).empty())
  {
     G4cout << "  PROBLEM: Undefined Simulation. Can't Run. grep for \"PROBLEM\" "
            << " in the output of this simulation to debug the problem." << G4endl;
     ND -> SetOkToRun(false);
  }
  else if( (ND -> GetSimulation()) == "Neutrino Beam Simulation" &&
           (ND -> GetBeamConfig()).empty())
  {
        G4cout << "  PROBLEM: Undefined Beam Configuration. Can't Run. grep for \"PROBLEM\" "
               << " in the output of this simulation to debug the problem." << G4endl;
        ND -> SetOkToRun(false);
  }
  else
  {

     ND -> SetOkToRun(true);
  }
  //---------------------------------------

  
  //------------------------
  //this is bad and may be obsolete. Fix this. (Laura)
  //
  
  // Outputs whether the rock density and alcove walls have been
  // changed, and by what fraction of the sigma values.
  if ( ND->GetMaterialSigma() != 0 ){
     G4cout << " Material Sigma = " << ND->GetMaterialSigma() << G4endl;
     if ( ND->GetMaterialSigma() < 0 ) ND->SetGeometryTag("_Y");
     else if ( ND->GetMaterialSigma() > 0 ) ND->SetGeometryTag("_V");
  }
  //------------------------

  G4cout << G4endl;
  G4cout << G4endl;
  G4cout << "  Processing "<<pRunManager->GetNumberOfEvents()<<" particles"<<G4endl;
  //Book histograms and ntuples
  NumiAnalysis* analysis = NumiAnalysis::getInstance();
  analysis->book();

  
  G4cout << "********************************************************************" << G4endl;
  G4cout << "NumiRunAction::BeginOfRunAction - Completed." << G4endl;
  G4cout << "********************************************************************" << G4endl;
  G4cout << "********************************************************************" << G4endl;
  G4cout << G4endl;
  G4cout << G4endl;

}


void NumiRunAction::EndOfRunAction(const G4Run* aRun)
{

   G4cout << G4endl;
  G4cout << G4endl;
  G4cout << "********************************************************************" << G4endl;
  G4cout << "********************************************************************" << G4endl;
  G4cout << "NumiRunAction::EndOfRunAction..." << G4endl;
  G4cout << "********************************************************************" << G4endl;

  G4String randomFile="rndm/endOfRun_";
  char runN[4];
  sprintf(runN,"%04d",aRun->GetRunID());
  randomFile.append(runN);
  randomFile.append(".rndm");
  CLHEP::HepRandom::saveEngineStatus(randomFile);
  G4cout << "  Random engine status at the end of the run saved in "<<randomFile<<G4endl;
  NumiAnalysis* analysis = NumiAnalysis::getInstance();
  analysis->finish();

  G4cout << "********************************************************************" << G4endl;
  G4cout << "NumiRunAction::EndOfRunAction - Completed." << G4endl;
  G4cout << "********************************************************************" << G4endl;
  G4cout << "********************************************************************" << G4endl;
  G4cout << G4endl;
  G4cout << G4endl;
}




















