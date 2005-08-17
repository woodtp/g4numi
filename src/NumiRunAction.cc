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
#include "Randomize.hh"

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
  G4cout << "Starting run " << aRun->GetRunID()<<G4endl;
  G4cout<<"Random seed used for this run "<<HepRandom::getTheSeed();
  G4String randomFile="rndm/beginOfRun_";
  char runN[4];
  sprintf(runN,"%04d",aRun->GetRunID());
  randomFile.append(runN);
  randomFile.append(".rndm");
  HepRandom::saveEngineStatus(randomFile);
  G4cout << "; Random engine status saved in "<<randomFile<<G4endl;
  
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




















