//
// NumiRunAction.cc
//

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "NumiRunAction.hh"
#include "NumiAnalysis.hh"
#include "NumiTrajectory.hh"
#include "Randomize.hh"

NumiRunAction::NumiRunAction()
{ }


NumiRunAction::~NumiRunAction()
{ }


void NumiRunAction::BeginOfRunAction(const G4Run* aRun)
{
  HepRandom::restoreEngineStatus("random.rndm");
  HepRandom::showEngineStatus();
  //Book histograms and ntuples
  NumiAnalysis* analysis = NumiAnalysis::getInstance();
  analysis->book();
 
}


void NumiRunAction::EndOfRunAction(const G4Run* )
{
  HepRandom::showEngineStatus();
  HepRandom::saveEngineStatus("random.rndm");
  NumiAnalysis* analysis = NumiAnalysis::getInstance();
  analysis->finish();
}
























