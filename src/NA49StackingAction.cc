
#include "NA49StackingAction.hh"
#include "NA49Analysis.hh"
#include "G4Track.hh"
#include "G4HadronicProcessStore.hh"
#include "G4NistManager.hh"
#include "NA49EventAction.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NA49StackingAction::NA49StackingAction()
{
  killSecondary  = false;
  pname          = ""; 
  elm            = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NA49StackingAction::~NA49StackingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
NA49StackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
 G4ClassificationOfNewTrack status = fUrgent;
 EvtManager = G4EventManager::GetEventManager();
 NA49EvtAct = (NA49EventAction*)(EvtManager -> GetUserEventAction());
 const G4ParticleDefinition* pd = aTrack->GetDefinition();
 if (aTrack->GetTrackStatus() == fAlive) 
   if( (1 == aTrack->GetParentID()) && (aTrack->GetCreatorProcess()->GetProcessName()=="ProtonInelastic") )
    {NA49EvtAct->AddTrack(aTrack);}

  //stack or delete secondaries
 if (killSecondary)      status = fKill;

  return status;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
