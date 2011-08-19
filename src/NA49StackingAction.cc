
#include "NA49StackingAction.hh"

#include "NA49StackingMessenger.hh"
#include "NumiAnalysis.hh"
#include "G4Track.hh"
#include "G4HadronicProcessStore.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NA49StackingAction::NA49StackingAction()
{
  stackMessenger = new NA49StackingMessenger(this);
  killSecondary  = false;
  pname          = ""; 
  elm            = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NA49StackingAction::~NA49StackingAction()
{
  delete stackMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
NA49StackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
  const G4ParticleDefinition* pd = aTrack->GetDefinition();
  G4HadronicProcessStore* store=G4HadronicProcessStore::Instance();
  G4ClassificationOfNewTrack status = fUrgent;
  NA49Analysis* analysis = NA49Analysis::getInstance();
  primaryDef = pd;
  primaryTotalEnergy = aTrack->GetTotalEnergy()/GeV;
  // G4double xs=store->GetInelasticCrossSectionPerAtom(primaryDef,primaryTotalEnergy,elm)/barn*1000;//mb units  
  if (aTrack->GetTrackStatus() == fAlive) 
    if(1 == aTrack->GetParentID())analysis->FillNtuple(*aTrack); 

  const G4String name = aTrack->GetDefinition()->GetParticleName();
  //  if(aTrack->GetTrackID() == 1) return status;

  //stack or delete secondaries
  if (killSecondary)      status = fKill;
  else if(pname == name)  status = fKill; 

  return status;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
