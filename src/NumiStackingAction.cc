//
// NumiStackingAction.cc
//

#include "NumiStackingAction.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ios.hh"

NumiStackingAction::NumiStackingAction()
{ 
}

NumiStackingAction::~NumiStackingAction()
{
}

G4ClassificationOfNewTrack 
NumiStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  G4ClassificationOfNewTrack classification = fUrgent;
  G4ParticleDefinition * particleType = aTrack->GetDefinition();

  // Discard Gammas, Electrons, ...
  if ((particleType==G4Gamma::GammaDefinition())||
      (particleType==G4Electron::ElectronDefinition())||
      (particleType==G4Positron::PositronDefinition()))
    {classification = fKill;}

  //Discard particles with pz<0
  G4ThreeVector momentum=aTrack->GetMomentumDirection();
  if (momentum[2]<0) 
    {classification = fKill;}
  
  //Discard particles with kinetic energy < 0.5GeV
  G4double energy = aTrack->GetKineticEnergy();
  if (energy < 0.5*GeV) 
    {classification = fKill;}
  
  return classification;
}

void NumiStackingAction::NewStage()
{
  // stackManager->ReClassify();
  //  return;
}
    
void NumiStackingAction::PrepareNewEvent()
{ 
}


