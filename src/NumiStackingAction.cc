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
#include "NumiImpWeight.hh"
#include "NumiTrackInformation.hh"
#include "NumiDataInput.hh"

NumiStackingAction::NumiStackingAction()
{ 
  NumiData=NumiDataInput::GetNumiDataInput();
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
      (particleType==G4Positron::PositronDefinition())&&
      (classification != fKill))
    {classification = fKill;}

  //Discard particles with pz<0
  G4ThreeVector momentum=aTrack->GetMomentumDirection();
  if (momentum[2]<0&&(classification != fKill)) 
    {classification = fKill;}
  
  //Discard particles with kinetic energy < 0.5GeV
  G4double energy = aTrack->GetKineticEnergy();
  if (energy < 0.5*GeV&&(classification != fKill)) 
    {classification = fKill;} 
  
  //If importance weighting is on:
  if (NumiData->NImpWeightOn&&(classification != fKill)){
    NumiTrackInformation* oldinfo=(NumiTrackInformation*)(aTrack->GetUserInformation());  
    if (oldinfo!=0) {
      if (oldinfo->GetNImpWt()==1.){                                     //Check if it already has weight (because one of the parents was weighted)
	G4double Nimpweight=NumiImpWeight::CalculateImpWeight(aTrack);

	if(Nimpweight==0)
	  {classification = fKill;} 
	else {	  
	  NumiTrackInformation* oldinfo=(NumiTrackInformation*)(aTrack->GetUserInformation());  
	  if (oldinfo!=0) oldinfo->SetNImpWt(Nimpweight);
	  // only primary protons don't have TrackInformation already     
	  // all others have some info set by NumiTrackingAction
	}
      }
    }
  }

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


