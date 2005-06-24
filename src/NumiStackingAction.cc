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
#include "NumiTrajectory.hh"
#include "NumiAnalysis.hh"

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
  
  //Discard particles with kinetic energy < 1.GeV (that are not neutrinos)
  if ((particleType!=G4NeutrinoE::NeutrinoEDefinition())&&
      (particleType!=G4NeutrinoMu::NeutrinoMuDefinition())&&
      (particleType!=G4NeutrinoTau::NeutrinoTauDefinition())&&
      (particleType!=G4AntiNeutrinoE::AntiNeutrinoEDefinition())&&
      (particleType!=G4AntiNeutrinoMu::AntiNeutrinoMuDefinition())&&
      (particleType!=G4AntiNeutrinoTau::AntiNeutrinoTauDefinition()))
    {
      G4double energy = aTrack->GetKineticEnergy();
      if (energy < 1.*GeV&&(classification != fKill)) 
	{classification = fKill;} 
    }
  else {
    NumiAnalysis *NA=NumiAnalysis::getInstance();
    NumiTrajectory* NuParentTrack=NA->GetParentTrajectory(aTrack->GetParentID());
    G4int noPoints=NuParentTrack->GetPointEntries();
    G4ThreeVector ParentMomentumFinal=NuParentTrack->GetMomentum(noPoints-1);
    if (ParentMomentumFinal*ParentMomentumFinal<1.*GeV){
      classification = fKill;
      // this kills all of the neutrinos that come from parents that are stopped
      // probably it would be better to kill all particles below some energy during tracking
      // that would speed up things 
    }
}
      
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


