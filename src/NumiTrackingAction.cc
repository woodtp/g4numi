// 
// NumiTrackingAction.cc
//
#include "NumiTrackInformation.hh"
#include "NumiTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4Trajectory.hh"
#include "G4RunManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "NumiAnalysis.hh"
#include "NumiTrajectory.hh"


NumiTrackingAction::NumiTrackingAction()
{;}

NumiTrackingAction::~NumiTrackingAction()
{;}

void NumiTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  //set tgen=1 for primary protons 
  if (aTrack->GetTrackID()==1) 
    {   
      G4RunManager* pRunManager=G4RunManager::GetRunManager();
      G4int p_no=pRunManager->GetCurrentEvent()->GetEventID();
      if (p_no%500==0) G4cout<<"Processing proton #: "<<p_no<<" to "<< p_no+500<<G4endl;
      NumiTrackInformation* info=new NumiTrackInformation();
      info->Settgen(1); 
      G4Track* theTrack = (G4Track*)aTrack;
      theTrack->SetUserInformation(info);
    }
 
   //Store tracks in trajectory container
   fpTrackingManager->SetStoreTrajectory(true);
   fpTrackingManager->SetTrajectory(new NumiTrajectory(aTrack));
   
   //if a particle is a neutrino then analyse and store in ntuple
   G4ParticleDefinition * particleType = aTrack->GetDefinition();
   if ((particleType==G4NeutrinoE::NeutrinoEDefinition())||
       (particleType==G4NeutrinoMu::NeutrinoMuDefinition()) ||
       (particleType==G4NeutrinoTau::NeutrinoTauDefinition()) ||
       (particleType==G4AntiNeutrinoE::AntiNeutrinoEDefinition()) ||
       (particleType==G4AntiNeutrinoMu::AntiNeutrinoMuDefinition()) ||
       (particleType==G4AntiNeutrinoTau::AntiNeutrinoTauDefinition()))
    {
      NumiAnalysis* analysis = NumiAnalysis::getInstance();
      analysis->analyseStepping(*aTrack);
    }
  
}

void NumiTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  // set tgen(secondary)=tgen(parent)+1
  NumiTrackInformation* info=(NumiTrackInformation*)(aTrack->GetUserInformation());
    if (info!=0) {
       G4int tgen=info->Gettgen();
       G4TrackVector* SecVector=fpTrackingManager->GimmeSecondaries();

       for (size_t ii=0;ii<(*SecVector).size();ii++){
	 NumiTrackInformation* newinfo=new NumiTrackInformation();
	 newinfo->Settgen(tgen+1);
	 (*SecVector)[ii]->SetUserInformation(newinfo);
       }
    }
}






