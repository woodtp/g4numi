//
// NumiSteppingAction.cc
//

#include "NumiSteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "NumiTrajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4RunManager.hh"
#include "NumiTrackInformation.hh"
#include "NumiAnalysis.hh"

NumiSteppingAction::NumiSteppingAction()
{  
}

NumiSteppingAction::~NumiSteppingAction()
{
}

void NumiSteppingAction::UserSteppingAction(const G4Step * theStep)
{
  // Check if the Pi+, Pi-, K+, K-, K0L, mu+ or mu- decayed and set Ndecay code:
  // 1  K0L -> nu_e pi- e+
  // 2  K0L -> anti_nu_e pi+ e-
  // 3  K0L -> nu_mu pi- mu+
  // 4  K0L -> anti_nu_mu pi+ mu-
  // 5  K+  -> nu_mu mu+
  // 6  K+  -> nu_e pi0 e+
  // 7  K+  -> nu_mu pi0 mu+
  // 8  K-  -> anti_nu_mu mu-
  // 9  K-  -> anti_nu_e pi0 e-
  // 10 K-  -> anti_nu_mu pi0 mu-
  // 11 mu+ -> anti_nu_mu nu_e e+
  // 12 mu- -> nu_mu anti_nu_e e-
  // 13 pi+ -> nu_mu mu+
  // 14 pi- -> anti_nu_mu mu-
  G4Track * theTrack = theStep->GetTrack();
  G4ParticleDefinition * particleType = theTrack->GetDefinition();
  
  //check if the particle is at hadmon or mumon
  if (theStep->GetPostStepPoint()->GetPhysicalVolume()!=NULL){
    if ((theStep->GetPostStepPoint()->GetPhysicalVolume()->GetName()=="PVHadMon")||
	(theStep->GetPostStepPoint()->GetPhysicalVolume()->GetName()=="PVMuMonA1")){
      
      NumiAnalysis* analysis=NumiAnalysis::getInstance();
      analysis->FillHadmmNtuple(*theTrack);
    }
  }
  if (theStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL){
    G4int decay_code=0;
    if (theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="Decay"){
      G4int nSecAtRest = fpSteppingManager->GetfN2ndariesAtRestDoIt();
      G4int nSecAlong  = fpSteppingManager->GetfN2ndariesAlongStepDoIt();
      G4int nSecPost   = fpSteppingManager->GetfN2ndariesPostStepDoIt();
      G4int nSecTotal  = nSecAtRest+nSecAlong+nSecPost;
      G4TrackVector* secVec = fpSteppingManager->GetfSecondary();
      
      if (particleType==G4PionPlus::PionPlusDefinition()) {
	for (size_t partno=(*secVec).size()-nSecTotal;partno<(*secVec).size();partno++)
	  {
	    if ((*secVec)[partno]->GetDefinition()->GetParticleName()=="nu_mu")
	      decay_code=13;
	  }
      }
      if (particleType==G4PionMinus::PionMinusDefinition()) {
	for (size_t partno=(*secVec).size()-nSecTotal;partno<(*secVec).size();partno++)
	  {
	    if ((*secVec)[partno]->GetDefinition()->GetParticleName()=="anti_nu_mu")
	      decay_code=14;
	  }
      }
      if (particleType==G4KaonPlus::KaonPlusDefinition()) {
	for (size_t partno=(*secVec).size()-nSecTotal;partno<(*secVec).size();partno++)
	  {
	    if ((*secVec)[partno]->GetDefinition()->GetParticleName()=="nu_mu")
	      {if (nSecTotal==2) decay_code=5;
	      if (nSecTotal==3) decay_code=7;}
	    if ((*secVec)[partno]->GetDefinition()->GetParticleName()=="nu_e")
	      decay_code=6;
	  }
      }
      if (particleType==G4KaonMinus::KaonMinusDefinition()) {
	for (size_t partno=(*secVec).size()-nSecTotal;partno<(*secVec).size();partno++)
	  {
	    if ((*secVec)[partno]->GetDefinition()->GetParticleName()=="anti_nu_mu")
	      {if (nSecTotal==2) decay_code=8;
	      if (nSecTotal==3) decay_code=10;}
	    if ((*secVec)[partno]->GetDefinition()->GetParticleName()=="anti_nu_e")
	      decay_code=9;
	  }
      }
      if (particleType==G4KaonZeroLong::KaonZeroLongDefinition()) {
	for (size_t partno=(*secVec).size()-nSecTotal;partno<(*secVec).size();partno++)
	  {
	    if ((*secVec)[partno]->GetDefinition()->GetParticleName()=="nu_e")
	      decay_code=1;
	    if ((*secVec)[partno]->GetDefinition()->GetParticleName()=="anti_nu_e")
	      decay_code=2;
	    if ((*secVec)[partno]->GetDefinition()->GetParticleName()=="nu_mu")
	      decay_code=3;
	    if ((*secVec)[partno]->GetDefinition()->GetParticleName()=="anti_nu_mu")
	      decay_code=4;	    
	  }
      }
      if (particleType==G4MuonPlus::MuonPlusDefinition()) {
	for (size_t partno=(*secVec).size()-nSecTotal;partno<(*secVec).size();partno++)
	  {
	    if ((*secVec)[partno]->GetDefinition()->GetParticleName()=="anti_nu_mu")
	      decay_code=11;	
	  }
      }
      if (particleType==G4MuonMinus::MuonMinusDefinition()) {
	for (size_t partno=(*secVec).size()-nSecTotal;partno<(*secVec).size();partno++)
	  {
	    if ((*secVec)[partno]->GetDefinition()->GetParticleName()=="nu_mu")
	      decay_code=12;	
	  }
      }

      NumiTrackInformation* oldinfo=(NumiTrackInformation*)(theTrack->GetUserInformation()); 
      if (oldinfo!=0) {
	oldinfo->SetDecayCode(decay_code);                                                      
	theTrack->SetUserInformation(oldinfo); 
      }
      else {
	NumiTrackInformation* newinfo=new NumiTrackInformation(); 
	newinfo->SetDecayCode(decay_code);                                                       
	theTrack->SetUserInformation(newinfo); 
      }
      /*
      NumiTrackInformation* oldinfo=(NumiTrackInformation*)(theTrack->GetUserInformation());
      NumiTrackInformation* newinfo=new NumiTrackInformation();
      if (oldinfo!=0) newinfo=oldinfo;
      newinfo->SetDecayCode(decay_code);
      theTrack->SetUserInformation(newinfo);
      */
    }
  }
}


