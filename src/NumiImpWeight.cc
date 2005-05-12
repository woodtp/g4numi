#include "NumiTrackInformation.hh"
#include "NumiDataInput.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"
#include "NumiImpWeight.hh"
#include "Randomize.hh"
#include "G4Track.hh"
#include "NumiTrajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "NumiAnalysis.hh"

NumiImpWeight::NumiImpWeight()
{ 
}

NumiImpWeight::~NumiImpWeight()
{
}
G4double NumiImpWeight::CalculateImpWeight(const G4Track *aTrack)
{
  //Set weight for primary proton 
   if (aTrack->GetTrackID()==0) return 1;
 
   G4int IMPWEIGHT_NBINS=40;
   G4double impwt=0.;
   G4double psave; 
   G4ParticleDefinition * particleType = aTrack->GetDefinition();
   G4double totE=aTrack->GetTotalEnergy()/GeV;
  
   //Check if particle is neutrino
   if (particleType==G4NeutrinoE::NeutrinoEDefinition()||  
       particleType==G4AntiNeutrinoE::AntiNeutrinoEDefinition()||  
       particleType==G4NeutrinoMu::NeutrinoMuDefinition()||  
       particleType==G4AntiNeutrinoMu::AntiNeutrinoMuDefinition()|| 
       particleType==G4NeutrinoTau::NeutrinoTauDefinition()||  
       particleType==G4AntiNeutrinoTau::AntiNeutrinoTauDefinition())
     {
       return 1.;
     }
   /*
   //Find particle parent and check its weight
   // If parent is weighted, return weight 1 
   G4TrajectoryContainer* container = 
     G4RunManager::GetRunManager()->GetCurrentEvent()->GetTrajectoryContainer();
   if(container!=0){
   
     TrajectoryVector* vect = container->GetVector();
     G4VTrajectory** tr = vect->begin();
  
     while(tr!=vect->end())
       { 
	 NumiTrajectory* tr1 = (NumiTrajectory*)(*tr);
	 if(tr1->GetTrackID()==aTrack->GetParentID()) {
	   if(tr1->GetNImpWt()!=1.) return 1;
	   else break;
	 }
	 tr++;
       }
   }
   */
   //From gnumi/Impwght.F and beam/ffkeyImpWeight.F
   G4int indx = (G4int)(totE+0.5);      // Bins are 1 GeV
   if (indx<1) indx = 1;
   else if (indx>IMPWEIGHT_NBINS) indx = IMPWEIGHT_NBINS;
   
   G4double wEM  = 0.01;
   G4double wHad = (G4double)(indx)/30.0;
   if (wHad>1.) wHad = 1.000;
   if (wHad<0.001) wHad = 0.001;
   //     
   //     Get the survival probability from look up tables
   //     
   if (particleType==G4Gamma::GammaDefinition())                      psave = wEM;
   else if (particleType==G4Electron::ElectronDefinition())           psave = wEM;
   else if (particleType==G4Positron::PositronDefinition())           psave = wEM;
   else if (particleType==G4MuonPlus::MuonPlusDefinition())           psave = 1.;
   else if (particleType==G4MuonMinus::MuonMinusDefinition())         psave = 1.;
   else if (particleType==G4PionZero::PionZeroDefinition())           psave = wHad;
   else if (particleType==G4PionPlus::PionPlusDefinition())           psave = wHad;
   else if (particleType==G4PionMinus::PionMinusDefinition())         psave = wHad;
   else if (particleType==G4KaonZero::KaonZeroDefinition())           psave = wHad;
   else if (particleType==G4KaonPlus::KaonPlusDefinition())           psave = wHad;
   else if (particleType==G4KaonMinus::KaonMinusDefinition())         psave = wHad;
   else if (particleType==G4Neutron::NeutronDefinition())             psave = wHad;
   else if (particleType==G4Proton::ProtonDefinition())               psave = wHad;
   else if (particleType==G4AntiProton::AntiProtonDefinition())       psave = wHad;
   else if (particleType==G4KaonZeroShort::KaonZeroShortDefinition()) psave = wHad;
   else    	                                                      psave = wHad;

   
   
   if (psave > 1.) psave = 1.0;
   if (psave > 0.0) impwt = 1.0/psave;
    if (G4UniformRand()<psave) { //G4UniformRand() returns random number between 0 and 1
     return impwt;
   }

   else return 0.;
   
}
      

