//
// NumiRunManager.cc
//

#include "G4RunManager.hh"
#include "NumiRunManager.hh"
#include "NumiDataInput.hh"
#include "TStopwatch.h"
#include "TTime.h"
#include "NumiPrimaryGeneratorAction.hh"

NumiRunManager::NumiRunManager()
  :primaryGeneratorAction(0)
{;}

NumiRunManager::~NumiRunManager()
{;}
void NumiRunManager::BeamOn(G4int n_event,const char* macroFile,G4int n_select)
{
  nEvents = n_event;
  G4bool runOn(true);
  G4bool cond = ConfirmBeamOnCondition();
   if(cond)
     {
       NumiDataInput *ND = NumiDataInput::GetNumiDataInput();
       //       ND->SetRunNumber(macroFile);
       primaryGeneratorAction = (NumiPrimaryGeneratorAction*)(this)->userPrimaryGeneratorAction;
       if(ND->useFlukaInput || ND->useMarsInput){
	 runOn=primaryGeneratorAction->OpenNtuple(ND->GetExtNtupleFileName());
  	 n_event = primaryGeneratorAction->GetNoOfPrimaries();
       }
       else if(ND->useMuonBeam){   
	 G4cout << "Muon Beam is set" << G4endl;  	 
	 primaryGeneratorAction->SetMuonBeam();
       }
       else{
     	 primaryGeneratorAction->SetProtonBeam();
       }
       if (runOn){
	 TStopwatch *sWatch=new TStopwatch();
	 sWatch->Start(true);
	 numberOfEventToBeProcessed = n_event;
	 RunInitialization();
	 if(n_event>0) DoEventLoop(n_event,macroFile,n_select);
	 RunTermination();
	 
	 if(ND->useFlukaInput){
	   primaryGeneratorAction->CloseNtuple();
	 }
	 
	 sWatch->Stop();
	 G4cout <<"Processed "<<n_event<<" particles in ";
	 G4cout.precision(3);
	 G4cout << sWatch->RealTime()
		<< " s" << G4endl;
	 delete sWatch;
       }
     }
}

