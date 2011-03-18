//----------------------------------------------------------------------
// NumiSteppingAction.cc
// $Id: NumiSteppingAction.cc,v 1.16.4.3 2011/03/18 18:31:12 loiacono Exp $
//----------------------------------------------------------------------

//C++
#include <string>

#include "G4ProcessManager.hh"
#include "NumiSteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "NumiTrajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4RunManager.hh"
#include "NumiTrackInformation.hh"
#include "NumiDataInput.hh"
#include "NumiAnalysis.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Event.hh"
#include "NumiTrackInformation.hh"
#include "NumiAnalysis.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4EventManager.hh"
#include "NumiEventAction.hh"
#include "NumiRunManager.hh"

NumiSteppingAction::NumiSteppingAction()
   :fPrintAllSteps(false),
    fPrintSplitting(false),
    fPrintMuAlcove1(false),
    fPrintMuAlcove2(false),
    fPrintMuAlcove3(false),
    fPrintDeltaAlcove1(false),
    fPrintDeltaAlcove2(false),
    fPrintDeltaAlcove3(false),
    fPrintProcesses(false),
    fPrintTouchableHistory(false)
{
  NDI = NumiDataInput::GetNumiDataInput();
  pRunManager=(NumiRunManager*)NumiRunManager::GetRunManager();

  
}

NumiSteppingAction::~NumiSteppingAction()
{
}


void NumiSteppingAction::UserSteppingAction(const G4Step * theStep)
{

   if(NDI->GetDebugLevel() > 3)
   {
      G4int evtno = pRunManager->GetCurrentEvent()->GetEventID();
      std::cout << "Event " << evtno << ": NumiSteppingAction::UserSteppingAction() Called." << std::endl;
   }

   
   NumiAnalysis* analysis = NumiAnalysis::getInstance();
   analysis->FillBXDRAW(theStep);
   
   G4Track * theTrack = theStep->GetTrack();
   G4ParticleDefinition * particleDefinition = theTrack->GetDefinition();
   G4String particleType = particleDefinition -> GetParticleType();
   
   
   
   
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
   
   
   
   //*****************************************************************
   //*****************************************************************
   //This is for hadron production simulation
   //
   //
   if(NDI->createTarNtuple){
      if( (theStep->GetPreStepPoint()->GetPhysicalVolume() != NULL) && (theStep->GetPostStepPoint()->GetPhysicalVolume() != NULL)){
         G4String testNamePre = theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
         G4String testNamePost = theStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
         if( (testNamePre.contains("TGTExit")) && (testNamePost.contains("TargetMother"))){
            NumiAnalysis* analysis=NumiAnalysis::getInstance();
            analysis->FillTarNtuple(*theTrack);
            theTrack->SetTrackStatus(fStopAndKill);
         }      
      }
   }
   //
   //
   //*****************************************************************
   //*****************************************************************




           
  
  if (!(NDI->useMuonBeam) && (NDI->GetKillTracking() && theTrack->GetKineticEnergy() < NDI->GetKillTrackingThreshold() ) &&
      (particleDefinition != G4NeutrinoE::NeutrinoEDefinition())&&
      (particleDefinition != G4NeutrinoMu::NeutrinoMuDefinition()) &&
      (particleDefinition != G4NeutrinoTau::NeutrinoTauDefinition()) &&
      (particleDefinition != G4AntiNeutrinoE::AntiNeutrinoEDefinition()) &&
      (particleDefinition != G4AntiNeutrinoMu::AntiNeutrinoMuDefinition()) &&
      (particleDefinition != G4AntiNeutrinoTau::AntiNeutrinoTauDefinition()))
    {
      theTrack->SetTrackStatus(fStopAndKill);
    }

  
  //*****************************************************************
  //*****************************************************************
  //This is for muon monitor stuff 
  //
  //

   if(fPrintAllSteps)
   {
      if(analysis->GetEntry() == 107 &&
         (theTrack -> GetTrackID() == 1828)
         )
      {
         G4cout << "   Parent muon entry # = " << analysis->GetEntry()
                << "   track ID = " << theTrack->GetTrackID()
                << "   track weight = " << theTrack->GetWeight()
                << "   post vol = " << theStep->GetPostStepPoint()->GetPhysicalVolume()->GetName()
                << "   pre vol = " << theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() << G4endl;
         
      }
   }
  
   if(NDI->useMuonBeam)
   {
      
      G4double newWeight = 1.0;
      
      G4bool reWeightedDelta = false;
      
      G4EventManager* eventManager    = G4EventManager::GetEventManager();
      G4StackManager* stackManager    = eventManager -> GetStackManager();
      G4TrackingManager* trackManager = eventManager -> GetTrackingManager();
      G4SteppingManager* stepManager  = trackManager -> GetSteppingManager();
      
      NumiEvtAct = (NumiEventAction*)(eventManager -> GetUserEventAction());
      
      if(NDI->GetSimDRays() && NDI->GetReWeightDeltas())
      {
         
         
         G4int nSplitDeltas = NDI->GetNSplitDeltas();
         
         if(nSplitDeltas < 2)
            G4cout << "NumiSteppingAction::UserSteppingAction - PROBLEM: requesting to split deltas but number of spilts is < 2" << G4endl;
         
         if(theStep->GetPreStepPoint()->GetPhysicalVolume() != NULL &&
            theStep->GetPostStepPoint()->GetPhysicalVolume() != NULL)
         {
            
            std::string preStepPointName  = theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
            std::string postStepPointName = theStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
            if((particleDefinition == G4Electron::ElectronDefinition() || 
                particleDefinition == G4Positron::PositronDefinition() ||
                particleDefinition == G4Gamma::Gamma()) &&
               preStepPointName.find("MuCell") == std::string::npos &&
               (postStepPointName.find("MuMonAlcvShot_") != std::string::npos ||
                postStepPointName.find("MuMonAlcv_") != std::string::npos ||
                postStepPointName.find("MuMon_") != std::string::npos) &&
               theStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL &&
               theStep->GetPostStepPoint()->GetProcessDefinedStep() -> GetProcessName() == "Transportation" &&
               !(theTrack->GetWeight() < 1.0))
            {
               if(fPrintSplitting)
               {
                  G4cout << "Spilling " << particleDefinition -> GetParticleName() << " : " 
                         << " parent muon entry # = " << analysis->GetEntry()
                         << " track id = " << theTrack ->GetTrackID()
                         << " parent id = " << theTrack ->GetParentID() 
                         << " post step point = " << postStepPointName << G4endl;
               }
               
               const G4Track * theTempTrack = theStep->GetTrack();
               
               newWeight = newWeight/(double)nSplitDeltas;
               
               reWeightedDelta = true;
               
               for(G4int isplit = 1; isplit < nSplitDeltas; ++isplit)
               {
                  //
                  //create a copy of the dynamic particle
                 //
                 //G4DynamicParticle* newDP = new G4DynamicParticle(theTrack->GetDynamicParticle());
                  
                 //
                 //create new track
                 //
                 //G4Track* newTrack = new G4Track(newDP, theTrack->GetGlobalTime(), theTrack->GetPosition());   
                  G4Track* newTrack = new G4Track(*theTempTrack);
                 //
                 //   Touchable handle is copied to keep the pointer
                 //
                 newTrack->SetTouchableHandle(theTempTrack->GetTouchableHandle());
                 newTrack->SetCreatorProcess(theTempTrack->GetCreatorProcess());
                 
                 //
                 //set track ID to be the SAME!!!
                 //set parent ID to the parent of theTempTrack
              //
                 newTrack->SetParentID( theTempTrack->GetParentID() );
                 //newTrack->SetTrackID( theTempTrack->GetTrackID() );
                 NumiEvtAct -> IncrementDRayTrackID();
                 newTrack->SetTrackID( NumiEvtAct -> GetDRayTrackID() );
                 
                 newTrack->SetWeight(newWeight);
                 
                 
                 if(fPrintSplitting)
                 {
                    G4cout << " created spilt delta ray with trackID = " << newTrack->GetTrackID()
                        << " and weight = " << newTrack->GetWeight() << G4endl;
                 }  
                 //(trackManager -> GetTrack()) -> SetWeight(0.666666665);
                 
                 
                 
                 //
                 //I could probably set the weight of theTempTrack
                 //here because if the postStepPoint = MuMonAlcv then
                 //there is no chance that theTempTrack could make it into the alcove on this loop.
                 //But to be safe flag it to set the weight at the end of the loop.
                 //
              
                 
                 G4int nstep = theTempTrack->GetCurrentStepNumber();
                 for(G4int istep = 0; istep < nstep; ++istep)
                    newTrack->IncrementCurrentStepNumber();
                 
                 
                 //G4Step* newStep = new G4Step(theStep);
                 //newTrack->SetStep(newStep);
                 //G4TrackVector *trackVec = new G4TrackVector();
                 //trackVec -> push_back(newTrack);
                 //eventManager -> StackTracks(trackVec);
                 
                 stackManager -> PushOneTrack(newTrack);
              }
           }
         }
      }//end if sim drays and reweight deltas
      
           
     
     if ( (NDI->GetKillTracking() && theTrack->GetKineticEnergy() < NDI->GetKillTrackingThreshold()) )
     {
	theTrack->SetTrackStatus(fStopAndKill);
     }
     
     if( particleDefinition != G4MuonPlus::MuonPlusDefinition()
         && particleDefinition != G4MuonMinus::MuonMinusDefinition() 
         && particleDefinition != G4Gamma::GammaDefinition() 
         && particleDefinition != G4Electron::ElectronDefinition() 
         && particleDefinition != G4Positron::PositronDefinition())      
     {
	//NumiAnalysis* analysis = NumiAnalysis::getInstance();
	//G4cout << "Event " << analysis->GetEntry() << ", Process = "
        //<< theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl;
	//G4cout<< "Event " << analysis->GetEntry() << ", KE = " << theTrack->GetKineticEnergy()
        //<< ", threshold = " << NDI->GetKillTrackingThreshold()  << ", type = "
        //<< particleDefinition << ", zpos = " << theTrack->GetPosition()[2] << G4endl; 
	
	theTrack->SetTrackStatus(fStopAndKill);
        if( !(NDI->GetSimDRays()) && particleDefinition != G4MuonPlus::MuonPlusDefinition()
            && particleDefinition != G4MuonMinus::MuonMinusDefinition() )      
        {
           theTrack->SetTrackStatus(fStopAndKill);
        }
     }
     else if( !(NDI->GetSimDRays()) && theStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL 
              && theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="Decay")
     {
        //NumiAnalysis* analysis = NumiAnalysis::getInstance();
	
        //G4cout<< "Event " << analysis->GetEntry() << ": Particle decayed" << ", zpos = "
        //<< theTrack->GetPosition()[2] << G4endl; 
	theTrack->SetTrackStatus(fStopAndKill);
     }
     
     /* if(theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()!="Transportation" &&
        theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()!="muIoni" &&
        theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()!="muBrems" &&
        theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()!="msc" &&
        theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()!="muPairProd")*/
     {
        
     }

     if(fPrintProcesses && particleDefinition == G4Gamma::Gamma())
     {
        G4ProcessManager* pm = particleDefinition -> GetProcessManager();
        G4cout << "particle = " << particleDefinition -> GetParticleName() << G4endl;
        G4ProcessVector* pv = pm -> GetProcessList();
        G4cout << " process list: " ;
        for(int i = 0; i< pv->size(); ++i)
        {
           G4VProcess* proc = (*pv)[i];
           //proc -> SetVerboseLevel(2);
           G4cout <<  proc->GetProcessName() << ", ";
        }
        G4cout << G4endl;
        for(int i = 0; i< pv->size(); ++i)
        {
           G4VProcess* proc = (*pv)[i];
           G4cout <<  proc->GetProcessTypeName(proc->GetProcessType()) << ", ";
        }
        G4cout << G4endl;
        G4ProcessVector* pv_along = pm -> GetAlongStepProcessVector();
        G4cout << " along step process list: " ;
        for(int i = 0; i< pv_along->size(); ++i)
        {
           G4VProcess* proc = (*pv_along)[i];
           G4cout <<  proc->GetProcessName() << ", ";
        }
        G4cout << G4endl;


        if(theStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL )
           G4cout << "Process = " << theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl;
     }
     


     
     if( !(NDI->GetSimDRays()) &&  (!(NDI->GetMaterialSigma() < 0) && !(NDI->GetMaterialSigma() > 0)) )
     {
        if(particleDefinition == G4MuonPlus::MuonPlusDefinition()
           || particleDefinition == G4MuonMinus::MuonMinusDefinition())
        {
           //
           //To see if the particle has just past the point that is 0.5m upstream into the rock before
           //monitor. Record the particle momentum in order to apply delta ray corrections. --LL
           //
           //for alc0; info will only be recorded the first time this executes
           //
           if( theTrack->GetPosition()[2] > 729811.) 
           {
              NumiAnalysis* analysis = NumiAnalysis::getInstance();
              analysis->FillAlcEdepInfo(*theTrack, 0);
           }
           //
           //for alc1; info will only be recorded the first time this executes
           //
           if( theTrack->GetPosition()[2] > 750941.) 
           {
              NumiAnalysis* analysis = NumiAnalysis::getInstance();
              analysis->FillAlcEdepInfo(*theTrack, 1);
           }
           //
           //for alc2; info will only be recorded the first time this executes
           //
           if( theTrack->GetPosition()[2] > 772062.) 
           {
              NumiAnalysis* analysis = NumiAnalysis::getInstance();
              analysis->FillAlcEdepInfo(*theTrack, 2);
           }
        }
     }
     
     //
     // Checks to see whether the particle has entered the Hadron
     // or muon monitors, and if so calls the NumiAnalysis class
     // to record the particle properties through the monitors.-DJK
     //
     if (theStep->GetPostStepPoint()->GetPhysicalVolume()==NULL
         && theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() != "ROCK")
     {
        G4cerr << "Problem in NumiSteppingAction::UserSteppingAction - PostStepVolume is Null"
               <<" PreStepVolume is " << theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() << G4endl;
     }
     
     if (NDI->createHadmmNtuple && theStep->GetPostStepPoint()->GetPhysicalVolume()!=NULL)
     {
        
        if(theStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "PVHadMon"
           && theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "ConcShield"
           && (particleDefinition == G4MuonPlus::MuonPlusDefinition()
               || particleDefinition == G4MuonMinus::MuonMinusDefinition()) )
        {
           //NumiAnalysis* analysis = NumiAnalysis::getInstance();
           //analysis->FillHadmmNtuple(*theTrack, 4, 0);
        }
        //else if( theStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "MuCell" ||
        //         (theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "MuCell"
        //         && theStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "ChamberLayer") )
        else if( theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "MuCell")
        {
           
           
           G4TouchableHistory* theTouchable = (G4TouchableHistory*)(theStep->GetPostStepPoint()->GetTouchable());
           std::string CurrVol = theTouchable->GetVolume(0) -> GetName();

           
           
           int cellNum = -99;
           std::string Alcove = "";

           cellNum = theStep->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();
           

           Alcove  = theTouchable -> GetVolume(1) -> GetName();
           
           


           //////////////////////////////////////////////////////
           

           if(fPrintTouchableHistory)
           {
              G4int histDepth = theTouchable->GetHistoryDepth();       // depth of enclosion
              G4cout << G4endl;
              
              G4cout << "parent muon entry # = " << analysis->GetEntry()
                     << " particle: " << particleDefinition->GetParticleName()
                     << " track id = " << theTrack ->GetTrackID()
                     << " track weight = " << theTrack ->GetWeight() 
                     << "   alcove = " << Alcove  << "  cellNum = " << cellNum << G4endl;
              G4cout << "   edep = " <<   theStep->GetTotalEnergyDeposit() << G4endl;
              G4cout << "   post vol = " << theStep->GetPostStepPoint()->GetPhysicalVolume()->GetName()
                     << "   pre vol = " << theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() << G4endl;
              G4cout << "    depth = " << histDepth;
              for(int i= 0; i <= histDepth; ++i)
              {
                 G4cout << " " << theTouchable->GetVolume(i) -> GetName();
              }
              G4cout << G4endl;
              
           }
           
           ////////////////////////////////////////////////////
        


           if(Alcove.find("MuMon") == std::string::npos)
           {
              G4cerr << "Problem in NumiSteppingAction::UserSteppingAction - Haveing a problem"
                     << " getting alcove from Touchable" << G4endl;
           }

           NumiAnalysis* analysis = NumiAnalysis::getInstance();
///muons
           if(particleDefinition == G4MuonPlus::MuonPlusDefinition()
              || particleDefinition == G4MuonMinus::MuonMinusDefinition())
           {
              //
              //The first time through the preStep point has to be MuCellLayer
              //when it is exiting we don't care as far as the mmp and mmpos is
              //concerned because we only fill thes the first time it is in MuCell
              //
              
              if(Alcove.find("_0") != std::string::npos)
              {
                 analysis->FillHadmmNtuple(*theTrack, 
                                           0, 
                                           cellNum);
                 
                 analysis->FillEdep(theStep->GetTotalEnergyDeposit(), particleDefinition, 0, 1, -99999.0,
                                    cellNum, theTrack->GetTrackID(), theTrack->GetWeight());

                 if(fPrintMuAlcove1)
                 {
                    /*
                    G4ProcessManager* pm = particleDefinition -> GetProcessManager();
                    G4cout << "particle = " << particleDefinition -> GetParticleName() << G4endl;
                    G4ProcessVector* pv = pm -> GetProcessList();
                    G4cout << " process list: " ;
                    for(int i = 0; i< pv->size(); ++i)
                    {
                       G4VProcess* proc = (*pv)[i];
                       G4cout <<  proc->GetProcessName() << ", ";
                    }
                    G4cout << G4endl;
                    G4ProcessVector* pv_along = pm -> GetAlongStepProcessVector();
                    G4cout << " along step process list: " ;
                    for(int i = 0; i< pv_along->size(); ++i)
                    {
                       G4VProcess* proc = (*pv_along)[i];
                       G4cout <<  proc->GetProcessName() << ", ";
                    }
                    G4cout << G4endl;
                    */
                    
                    G4cout << "parent muon entry # = " << analysis->GetEntry()
                           << " Alcove = " << Alcove
                           << " track id = " << theTrack ->GetTrackID()
                           << " track weight = " << theTrack ->GetWeight() << G4endl;
                    G4cout << " post process = " ;
                    if(theStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL)
                       G4cout << theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
                    G4cout << " pre process = " ;
                    if(theStep->GetPreStepPoint()->GetProcessDefinedStep() != NULL )   
                       G4cout << theStep->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName();
                    G4cout << " energy dep = " << theStep->GetTotalEnergyDeposit()
                           << " cell = " << cellNum << G4endl;
                 }
              }
              else if(Alcove.find("_1") != std::string::npos)
              {
                 analysis->FillHadmmNtuple(*theTrack,
                                           1,
                                           cellNum);
                 
                 analysis->FillEdep(theStep->GetTotalEnergyDeposit(), particleDefinition, 1, 1, -99999.0,
                                    cellNum, theTrack->GetTrackID(), theTrack->GetWeight());

                 if(fPrintMuAlcove2)
                 {
                    /*
                    G4ProcessManager* pm = particleDefinition -> GetProcessManager();
                    G4cout << "particle = " << particleDefinition -> GetParticleName() << G4endl;
                    G4ProcessVector* pv = pm -> GetProcessList();
                    G4cout << " process list: " ;
                    for(int i = 0; i< pv->size(); ++i)
                    {
                       G4VProcess* proc = (*pv)[i];
                       G4cout <<  proc->GetProcessName() << ", ";
                    }
                    G4cout << G4endl;
                    G4ProcessVector* pv_along = pm -> GetAlongStepProcessVector();
                    G4cout << " along step process list: " ;
                    for(int i = 0; i< pv_along->size(); ++i)
                    {
                       G4VProcess* proc = (*pv_along)[i];
                       G4cout <<  proc->GetProcessName() << ", ";
                    }
                    G4cout << G4endl;
                    */
                    G4cout << "parent muon entry # = " << analysis->GetEntry()
                           << " Alcove = " << Alcove
                           << " track id = " << theTrack ->GetTrackID()
                           << " track weight = " << theTrack ->GetWeight() << G4endl;
                    G4cout << " post process = " ;
                    if(theStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL)
                       G4cout << theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
                    G4cout << " pre process = " ;
                    if(theStep->GetPreStepPoint()->GetProcessDefinedStep() != NULL )   
                       G4cout << theStep->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName();
                    G4cout << " energy dep = " << theStep->GetTotalEnergyDeposit()
                           << " cell = " << cellNum << G4endl;
                 }
              }
              else if(Alcove.find("_2") != std::string::npos)
              {
                 analysis->FillHadmmNtuple(*theTrack,
                                           2,
                                           cellNum);

                 analysis->FillEdep(theStep->GetTotalEnergyDeposit(), particleDefinition, 2, 1, -99999.0,
                                    cellNum, theTrack->GetTrackID(), theTrack->GetWeight());

                 if(fPrintMuAlcove3)
                 {
                    /*
                    G4ProcessManager* pm = particleDefinition -> GetProcessManager();
                    G4cout << "particle = " << particleDefinition -> GetParticleName() << G4endl;
                    G4ProcessVector* pv = pm -> GetProcessList();
                    G4cout << " process list: " ;
                    for(int i = 0; i< pv->size(); ++i)
                    {
                       G4VProcess* proc = (*pv)[i];
                       G4cout <<  proc->GetProcessName() << ", ";
                    }
                    G4cout << G4endl;
                    G4ProcessVector* pv_along = pm -> GetAlongStepProcessVector();
                    G4cout << " along step process list: " ;
                    for(int i = 0; i< pv_along->size(); ++i)
                    {
                       G4VProcess* proc = (*pv_along)[i];
                       G4cout <<  proc->GetProcessName() << ", ";
                    }
                    G4cout << G4endl;
                    */
                    G4cout << "parent muon entry # = " << analysis->GetEntry()
                           << " Alcove = " << Alcove
                           << " track id = " << theTrack ->GetTrackID()
                           << " track weight = " << theTrack ->GetWeight() << G4endl;
                    G4cout << " post process = " ;
                    if(theStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL)
                       G4cout << theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
                    G4cout << " pre process = " ;
                    if(theStep->GetPreStepPoint()->GetProcessDefinedStep() != NULL )   
                       G4cout << theStep->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName();
                    G4cout << " energy dep = " << theStep->GetTotalEnergyDeposit()
                           << " cell = " << cellNum << G4endl;
                 }
              }
              else
              {
                 G4cout << "NumiSteppingAction::UserSteppingAction - PROBLEM: CAN'T GET ALCOVE FROM " << Alcove << G4endl;
              }
           }
///deltas
           if(particleDefinition == G4Gamma::Gamma())
           {
              G4double edep = theStep->GetTotalEnergyDeposit();
              if(edep > 0.0)
                 G4cout << "NumiSteppingAction::UserSteppingAction - Gamma with edep = " << edep
                        << " and weight = " << theTrack->GetWeight() << " in alcove " << Alcove << G4endl;

           }
           
           if(particleDefinition == G4Electron::ElectronDefinition()
              || particleDefinition == G4Positron::PositronDefinition())
           {

              if(reWeightedDelta)
              {
                 G4cout << "NumiSteppingAction::UserSteppingAction - PROBLEM: reweighted this delta on this loop "
                        << " and it is making it into a pixel. This could be a problem. I don't think it should happen." << G4endl;
              }
              
              
              G4int IntExt = NumiEvtAct -> GetIntExt(theTrack);
             
              if(Alcove.find("_0") != std::string::npos)
              {
                 analysis->FillEdep(theStep->GetTotalEnergyDeposit(), particleDefinition, 0, IntExt,
                                    NumiEvtAct->GetZPos(theTrack), cellNum, theTrack->GetTrackID(), theTrack->GetWeight());

                 if(fPrintDeltaAlcove1)
                 {
                    /*
                    G4ProcessManager* pm = particleDefinition -> GetProcessManager();
                    G4cout << "particle = " << particleDefinition -> GetParticleName() << G4endl;
                    G4ProcessVector* pv = pm -> GetProcessList();
                    G4cout << " process list: " ;
                    for(int i = 0; i< pv->size(); ++i)
                    {
                       G4VProcess* proc = (*pv)[i];
                       G4cout <<  proc->GetProcessName() << ", ";
                    }
                    G4cout << G4endl;
                    G4ProcessVector* pv_along = pm -> GetAlongStepProcessVector();
                    G4cout << " along step process list: " ;
                    for(int i = 0; i< pv_along->size(); ++i)
                    {
                       G4VProcess* proc = (*pv_along)[i];
                       G4cout <<  proc->GetProcessName() << ", ";
                    }
                    G4cout << G4endl;
                    */
                    G4cout << "parent muon entry # = " << analysis->GetEntry()
                           << " Alcove = " << Alcove
                           << " track id = " << theTrack ->GetTrackID()
                           << " parent id = " << theTrack ->GetParentID()
                           << " track weight = " << theTrack ->GetWeight() << G4endl;
                    G4cout << " post process = " ;
                    if(theStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL)
                       G4cout << theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
                    G4cout << " pre process = " ;
                    if(theStep->GetPreStepPoint()->GetProcessDefinedStep() != NULL )   
                       G4cout << theStep->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName();
                    G4cout << " energy dep = " << theStep->GetTotalEnergyDeposit()
                           << " cell = " << cellNum << G4endl;
                         
                 }

              }
              else if(Alcove.find("_1") != std::string::npos)
              {
                 analysis->FillEdep(theStep->GetTotalEnergyDeposit(), particleDefinition, 1, IntExt,
                                    NumiEvtAct->GetZPos(theTrack), cellNum, theTrack->GetTrackID(), theTrack->GetWeight());
                          
                 if(fPrintDeltaAlcove2)
                 {
                    /*
                    G4ProcessManager* pm = particleDefinition -> GetProcessManager();
                    G4cout << "particle = " << particleDefinition -> GetParticleName() << G4endl;
                    G4ProcessVector* pv = pm -> GetProcessList();
                    G4cout << " process list: " ;
                    for(int i = 0; i< pv->size(); ++i)
                    {
                       G4VProcess* proc = (*pv)[i];
                       G4cout <<  proc->GetProcessName() << ", ";
                    }
                    G4cout << G4endl;
                    G4ProcessVector* pv_along = pm -> GetAlongStepProcessVector();
                    G4cout << " along step process list: " ;
                    for(int i = 0; i< pv_along->size(); ++i)
                    {
                       G4VProcess* proc = (*pv_along)[i];
                       G4cout <<  proc->GetProcessName() << ", ";
                    }
                    G4cout << G4endl;
                    */
                    G4cout << "parent muon entry # = " << analysis->GetEntry()
                           << " Alcove = " << Alcove
                           << " track id = " << theTrack ->GetTrackID()
                           << " parent id = " << theTrack ->GetParentID()
                           << " track weight = " << theTrack ->GetWeight() << G4endl;
                    G4cout << " post process = " ;
                    if(theStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL)
                       G4cout << theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
                    G4cout << " pre process = " ;
                    if(theStep->GetPreStepPoint()->GetProcessDefinedStep() != NULL )   
                       G4cout << theStep->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName();
                    G4cout << " energy dep = " << theStep->GetTotalEnergyDeposit()
                           << " cell = " << cellNum << G4endl;
                 }
              }
              else if(Alcove.find("_2") != std::string::npos)
              {
                 analysis->FillEdep(theStep->GetTotalEnergyDeposit(), particleDefinition, 2, IntExt,
                                    NumiEvtAct->GetZPos(theTrack), cellNum, theTrack->GetTrackID(), theTrack->GetWeight());

                 if(fPrintDeltaAlcove3)
                 {
                    /*
                    G4ProcessManager* pm = particleDefinition -> GetProcessManager();
                    G4cout << "particle = " << particleDefinition -> GetParticleName() << G4endl;
                    G4ProcessVector* pv = pm -> GetProcessList();
                    G4cout << " process list: " ;
                    for(int i = 0; i< pv->size(); ++i)
                    {
                       G4VProcess* proc = (*pv)[i];
                       G4cout <<  proc->GetProcessName() << ", ";
                    }
                    G4cout << G4endl;
                    G4ProcessVector* pv_along = pm -> GetAlongStepProcessVector();
                    G4cout << " along step process list: " ;
                    for(int i = 0; i< pv_along->size(); ++i)
                    {
                       G4VProcess* proc = (*pv_along)[i];
                       G4cout <<  proc->GetProcessName() << ", ";
                    }
                    G4cout << G4endl;
                    */
                    G4cout << "parent muon entry # = " << analysis->GetEntry()
                           << " Alcove = " << Alcove
                           << " track id = " << theTrack ->GetTrackID()
                           << " parent id = " << theTrack ->GetParentID()
                           << " track weight = " << theTrack ->GetWeight() << G4endl;
                    G4cout << " post process = " ;
                    if(theStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL)
                       G4cout << theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
                    G4cout << " pre process = " ;
                    if(theStep->GetPreStepPoint()->GetProcessDefinedStep() != NULL )   
                       G4cout << theStep->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName();
                    G4cout << " energy dep = " << theStep->GetTotalEnergyDeposit()
                           << " cell = " << cellNum << G4endl;
                 }
              }
              else
              {
                 G4cout << "NumiSteppingAction::UserSteppingAction - PROBLEM: CAN'T GET ALCOVE FROM " << Alcove << G4endl;
              }
           }

        }
        else if((theStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "ROCK"
                 && theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "MuMonAlcvShot_2_Down"))
        {
           theTrack->SetTrackStatus(fStopAndKill);   
        }
     }

     //
     //for spiltting deltas. see above
     //
     if(reWeightedDelta) (stepManager -> GetfPostStepPoint()) -> SetWeight(newWeight);
     
  }//end if(NDI->useMuonBeam)
   //
  //
  //*****************************************************************
  //*****************************************************************

  
  //*****************************************************************
  //*****************************************************************
  //This is for the muon monitor Absorber background simulation
  //
  if(NDI->GetSimAbsBkg())
  {
     if (theStep->GetPostStepPoint()->GetPhysicalVolume() != NULL)
     {
        if(theStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "HadronAbsorber" &&
           theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "ConcShield")
        {
           if(particleType != "lepton")
           {
              NumiAnalysis* analysis = NumiAnalysis::getInstance();
              analysis -> FillAbsorberBkgrdNtuple(*theTrack);
           }
           theTrack->SetTrackStatus(fStopAndKill);   
        }
     }
  }
  //
  //
  //*****************************************************************
  //*****************************************************************

  //================================
  //=======for Raytracing===========
  //--------------------------------
  if(NDI->raytracing || NDI->createZpNtuple){
     for(G4int in = 0; in<(NDI->NZpoint);in++)//---Jasmine Added
     {
        if((theStep->GetPostStepPoint()->GetPosition()[2]>=NDI->Zpoint[in]&&
            theStep->GetPreStepPoint()->GetPosition()[2]<NDI->Zpoint[in]))
        {
           NumiAnalysis* analysis=NumiAnalysis::getInstance();
           analysis->FillZpNtuple(*theTrack,in);
        }
     }
  }
  //================================
  //--------------------------------

  if (theStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL)
  {
     G4int decay_code=0;
     if (theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="Decay")
     {
	G4int nSecAtRest = fpSteppingManager->GetfN2ndariesAtRestDoIt();
	G4int nSecAlong  = fpSteppingManager->GetfN2ndariesAlongStepDoIt();
	G4int nSecPost   = fpSteppingManager->GetfN2ndariesPostStepDoIt();
	G4int nSecTotal  = nSecAtRest+nSecAlong+nSecPost;
	G4TrackVector* secVec = fpSteppingManager->GetfSecondary();
	
	if (particleDefinition==G4PionPlus::PionPlusDefinition())
        {
           for (size_t partno=(*secVec).size()-nSecTotal;partno<(*secVec).size();partno++)
           {
	      if ((*secVec)[partno]->GetDefinition()->GetParticleName()=="nu_mu")
                 decay_code=13;
           }
	}
	if (particleDefinition==G4PionMinus::PionMinusDefinition())
        {
           for (size_t partno=(*secVec).size()-nSecTotal;partno<(*secVec).size();partno++)
           {
	      if ((*secVec)[partno]->GetDefinition()->GetParticleName()=="anti_nu_mu")
                 decay_code=14;
           }
	}
	if (particleDefinition==G4KaonPlus::KaonPlusDefinition())
        {
           for (size_t partno=(*secVec).size()-nSecTotal;partno<(*secVec).size();partno++)
           {
	      if ((*secVec)[partno]->GetDefinition()->GetParticleName()=="nu_mu")
              {if (nSecTotal==2) decay_code=5;
		if (nSecTotal==3) decay_code=7;}
	      if ((*secVec)[partno]->GetDefinition()->GetParticleName()=="nu_e")
                 decay_code=6;
           }
	}
	if (particleDefinition==G4KaonMinus::KaonMinusDefinition())
        {
           for (size_t partno=(*secVec).size()-nSecTotal;partno<(*secVec).size();partno++)
           {
              if ((*secVec)[partno]->GetDefinition()->GetParticleName()=="anti_nu_mu")
	      {if (nSecTotal==2) decay_code=8;
                 if (nSecTotal==3) decay_code=10;}
              if ((*secVec)[partno]->GetDefinition()->GetParticleName()=="anti_nu_e")
                 decay_code=9;
           }
	}
	if (particleDefinition==G4KaonZeroLong::KaonZeroLongDefinition())
        {
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
	if (particleDefinition==G4MuonPlus::MuonPlusDefinition())
        {
           for (size_t partno=(*secVec).size()-nSecTotal;partno<(*secVec).size();partno++)
           {
	      if ((*secVec)[partno]->GetDefinition()->GetParticleName()=="anti_nu_mu")
                 decay_code=11;	
           }
	}
	if (particleDefinition==G4MuonMinus::MuonMinusDefinition())
        {
           for (size_t partno=(*secVec).size()-nSecTotal;partno<(*secVec).size();partno++)
           {
	      if ((*secVec)[partno]->GetDefinition()->GetParticleName()=="nu_mu")
                 decay_code=12;	
           }
        }
	
	NumiTrackInformation* oldinfo=(NumiTrackInformation*)(theTrack->GetUserInformation()); 
	if (oldinfo!=0)
        {
           oldinfo->SetDecayCode(decay_code);                                                      
	  theTrack->SetUserInformation(oldinfo); 
	}
	else
        {
           NumiTrackInformation* newinfo=new NumiTrackInformation(); 
	  newinfo->SetDecayCode(decay_code);                                                       
	  theTrack->SetUserInformation(newinfo); 
	}
     }
  }

}


