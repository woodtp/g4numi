//----------------------------------------------------------------------
// NumiAnalysis.cc
//
// $Id: NumiAnalysis.cc,v 1.15 2008/02/14 19:30:20 koskinen Exp $
//----------------------------------------------------------------------

#include <fstream>
#include <iomanip>
#include <stdlib.h>

//Root 
#include <TSystem.h>        // ROOT head file for a generic interface to the OS
#include <TStopwatch.h>     // ROOT head file for stopwatch for real and cpu time
#include <TFile.h>          
#include <TTree.h>

//GEANT4 
#include "globals.hh"
#include "G4ios.hh"
#include "G4Track.hh"
#include "G4SteppingManager.hh"
#include "G4ThreeVector.hh"
#include "G4TrajectoryContainer.hh"
#include "G4RunManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4Run.hh"

//g4numi 
#include "data_t.hh"
#include "hadmmtuple_t.hh"
#include "NumiParticleCode.hh"
#include "NumiAnalysis.hh"
#include "NumiTrackInformation.hh"
#include "NumiPrimaryGeneratorAction.hh"
#include "NumiDataInput.hh"
#include "NumiNuWeight.hh"

using namespace std;

NumiAnalysis* NumiAnalysis::instance = 0;

NumiAnalysis::NumiAnalysis()
{
  NumiData = NumiDataInput::GetNumiDataInput();
#ifdef G4ANALYSIS_USE
#endif

  g4data = new data_t();
  g4hmmdata = new hadmmtuple_t();

  // part of the nasty hack the will be expunged later. -J

  g4hmmdata->run = -81579;
  g4hmmdata->mtgthsig = -81579; 
  g4hmmdata->mtgtvsig = -81579; 
  g4hmmdata->mtgthpos = -81579; 
  g4hmmdata->mtgtvpos = -81579; 
  g4hmmdata->evtno = -81579; 
  g4hmmdata->ptype = -81579;
  g4hmmdata->hmmenergy = -81579;

  g4hmmdata->hmmxpos = -81579;
  g4hmmdata->hmmypos = -81579;
  g4hmmdata->hmmzpos = -81579;
  g4hmmdata->hmmpx = -81579;
  g4hmmdata->hmmpy = -81579;
  g4hmmdata->hmmpz = -81579;
  for(Int_t i=0;i<3;i++){
    g4hmmdata->mmxpos[i] = -81579;
    g4hmmdata->mmpx[i] = -81579;
    g4hmmdata->mmypos[i] = -81579;
    g4hmmdata->mmpy[i] = -81579;
    g4hmmdata->mmzpos[i] = -81579;
    g4hmmdata->mmpz[i] = -81579; 
  }
}

NumiAnalysis::~NumiAnalysis()
{ 
#ifdef G4ANALYSIS_USE
  // delete things
#endif
}

NumiAnalysis* NumiAnalysis::getInstance()
{
  if (instance == 0) instance = new NumiAnalysis;
  return instance;
}

void NumiAnalysis::book()
{

  G4RunManager* pRunManager = G4RunManager::GetRunManager();
  if (NumiData->createNuNtuple){
    sprintf(nuNtupleFileName,"%s_%04d%s.root",(NumiData->nuNtupleName).c_str(),pRunManager->GetCurrentRun()->GetRunID(), (NumiData->geometry).c_str());
    nuNtuple = new TFile(nuNtupleFileName,"RECREATE","root ntuple");
    G4cout << "Creating neutrino ntuple: "<<nuNtupleFileName<<G4endl;
    tree = new TTree("nudata","g4numi Neutrino ntuple");
    tree->Branch("data","data_t",&g4data,32000,1);
  }
 
  if (NumiData->createHadmmNtuple){
    sprintf(hadmmNtupleFileName,"%s_%04d%s.root",(NumiData->hadmmNtupleName).c_str(),pRunManager->GetCurrentRun()->GetRunID(), (NumiData->geometry).c_str());
    G4cout << "Creating hadron and muon monitors ntuple: "<<hadmmNtupleFileName<<G4endl;
    hadmmNtuple = new TFile(hadmmNtupleFileName, "RECREATE","hadmm ntuple");
    hadmmtree = new TTree("hadmm","g4numi Hadron and muon monitor ntuple");
    hadmmtree->Branch("hadmmdata","hadmmtuple_t",&g4hmmdata,32000,1);
  }

  if (NumiData->createASCII) {
    G4RunManager* pRunManager = G4RunManager::GetRunManager();
    sprintf(asciiFileName,"%s%s%04d%s",(NumiData->asciiName).c_str(),"_",pRunManager->GetCurrentRun()->GetRunID(),".txt");
    G4cout << "Creating ASCII output file : "<<asciiFileName<<G4endl;
    std::ofstream asciiFile(asciiFileName);
  }

  //book histograms
}

void NumiAnalysis::finish()
{  
  if (NumiData->createNuNtuple){
    nuNtuple->cd();
    tree->Write();
    nuNtuple->Close();
    delete nuNtuple;
  }

  if (NumiData->createHadmmNtuple){
    hadmmNtuple->cd();
    hadmmtree->Write();
    hadmmNtuple->Close();
    delete hadmmNtuple;
  }
}

void NumiAnalysis::FillHadmmNtuple(const G4Track& track, Int_t hmm_num)
{
  if (!NumiData->createHadmmNtuple) return;

  G4RunManager* pRunManager = G4RunManager::GetRunManager();

  g4hmmdata->run = pRunManager->GetCurrentRun()->GetRunID();
  g4hmmdata->mtgthsig = NumiData->beamSigmaX/cm;
  g4hmmdata->mtgtvsig = NumiData->beamSigmaY/cm;
  g4hmmdata->mtgthpos = NumiData->beamPosition[0]/cm;
  g4hmmdata->mtgtvpos = NumiData->beamPosition[1]/cm;
  g4hmmdata->evtno = pRunManager->GetCurrentEvent()->GetEventID();

  G4ParticleDefinition* particleDefinition = track.GetDefinition();

  g4hmmdata->ptype = NumiParticleCode::AsInt(NumiParticleCode::StringToEnum(particleDefinition->GetParticleName()));
  g4hmmdata->hmmenergy = track.GetTotalEnergy();

  if(hmm_num == 4){
    g4hmmdata->hmmxpos = track.GetPosition()[0];
    g4hmmdata->hmmypos = track.GetPosition()[1];
    g4hmmdata->hmmzpos = track.GetPosition()[2];
    g4hmmdata->hmmpx = track.GetMomentum()[0];
    g4hmmdata->hmmpy = track.GetMomentum()[1];
    g4hmmdata->hmmpz = track.GetMomentum()[2]; 
  }
  else{
    g4hmmdata->mmxpos[hmm_num] = track.GetPosition()[0];
    g4hmmdata->mmpx[hmm_num] = track.GetMomentum()[0];
    g4hmmdata->mmypos[hmm_num] = track.GetPosition()[1];
    g4hmmdata->mmpy[hmm_num] = track.GetMomentum()[1];
    g4hmmdata->mmzpos[hmm_num] = track.GetPosition()[2];
    g4hmmdata->mmpz[hmm_num] = track.GetMomentum()[2]; 
  } 
}


void NumiAnalysis::WriteHadmmNtuple(){

  if (!NumiData->createHadmmNtuple) return;
  
  hadmmtree->Fill(); 
  
  g4hmmdata->hmmxpos = -81579;
  g4hmmdata->hmmypos = -81579;
  g4hmmdata->hmmzpos = -81579;
  g4hmmdata->hmmpx = -81579;
  g4hmmdata->hmmpy = -81579;
  g4hmmdata->hmmpz = -81579;

  for(G4int i=0; i<3; i++){
    g4hmmdata->mmxpos[i] = -81579;
    g4hmmdata->mmpx[i] = -81579;
    g4hmmdata->mmypos[i] = -81579;
    g4hmmdata->mmpy[i] = -81579;
    g4hmmdata->mmzpos[i] = -81579;
    g4hmmdata->mmpz[i ] = -81579; 
  }
}

void NumiAnalysis::FillNeutrinoNtuple(const G4Track& track)
{

 
  if (!NumiData->createNuNtuple) return;
    
  
  //Neutrino vertex position and momentum
  G4ThreeVector pos = track.GetPosition()/mm; 
  x = pos.x();
  y = pos.y();
  z = pos.z();
  G4ThreeVector NuMomentum = track.GetMomentum();
  G4int parentID = track.GetParentID();
  
  if (parentID == 0) return; //I have to make some changes so that neutrinos in fluka/mars ntuples don't crash

  NumiTrajectory* NuParentTrack = GetParentTrajectory(parentID);
  G4int point_no = NuParentTrack->GetPointEntries();
  G4ThreeVector ParentMomentumFinal = NuParentTrack->GetMomentum(point_no-1);
  G4ThreeVector vertex_r = (NuParentTrack->GetPoint(point_no-1)->GetPosition()/m)*m; //Should be the same as Neutrino vertex
  G4String parent_name = NuParentTrack->GetParticleName();
  G4double Parent_mass = NuParentTrack->GetMass();
  G4double gamma = sqrt(ParentMomentumFinal*ParentMomentumFinal+Parent_mass*Parent_mass)/Parent_mass; 
  G4double Parent_energy = gamma*Parent_mass;
  G4ThreeVector beta_vec = ParentMomentumFinal/Parent_energy;
  G4double partial = gamma*(beta_vec*NuMomentum);
 
  G4double enuzr = gamma*(track.GetTotalEnergy())-partial; //neutrino energy in parent rest frame

  //fill histograms, ntuples,...
  G4RunManager* pRunManager = G4RunManager::GetRunManager();
  NumiPrimaryGeneratorAction *NPGA = (NumiPrimaryGeneratorAction*)(pRunManager)->GetUserPrimaryGeneratorAction();
 
  g4data->run = pRunManager->GetCurrentRun()->GetRunID();
  g4data->evtno = pRunManager->GetCurrentEvent()->GetEventID();
  g4data->beamHWidth = NumiData->beamSigmaX/cm;
  g4data->beamVWidth = NumiData->beamSigmaY/cm;
  g4data->beamX = NumiData->beamPosition[0]/cm;
  g4data->beamY = NumiData->beamPosition[1]/cm;
 
  G4int particleID = track.GetParentID();
  G4ThreeVector protonOrigin = NPGA->GetProtonOrigin();
  g4data->protonX = protonOrigin[0];
  g4data->protonY = protonOrigin[1];
  g4data->protonZ = protonOrigin[2];

  G4ThreeVector protonMomentum = NPGA->GetProtonMomentum();
  g4data->protonPx = protonMomentum[0];
  g4data->protonPy = protonMomentum[1];
  g4data->protonPz = protonMomentum[2];

  g4data->nuTarZ = NumiData->TargetZ0;
  g4data->hornCurrent = NumiData->HornCurrent/ampere/1000.;

  // Random decay - these neutrinos rarely hit any of the detectors
  g4data->Ndxdz = NuMomentum[0]/NuMomentum[2];
  g4data->Ndydz = NuMomentum[1]/NuMomentum[2];
  g4data->Npz = NuMomentum[2]/GeV;
  g4data->Nenergy = track.GetTotalEnergy()/GeV;

   //other info
  // Neutrino origin:
  // 3 From muon decay
  // 1 From particle from target
  // 2 From scraping
  //check if nu is from muon decay or from a particle from target, otherwise Norig = 2
  G4int Norig = 2;
  if ((parent_name=="mu+") || (parent_name=="mu-")) Norig = 3;
  G4String firstvolname = NuParentTrack->GetPreStepVolumeName(0);
  if (firstvolname.contains("Baffle") || firstvolname.contains("TGT")) Norig = 1;

  g4data->Norig = Norig;
  g4data->Ndecay = NuParentTrack->GetDecayCode();

  G4ParticleDefinition * particleType = track.GetDefinition();
  G4int ntype = NumiParticleCode::AsInt(NumiParticleCode::StringToEnum(particleType->GetParticleName()));
  g4data->Ntype = ntype;
  g4data->Vx = x/cm;
  g4data->Vy = y/cm;
  g4data->Vz = z/cm;
  g4data->pdPx = ParentMomentumFinal[0]/GeV;
  g4data->pdPy = ParentMomentumFinal[1]/GeV;
  g4data->pdPz = ParentMomentumFinal[2]/GeV;

  G4ThreeVector ParentMomentumProduction = NuParentTrack->GetMomentum(0);
  g4data->ppdxdz = ParentMomentumProduction[0]/ParentMomentumProduction[2];
  g4data->ppdydz = ParentMomentumProduction[1]/ParentMomentumProduction[2];
  g4data->pppz = ParentMomentumProduction[2]/GeV; 

  G4double parentp = sqrt(ParentMomentumProduction*ParentMomentumProduction);

  g4data->ppenergy = sqrt((parentp*parentp-Parent_mass*Parent_mass))/GeV;

  g4data->ppmedium = 0.; //this is still empty

  g4data->ptype = NumiParticleCode::AsInt(NumiParticleCode::StringToEnum(parent_name));
 
  G4ThreeVector production_vertex = (NuParentTrack->GetPoint(0)->GetPosition()/m)*m; 
  g4data->ppvx = production_vertex[0]/cm;
  g4data->ppvy = production_vertex[1]/cm;
  g4data->ppvz = production_vertex[2]/cm;
  
  //if nu parent is a muon then find muon parent info
  if ((parent_name=="mu+" || parent_name=="mu-") && NuParentTrack->GetParentID()!=0)
    {
      G4int mupar = NuParentTrack->GetParentID();
      NumiTrajectory* MuParentTrack = GetParentTrajectory(mupar);
      G4int nopoint_mupar = MuParentTrack->GetPointEntries();
      G4ThreeVector muparp = MuParentTrack->GetMomentum(nopoint_mupar-1);
      G4double muparm = MuParentTrack->GetMass();
      g4data->muparpx = muparp[0]/GeV; // vector of hadron parent of muon
      g4data->muparpy = muparp[1]/GeV; // 
      g4data->muparpz = muparp[2]/GeV;
      g4data->mupare = (sqrt(muparp*muparp+muparm*muparm))/GeV;
    }
  else
    {
      g4data->muparpx = -999999.;  
      g4data->muparpy = -999999.;
      g4data->muparpz = -999999.;
      g4data->mupare = -999999.;
    }

  g4data->Necm = enuzr/GeV; // Neutrino energy in parent rest frame
  NumiTrackInformation* info = (NumiTrackInformation*)(track.GetUserInformation());
  g4data->Nimpwt = info->GetNImpWt();  // Importance weight
  g4data->tgen = info->GetTgen()-1;

  g4data->xpoint = 0.;  // x, y, z of parent at user selected vol
  g4data->xpoint = 0.;
  g4data->xpoint = 0.;

  /*    
	tgen is is the "generation" number
	of the particle that makes it out of the target. Beam protons have
	tgen=1, any particle produced by a p-C interaction would have tgen=2,
	particles produced from interactions of those products have tgen=3 etc.
	etc. until the cascade exiting the target core.
  */ 
  
  if(!NumiData->useFlukaInput && !NumiData->useMarsInput) //if not using external ntuple then need to find the particle that exited the target
    {
      G4bool findTarget = false;
      G4ThreeVector ParticleMomentum = G4ThreeVector(-999999,-999999,-999999);
      G4ThreeVector ParticlePosition = G4ThreeVector(-999999,-999999,-999999);
      NumiTrajectory* PParentTrack = GetParentTrajectory(track.GetParentID());
      G4int tptype = NumiParticleCode::AsInt(NumiParticleCode::StringToEnum(PParentTrack->GetParticleName()));
      particleID = PParentTrack->GetTrackID();
      
      while (!findTarget && particleID!=1){
	G4int numberOfPoints = PParentTrack->GetPointEntries();
	for (G4int ii=0;ii<numberOfPoints-1;ii++){
	  G4String lastVolName = PParentTrack->GetPreStepVolumeName(ii);
	  G4String nextVolName = PParentTrack->GetPreStepVolumeName(ii+1);      
	  if (lastVolName.contains("TGTExit") && nextVolName.contains("TargetMother"))
	    {
	      ParticleMomentum = PParentTrack->GetMomentum(ii);              // tv_ and tp_ are equal to position and  
	      ParticlePosition = PParentTrack->GetPoint(ii)->GetPosition();  // momentum of the particle exiting the target (actually shell around target)
	      tptype = NumiParticleCode::AsInt(NumiParticleCode::StringToEnum(PParentTrack->GetParticleName()));
	      findTarget = true;
	    }
	}
	PParentTrack = GetParentTrajectory(PParentTrack->GetParentID());
	particleID = PParentTrack->GetTrackID();
      }
      g4data->tvx = ParticlePosition[0]/cm;
      g4data->tvy = ParticlePosition[1]/cm;
      g4data->tvz = ParticlePosition[2]/cm;
      g4data->tpx = ParticleMomentum[0]/GeV;
      g4data->tpy = ParticleMomentum[1]/GeV;
      g4data->tpz = ParticleMomentum[2]/GeV;
      g4data->tptype = tptype;
    }
  else          // using external ntuple, so set these to whatever comes from that ntuple
    {
      G4ThreeVector ParticlePosition = NPGA->GetParticlePosition();
      g4data->tvx = ParticlePosition[0]/cm;
      g4data->tvy = ParticlePosition[1]/cm;
      g4data->tvz = ParticlePosition[2]/cm;
      G4ThreeVector ParticleMomentum = NPGA->GetParticleMomentum();
      g4data->tpx = ParticleMomentum[0]/GeV;
      g4data->tpy = ParticleMomentum[1]/GeV;
      g4data->tpz = ParticleMomentum[2]/GeV;
      g4data->tptype = NPGA->GetParticleType();
    }   
  
  //set all trk_ & trkp_ to -999999  
  for (G4int ii=0;ii<10;ii++){
    g4data->trkx[ii] = -999999.;
    g4data->trky[ii] = -999999.;
    g4data->trkz[ii] = -999999.;
    g4data->trkpx[ii] = -999999.;
    g4data->trkpy[ii] = -999999.;
    g4data->trkpz[ii] = -999999.;
  }
  // Neutrino data at different points
  // need neutrino parent info to be filled in g4data by this point

  for (G4int ii=0;ii<NumiData->nNear;ii++){          // near detector 
    g4data->NdxdzNear[ii] = (x-NumiData->xdet_near[ii])/(z-NumiData->zdet_near[ii]);
    g4data->NdydzNear[ii] = (y-NumiData->ydet_near[ii])/(z-NumiData->zdet_near[ii]);
    
    NumiNuWeight nuwgh;
    G4double nu_wght;
    G4double nu_energy;
    std::vector<double> r_det;
    r_det.push_back(NumiData->xdet_near[ii]/cm);
    r_det.push_back(NumiData->ydet_near[ii]/cm);
    r_det.push_back(NumiData->zdet_near[ii]/cm);
    nuwgh.GetWeight(g4data, r_det,nu_wght,nu_energy);
    g4data->NenergyN[ii] = nu_energy; //in GeV
    g4data->NWtNear[ii]  = nu_wght;
  }

  for (G4int ii=0;ii<NumiData->nFar;ii++){         // far detector
    g4data->NdxdzFar[ii] = (x-NumiData->xdet_far[ii])/(z-NumiData->zdet_far[ii]);
    g4data->NdydzFar[ii] = (y-NumiData->ydet_far[ii])/(z-NumiData->zdet_far[ii]);
    
    NumiNuWeight nuwgh;
    G4double nu_wght;
    G4double nu_energy;
    std::vector<double> r_det;
    r_det.push_back(NumiData->xdet_far[ii]/cm);
    r_det.push_back(NumiData->ydet_far[ii]/cm);
    r_det.push_back(NumiData->zdet_far[ii]/cm);
    nuwgh.GetWeight(g4data, r_det,nu_wght,nu_energy);
    g4data->NenergyF[ii] = nu_energy; //in GeV
    g4data->NWtFar[ii]   = nu_wght;
   
  }

  //if(parentID!=0){
  G4ThreeVector ParentMomentum;
  G4ThreeVector ParentPosition;
    
  G4bool wasInHorn1 = false;
  G4bool wasInHorn2 = false;  
  for (G4int ii=0;ii<point_no-1;ii++){ 
    ParentMomentum = NuParentTrack->GetMomentum(ii);
    ParentPosition = (NuParentTrack->GetPoint(ii)->GetPosition()/m)*m;
    
    G4String postvolname = "";
    G4String prevolname = NuParentTrack->GetPreStepVolumeName(ii);
    if (ii<point_no-2) postvolname = NuParentTrack->GetPreStepVolumeName(ii+1);
      
    // parent created inside target
    if ((prevolname.contains("TGT")||prevolname.contains("Budal")) && ii==0){
      g4data->trkx[0] = ParentPosition[0]/cm;
      g4data->trky[0] = ParentPosition[1]/cm;
      g4data->trkz[0] = ParentPosition[2]/cm;
      g4data->trkpx[0] = ParentMomentum[0]/GeV;
      g4data->trkpy[0] = ParentMomentum[1]/GeV;
      g4data->trkpz[0] = ParentMomentum[2]/GeV;}
    //parent at exits target
    if ((prevolname.contains("TGTExit")) && postvolname.contains("TargetMother")){
      g4data->trkx[1] = ParentPosition[0]/cm;
      g4data->trky[1] = ParentPosition[1]/cm;
      g4data->trkz[1] = ParentPosition[2]/cm;
      g4data->trkpx[1] = ParentMomentum[0]/GeV;
      g4data->trkpy[1] = ParentMomentum[1]/GeV;
      g4data->trkpz[1] = ParentMomentum[2]/GeV;}
    //enter horn1
    if (prevolname.contains("TGAR") && postvolname.contains("Horn1")){
      g4data->trkx[2] = ParentPosition[0]/cm;
      g4data->trky[2] = ParentPosition[1]/cm;
      g4data->trkz[2] = ParentPosition[2]/cm;
      g4data->trkpx[2] = ParentMomentum[0]/GeV;
      g4data->trkpy[2] = ParentMomentum[1]/GeV;
      g4data->trkpz[2] = ParentMomentum[2]/GeV;
      wasInHorn1 = true;
    }
    //exit horn1
    if (prevolname.contains("Horn1") && postvolname.contains("TGAR")){
      g4data->trkx[3] = ParentPosition[0]/cm;
      g4data->trky[3] = ParentPosition[1]/cm;
      g4data->trkz[3] = ParentPosition[2]/cm;
      g4data->trkpx[3] = ParentMomentum[0]/GeV;
      g4data->trkpy[3] = ParentMomentum[1]/GeV;
      g4data->trkpz[3] = ParentMomentum[2]/GeV;
    }
    //enter horn2
    if (prevolname.contains("TGAR") && postvolname.contains("Horn2")){
      g4data->trkx[4] = ParentPosition[0]/cm;
      g4data->trky[4] = ParentPosition[1]/cm;
      g4data->trkz[4] = ParentPosition[2]/cm;
      g4data->trkpx[4] = ParentMomentum[0]/GeV;
      g4data->trkpy[4] = ParentMomentum[1]/GeV;
      g4data->trkpz[4] = ParentMomentum[2]/GeV;
      wasInHorn2 = true;
    }
    //exit horn2
    if (prevolname.contains("Horn2") && postvolname.contains("TGAR")){
      g4data->trkx[5] = ParentPosition[0]/cm;
      g4data->trky[5] = ParentPosition[1]/cm;
      g4data->trkz[5] = ParentPosition[2]/cm;
      g4data->trkpx[5] = ParentMomentum[0]/GeV;
      g4data->trkpy[5] = ParentMomentum[1]/GeV;
      g4data->trkpz[5] = ParentMomentum[2]/GeV;
    }
    //enter decay pipe
    if (prevolname.contains("DVOL") && (postvolname.contains("UpWn"))){
      g4data->trkx[6] = ParentPosition[0]/cm;
      g4data->trky[6] = ParentPosition[1]/cm;
      g4data->trkz[6] = ParentPosition[2]/cm;
      g4data->trkpx[6] = ParentMomentum[0]/GeV;
      g4data->trkpy[6] = ParentMomentum[1]/GeV;
      g4data->trkpz[6] = ParentMomentum[2]/GeV;}


      
    // check if the particle passes through the neck of the horn
    // if yes then set the trk_ to +999999
    // need to make this work for arbitrary horn position!!
    if ((ParentPosition[2]>0.&&ParentPosition[2]<3.*m)&&  // horn 1 position 0-3m
	(sqrt(ParentPosition[0]*ParentPosition[0]+ParentPosition[1]*ParentPosition[1])<5.*cm)&&
	(g4data->trkx[2]==-999999. || g4data->trkx[2]==999999.))
      {
	g4data->trkx[2] = 999999.;
	g4data->trky[2] = 999999.;
	g4data->trkz[2] = 999999.;  
      }
    if ((ParentPosition[2]>10.*m&&ParentPosition[2]<13.*m)&&  //horn 2 position 10-13m
	(sqrt(ParentPosition[0]*ParentPosition[0]+ParentPosition[1]*ParentPosition[1])<5.*cm)&&
	(g4data->trkx[4]==-999999. || g4data->trkx[4]==999999.))
      {
	g4data->trkx[4] = 999999.;
	g4data->trky[4] = 999999.;
	g4data->trkz[4] = 999999.;  
      }
  }
    
  ParentMomentum = NuParentTrack->GetMomentum(point_no-1);
  ParentPosition = (NuParentTrack->GetPoint(point_no-1)->GetPosition()/m)*m;
  g4data->trkx[7] = ParentPosition[0]/cm;
  g4data->trky[7] = ParentPosition[1]/cm;
  g4data->trkz[7] = ParentPosition[2]/cm;
  g4data->trkpx[7] = ParentMomentum[0]/GeV;
  g4data->trkpy[7] = ParentMomentum[1]/GeV;
  g4data->trkpz[7] = ParentMomentum[2]/GeV;
  // }
  
  tree->Fill();  


  // Write to file
  if (NumiData->createASCII) {
    std::ofstream asciiFile(asciiFileName, std::ios::app);
    if(asciiFile.is_open()) {
      asciiFile << g4data->Ntype<< " " << g4data->Nenergy << " " << g4data->NenergyN[0] << " " << g4data->NWtNear[0];
      asciiFile << " " << g4data->NenergyF[0] << " " << g4data->NWtFar[0] <<" "<<g4data->Nimpwt<< G4endl; 
      asciiFile.close();
    }
  }

}

NumiTrajectory* NumiAnalysis::GetParentTrajectory(G4int parentID)
{
  G4TrajectoryContainer* container = 
    G4RunManager::GetRunManager()->GetCurrentEvent()->GetTrajectoryContainer();
  if(container==0) return 0;

  TrajectoryVector* vect = container->GetVector();
  G4VTrajectory* tr;
  G4int ii = 0; 
  while (ii<G4int(vect->size())){  
    tr = (*vect)[ii]; 
    NumiTrajectory* tr1 = (NumiTrajectory*)(tr);  
    if(tr1->GetTrackID() == parentID) return tr1; 
    ii++; 
  }

  return 0;
}
