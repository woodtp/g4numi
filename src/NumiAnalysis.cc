//
// NumiAnalysis.cc
//

//Root files
#include "TROOT.h"          // Top level (or root) structure for all classes
#include "TApplication.h"   // ROOT head file for GUI application singleton
#include <TH1.h>            // ROOT head file for 1 dimension histograms
#include <TH2.h>            // ROOT head file for 2 dimension histograms
#include <TSystem.h>        // ROOT head file for a generic interface to the OS
#include <TStopwatch.h>     // ROOT head file for stopwatch for real and cpu time
#include <TStyle.h>         // ROOT head file for all graphics attributes
#include <TFile.h>          
#include <TText.h>
#include <TTree.h>
#include <TF1.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TPaveLabel.h>
#include <TPostScript.h>
#include <TGraph.h>

#include "NumiAnalysis.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>
#include "G4SteppingManager.hh"
#include "G4ThreeVector.hh"
#include "G4TrajectoryContainer.hh"
#include "G4RunManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4Run.hh"
#include "NumiTrackInformation.hh"
#include "NumiDataInput.hh"
#include <stdlib.h>
using namespace std;

NumiAnalysis* NumiAnalysis::instance = 0;

NumiAnalysis::NumiAnalysis()
{
NumiData=NumiDataInput::GetNumiDataInput();
#ifdef G4ANALYSIS_USE
#endif
 
 writeASCII=false;
 if (writeASCII) {
   asciiFileName="numi.out";
   std::ofstream asciiFile(asciiFileName);
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

  G4RunManager* pRunManager=G4RunManager::GetRunManager();
  //sprintf(filename,"%s%d%s","Numi_nt",pRunManager->GetCurrentRun()->GetRunID(),".root");
  char filename[]="Numi_nt.root";
  
  if (NumiData->CreateNuNtuple){
    //Check if file already existed
    if (pRunManager->GetCurrentRun()->GetRunID()>0){
      ntuple = new TFile(filename,"update","root ntuple");
      tree=(TTree*)(ntuple->Get("nudata"));
      tree->SetBranchAddress("data",&g4data);
    }
    else{
      ntuple = new TFile(filename,"recreate","root ntuple");
      tree = new TTree("nudata","g4numi Neutrino ntuple");
      tree->Branch("data",&g4data,"run/I:evtno:beamHWidth/D:beamVWidth:beamX:beamY:protonX:protonY:nuTarZ:hornCurrent:Ndxdz/D:Ndydz:Npz:Nenergy:NdxdzNea:NdydzNea:NenergyN[10]:NWtNear[10]:NdxdzFar:NdydzFar:NenergyF[10]:NWtFar[10]:Norig/I:Ndecay:Ntype:Vx/D:Vy:Vz:pdPx:pdPy:pdPz:ppdxdz:ppdydz:pppz:ppenergy:ppmedium:ptype/I:ppvx/D:ppvy:ppvz:muparpx:muparpy:muparpz:mupare:Necm:Nimpwt:xpoint:ypoint:zpoint:tvx:tvy:tvz:tpx:tpy:tpz:tptype/I:tgen/I:trkx[10]/D:trky[10]:trkz[10]:trkpx[10]:trkpy[10]:trkpz[10]");
    }
  }
  
  if (NumiData->CreateHadmmNtuple){
    char hmmfilename[]="hadmm.root";
    //if (hadmmntuple->GetSize()>1000){
    if(pRunManager->GetCurrentRun()->GetRunID()>0){
      hadmmntuple = new TFile(hmmfilename, "update","hadmm ntuple");
      hadmmtree=(TTree*)(hadmmntuple->Get("hadmm"));
      hadmmtree->SetBranchAddress("data",&g4hmmdata);
    }
    else{
      hadmmntuple = new TFile(hmmfilename, "recreate","hadmm ntuple");
      hadmmtree = new TTree("hadmm","g4numi Hadmmmon  ntuple");
      hadmmtree->Branch("data",&g4hmmdata,"run/I:evtno:mtgthpos/D:mtgtvpos:mtgthsig:mtgtvsig:ptype/I:hmenergy/D:hmxpos:hmypos:hmzpos:hmmpx:hmmpy:hmmpz");
    }
  }
  //book histograms
}

void NumiAnalysis::finish()
{
  if (NumiData->CreateNuNtuple){
    ntuple->cd();
    tree->Write();
    ntuple->Close();
    delete ntuple;
  }

  if (NumiData->CreateHadmmNtuple){
    hadmmntuple->cd();
    hadmmtree->Write();
    hadmmntuple->Close();
    delete hadmmntuple;
  }
}
void NumiAnalysis::FillHadmmNtuple(const G4Track& track)
{
  if (!NumiData->CreateHadmmNtuple) return;
  G4RunManager* pRunManager=G4RunManager::GetRunManager();
  g4hmmdata.run=pRunManager->GetCurrentRun()->GetRunID();
  g4hmmdata.evtno=pRunManager->GetCurrentEvent()->GetEventID();
  G4ParticleDefinition* particleDefinition=track.GetDefinition();
  g4hmmdata.ptype=GetParticleCode(particleDefinition->GetParticleName());
  g4hmmdata.hmmenergy=track.GetTotalEnergy();
  g4hmmdata.hmmxpos=track.GetPosition()[0];
  g4hmmdata.hmmypos=track.GetPosition()[1];
  g4hmmdata.hmmzpos=track.GetPosition()[2];
  g4hmmdata.hmmpx=track.GetMomentum()[0];
  g4hmmdata.hmmpy=track.GetMomentum()[1];
  g4hmmdata.hmmpz=track.GetMomentum()[2]; 
 
  hadmmtree->Fill(); 
 
}
void NumiAnalysis::FillNeutrinoNtuple(const G4Track& track)
{
  if (!NumiData->CreateNuNtuple) return;
  
  
  //Neutrino vertex position and momentum
  G4ThreeVector pos = track.GetPosition()/mm; 
  x = pos.x();
  y = pos.y();
  z = pos.z();
  G4ThreeVector NuMomentum = track.GetMomentum();
  G4int parentID=track.GetParentID();
  NumiTrajectory* NuParentTrack=GetParentTrajectory(parentID);
  G4int point_no=NuParentTrack->GetPointEntries();
  G4ThreeVector ParentMomentumFinal=NuParentTrack->GetMomentum(point_no-1);
  G4ThreeVector vertex_r=(NuParentTrack->GetPoint(point_no-1)->GetPosition()/m)*m; //Should be the same as Neutrino vertex
  G4String parent_name=NuParentTrack->GetParticleName();
  G4double Parent_mass=NuParentTrack->GetMass();
  G4double gamma=sqrt(ParentMomentumFinal*ParentMomentumFinal+Parent_mass*Parent_mass)/Parent_mass; 
  G4double Parent_energy=gamma*Parent_mass;
  G4double beta=sqrt((gamma*gamma-1.)/(gamma*gamma));
  G4ThreeVector beta_vec=ParentMomentumFinal/Parent_energy;
  G4double partial=gamma*(beta_vec*NuMomentum);
 
  G4double enuzr=gamma*(track.GetTotalEnergy())-partial; //neutrino energy in parent rest frame
 
  //fill histograms, ntuples,...

  G4RunManager* pRunManager=G4RunManager::GetRunManager();
  g4data.run=pRunManager->GetCurrentRun()->GetRunID();
  g4data.evtno=pRunManager->GetCurrentEvent()->GetEventID();
  g4data.beamHWidth=NumiData->beamSigmaX/cm;
  g4data.beamVWidth=NumiData->beamSigmaY/cm;
  g4data.beamX=NumiData->beamPosition[0]/cm;
  g4data.beamY=NumiData->beamPosition[1]/cm;
  
  G4int particleID=track.GetParentID();
  NumiTrajectory* dummyTrack=GetParentTrajectory(particleID);
  while (particleID!=1){ 
    particleID=dummyTrack->GetParentID();
    dummyTrack=GetParentTrajectory(particleID);
  }
  g4data.protonX=dummyTrack->GetVertexPosition()[0];
  g4data.protonY=dummyTrack->GetVertexPosition()[1];

  g4data.nuTarZ=NumiData->TargetZ0;
  g4data.hornCurrent=NumiData->HornCurrent/ampere/1000.;

  // Random decay - these neutrinos rarely hit any of the detectors
  g4data.Ndxdz=NuMomentum[0]/NuMomentum[2];
  g4data.Ndydz=NuMomentum[1]/NuMomentum[2];
  g4data.Npz=NuMomentum[2]/GeV;
  g4data.Nenergy=track.GetTotalEnergy()/GeV;

  // Neutrino data at different points
  for (G4int ii=0;ii<10;ii++){          // near detector 
  G4double xdet=NumiData->xdet_near[ii]; 
  G4double ydet=NumiData->ydet_near[ii];
  G4double zdet=NumiData->zdet_near[ii];
 
  g4data.NdxdzNea=(x-NumiData->xdet_near[0])/(z-NumiData->zdet_near[0]);
  g4data.NdydzNea=(y-NumiData->ydet_near[0])/(z-NumiData->zdet_near[0]);

  G4double theta=GetTheta(vertex_r,ParentMomentumFinal,xdet,ydet,zdet); 

  g4data.NenergyN[ii]=GetNuEnergy(enuzr,gamma,beta,theta)/GeV;
  g4data.NWtNear[ii]=GetWeight(track,enuzr,vertex_r,gamma,beta,theta,xdet,ydet,zdet);
  }
  
  
  for (G4int ii=0;ii<10;ii++){         // far detector
  G4double xdet=NumiData->xdet_far[ii]; 
  G4double ydet=NumiData->ydet_far[ii];
  G4double zdet=NumiData->zdet_far[ii];

  g4data.NdxdzFar=(x-NumiData->xdet_far[0])/(z-NumiData->zdet_far[0]);
  g4data.NdydzFar=(y-NumiData->ydet_far[0])/(z-NumiData->zdet_far[0]);

  G4double theta=GetTheta(vertex_r,ParentMomentumFinal,xdet,ydet,zdet);
 
  g4data.NenergyF[ii]=GetNuEnergy(enuzr,gamma,beta,theta)/GeV;
  g4data.NWtFar[ii]=GetWeight(track,enuzr,vertex_r,gamma,beta,theta,xdet,ydet,zdet);
  }
  
  //other info
  // Neutrino origin:
  // 3 From muon decay
  // 1 From particle from target
  // 2 From scraping
  //check if nu is from muon decay or from a particle from target, otherwise Norig=2
  G4int Norig=2;
  if ((parent_name=="mu+") || (parent_name=="mu-")) Norig=3;
  G4String firstvolname=NuParentTrack->GetPreStepVolumeName(0);
  if (firstvolname.contains("HPB") || firstvolname.contains("TGT")) Norig=1;
  
  g4data.Norig=Norig;
  g4data.Ndecay=NuParentTrack->GetDecayCode();

  G4ParticleDefinition * particleType = track.GetDefinition();
  G4int ntype=GetParticleCode(particleType->GetParticleName());
  g4data.Ntype=ntype;
  g4data.Vx=x/cm;
  g4data.Vy=y/cm;
  g4data.Vz=z/cm;
  g4data.pdPx=ParentMomentumFinal[0]/GeV;
  g4data.pdPy=ParentMomentumFinal[1]/GeV;
  g4data.pdPz=ParentMomentumFinal[2]/GeV;

  G4ThreeVector ParentMomentumProduction=NuParentTrack->GetMomentum(0);
  g4data.ppdxdz=ParentMomentumProduction[0]/ParentMomentumProduction[2];
  g4data.ppdydz=ParentMomentumProduction[1]/ParentMomentumProduction[2];
  g4data.pppz=ParentMomentumProduction[2]/GeV; 

  G4double parentp=sqrt(ParentMomentumProduction*ParentMomentumProduction);

  g4data.ppenergy=sqrt((parentp*parentp-Parent_mass*Parent_mass))/GeV;

  g4data.ppmedium=0.; //this is still empty

  g4data.ptype=GetParticleCode(parent_name);

  G4ThreeVector production_vertex=(NuParentTrack->GetPoint(0)->GetPosition()/m)*m; 
  g4data.ppvx=production_vertex[0]/cm;
  g4data.ppvy=production_vertex[1]/cm;
  g4data.ppvz=production_vertex[2]/cm;
  
  //if nu parent is a muon then find muon parent info
  if (parent_name=="mu+" || parent_name=="mu-")
    {
      G4int mupar=NuParentTrack->GetParentID();
      NumiTrajectory* MuParentTrack=GetParentTrajectory(mupar);
      G4int nopoint_mupar=MuParentTrack->GetPointEntries();
      G4ThreeVector muparp=MuParentTrack->GetMomentum(nopoint_mupar-1);
      G4double muparm=MuParentTrack->GetMass();
      g4data.muparpx=muparp[0]/GeV; // vector of hadron parent of muon
      g4data.muparpy=muparp[1]/GeV; // 
      g4data.muparpz=muparp[2]/GeV;
      g4data.mupare=(sqrt(muparp*muparp-muparm*muparm))/GeV;
    }
  else
    {
      g4data.muparpx=-999999.;  
      g4data.muparpy=-999999.;
      g4data.muparpz=-999999.;
      g4data.mupare=-999999.;
    }

  g4data.Necm=enuzr/GeV; // Neutrino energy in parent rest frame

  NumiTrackInformation* info=(NumiTrackInformation*)(track.GetUserInformation());
 
  g4data.Nimpwt=info->GetNImpWt();  // Importance weight


  g4data.xpoint=0.;  // x, y, z of parent at user selected vol
  g4data.xpoint=0.;
  g4data.xpoint=0.;

  /*    
	tgen is is the "generation" number
	of the particle that makes it out of the target. Beam protons have
	tgen=1, any particle produced by a p-C interaction would have tgen=2,
	particles produced from interactions of those products have tgen=3 etc.
	etc. until the cascade exiting the target core.
  */ 
  
  G4int tgen=0;    
  G4int tptype=0;
  G4bool findTarget=false;
  G4ThreeVector ParticleMomentum=G4ThreeVector(-999999,-999999,-999999);
  G4ThreeVector ParticlePosition=G4ThreeVector(-999999,-999999,-999999);
  NumiTrajectory* PParentTrack=GetParentTrajectory(track.GetParentID());
  particleID=PParentTrack->GetTrackID();
  
  while (!findTarget&&particleID!=1){
    G4int numberOfPoints=PParentTrack->GetPointEntries();
    for (G4int ii=0;ii<numberOfPoints-1;ii++){
      G4String lastVolName=PParentTrack->GetPreStepVolumeName(ii);
      G4String nextVolName=PParentTrack->GetPreStepVolumeName(ii+1);      
      if (lastVolName.contains("TGTExit")&&nextVolName.contains("TargetMother"))
	{
	  ParticleMomentum=PParentTrack->GetMomentum(ii);              // tv_ and tp_ are equal to position and  
	  ParticlePosition=PParentTrack->GetPoint(ii)->GetPosition();  // momentum of the particle exiting the target (actually shell around target)
	  NumiTrackInformation* info=(NumiTrackInformation*)(track.GetUserInformation());
	  tgen=info->Gettgen();
	  tptype=GetParticleCode(PParentTrack->GetParticleName());
	  findTarget=true;
	}
    }
    PParentTrack=GetParentTrajectory(PParentTrack->GetParentID());
    particleID=PParentTrack->GetTrackID();
  }

  g4data.tvx=ParticlePosition[0]/cm;
  g4data.tvy=ParticlePosition[1]/cm;
  g4data.tvz=ParticlePosition[2]/cm;
  g4data.tpx=ParticleMomentum[0]/GeV;
  g4data.tpy=ParticleMomentum[1]/GeV;
  g4data.tpz=ParticleMomentum[2]/GeV;
  g4data.tgen=tgen;
  
  g4data.tptype=tptype;
   
  //set all trk_ & trkp_ to -999999  
  for (G4int ii=0;ii<10;ii++){
	g4data.trkx[ii]=-999999.;
	g4data.trky[ii]=-999999.;
	g4data.trkz[ii]=-999999.;
	g4data.trkpx[ii]=-999999.;
	g4data.trkpy[ii]=-999999.;
	g4data.trkpz[ii]=-999999.;
  }
  G4ThreeVector ParentMomentum;
  G4ThreeVector ParentPosition;
  
  for (G4int ii=0;ii<point_no-1;ii++){ 
    ParentMomentum=NuParentTrack->GetMomentum(ii);
    ParentPosition=(NuParentTrack->GetPoint(ii)->GetPosition()/m)*m;
    // check if the particle passes through the neck of the horn
    // if yes then set the trk_ to +999999
    if ((ParentPosition[2]>0.&&ParentPosition[2]<3.*m)&&  // horn 1 position 0-3m
	(ParentPosition[0]*ParentPosition[0]+ParentPosition[1]*ParentPosition[1]<5.*cm)&&
	(g4data.trkx[2]==-999999. || g4data.trkx[2]==999999.))
      {
	g4data.trkx[2]=999999.;
	g4data.trky[2]=999999.;
	g4data.trkz[2]=999999.;  
      }
    if ((ParentPosition[2]>10.*m&&ParentPosition[2]<13.*m)&&  //horn 2 position 10-13m
	(ParentPosition[0]*ParentPosition[0]+ParentPosition[1]*ParentPosition[1]<5.*cm)&&
	(g4data.trkx[4]==-999999. || g4data.trkx[4]==999999.))
      {
	g4data.trkx[4]=999999.;
	g4data.trky[4]=999999.;
	g4data.trkz[4]=999999.;  
      }
      
    G4String postvolname="";
    G4String prevolname=NuParentTrack->GetPreStepVolumeName(ii);
    if (ii<point_no-2) postvolname=NuParentTrack->GetPreStepVolumeName(ii+1);
    
    // parent created inside target
    if ((prevolname.contains("TGT")) && ii==0){
      	g4data.trkx[0]=ParentPosition[0]/cm;
	g4data.trky[0]=ParentPosition[1]/cm;
	g4data.trkz[0]=ParentPosition[2]/cm;
	g4data.trkpx[0]=ParentMomentum[0]/GeV;
	g4data.trkpy[0]=ParentMomentum[1]/GeV;
	g4data.trkpz[0]=ParentMomentum[2]/GeV;}
    //parent at exits target
    if ((prevolname.contains("TGTExit")) && postvolname.contains("TargetMother")){
      	g4data.trkx[1]=ParentPosition[0]/cm;
	g4data.trky[1]=ParentPosition[1]/cm;
	g4data.trkz[1]=ParentPosition[2]/cm;
	g4data.trkpx[1]=ParentMomentum[0]/GeV;
	g4data.trkpy[1]=ParentMomentum[1]/GeV;
	g4data.trkpz[1]=ParentMomentum[2]/GeV;}
    //enter horn1
    if (prevolname.contains("TGAR") && postvolname.contains("Horn1")){
      g4data.trkx[2]=ParentPosition[0]/cm;
      g4data.trky[2]=ParentPosition[1]/cm;
      g4data.trkz[2]=ParentPosition[2]/cm;
      g4data.trkpx[2]=ParentMomentum[0]/GeV;
      g4data.trkpy[2]=ParentMomentum[1]/GeV;
      g4data.trkpz[2]=ParentMomentum[2]/GeV;}
    //exit horn1
    if (prevolname.contains("Horn1") && postvolname.contains("TGAR")){
      g4data.trkx[3]=ParentPosition[0]/cm;
      g4data.trky[3]=ParentPosition[1]/cm;
      g4data.trkz[3]=ParentPosition[2]/cm;
      g4data.trkpx[3]=ParentMomentum[0]/GeV;
      g4data.trkpy[3]=ParentMomentum[1]/GeV;
      g4data.trkpz[3]=ParentMomentum[2]/GeV;}
    //enter horn2
    if (prevolname.contains("TGAR") && postvolname.contains("Horn2")){
      g4data.trkx[4]=ParentPosition[0]/cm;
      g4data.trky[4]=ParentPosition[1]/cm;
      g4data.trkz[4]=ParentPosition[2]/cm;
      g4data.trkpx[4]=ParentMomentum[0]/GeV;
      g4data.trkpy[4]=ParentMomentum[1]/GeV;
      g4data.trkpz[4]=ParentMomentum[2]/GeV;}
    //exit horn2
    if (prevolname.contains("Horn2") && postvolname.contains("TGAR")){
      g4data.trkx[5]=ParentPosition[0]/cm;
      g4data.trky[5]=ParentPosition[1]/cm;
      g4data.trkz[5]=ParentPosition[2]/cm;
      g4data.trkpx[5]=ParentMomentum[0]/GeV;
      g4data.trkpy[5]=ParentMomentum[1]/GeV;
      g4data.trkpz[5]=ParentMomentum[2]/GeV;}
    //enter decay pipe
    if (prevolname.contains("DVOL")&&(postvolname.contains("UpWn"))){
      g4data.trkx[6]=ParentPosition[0]/cm;
      g4data.trky[6]=ParentPosition[1]/cm;
      g4data.trkz[6]=ParentPosition[2]/cm;
      g4data.trkpx[6]=ParentMomentum[0]/GeV;
      g4data.trkpy[6]=ParentMomentum[1]/GeV;
      g4data.trkpz[6]=ParentMomentum[2]/GeV;}
    }

    ParentMomentum=NuParentTrack->GetMomentum(point_no-1);
    ParentPosition=(NuParentTrack->GetPoint(point_no-1)->GetPosition()/m)*m;
    g4data.trkx[7]=ParentPosition[0]/cm;
    g4data.trky[7]=ParentPosition[1]/cm;
    g4data.trkz[7]=ParentPosition[2]/cm;
    g4data.trkpx[7]=ParentMomentum[0]/GeV;
    g4data.trkpy[7]=ParentMomentum[1]/GeV;
    g4data.trkpz[7]=ParentMomentum[2]/GeV;
  
  
  tree->Fill();  // since I already defined where does tree get the data


  // Write to file
  if (writeASCII) {
    std::ofstream asciiFile(asciiFileName, std::ios::app);
    if(asciiFile.is_open()) {
      asciiFile << g4data.Ntype<< " " << g4data.Nenergy << " " << g4data.NenergyN[0] << " " << g4data.NWtNear[0];
      asciiFile << " " << g4data.NenergyF[0] << " " << g4data.NWtFar[0] <<" "<<g4data.Nimpwt<< G4endl; 
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
  G4int ii=0; 
  while (ii<G4int(vect->size())){  
    tr=(*vect)[ii]; 
    NumiTrajectory* tr1=(NumiTrajectory*)(tr);  
    if(tr1->GetTrackID()==parentID) return tr1; 
    ii++; 
  }
  /*
  G4VTrajectory** tr = vect->begin();

  while(tr!=vect->end())
  { 
    NumiTrajectory* tr1 = (NumiTrajectory*)(*tr);
    if(tr1->GetTrackID()==parentID) return tr1;
    tr++;
  }
  */
  return 0;
}

G4double NumiAnalysis::GetWeight(const G4Track& nutrack,G4double enuzr,G4ThreeVector vertex_r,G4double gamma,G4ThreeVector beta_vec,G4double theta_pardet, G4double x_det, G4double y_det, G4double z_det)
{
  G4double beta=beta_vec.mag();
  G4double emrat= 1/(gamma*(1-beta*cos(theta_pardet)));
  
  // get solid angle/4pi for detector element
  G4double rdet=1.;//beam-line.input :FluxAreaR
  G4double sangdet=(rdet*rdet/((z_det/m-vertex_r[2]/m)*(z_det/m-vertex_r[2]/m)))/4.;

  G4double wgt=sangdet*emrat*emrat;

  // done for all except polarized muon decay
  // 
  G4int parentID=nutrack.GetParentID();
  NumiTrajectory* NuParentTrack=GetParentTrajectory(parentID);
  G4String parent_name=NuParentTrack->GetParticleName();

  if ((parent_name=="mu+") || (parent_name=="mu-")){
    G4double rad=sqrt((vertex_r[0]/m-x_det/m)*(vertex_r[0]/m-x_det/m)+(vertex_r[1]/m-y_det/m)*(vertex_r[1]/m-y_det/m)+(vertex_r[2]/m-z_det/m)*(vertex_r[2]/m-z_det/m))*m;
    G4double eneu=enuzr*emrat;
    G4ThreeVector P_nu=G4ThreeVector((x_det-vertex_r[0])/rad*eneu,(y_det-vertex_r[1])/rad*eneu,(z_det-vertex_r[2])/rad*eneu);

    G4double partial=gamma*(beta_vec*P_nu);
    partial=eneu-partial/(gamma+1.);

    G4ThreeVector P_dcm_nu=P_nu-beta_vec*gamma*partial;

    //boost parent of mu to mu production cm
    G4double mu_mass=NuParentTrack->GetMass();
    G4ThreeVector ParentMomentumInitial=NuParentTrack->GetMomentum(0);  
    G4double gamma_mu_init=sqrt(ParentMomentumInitial*ParentMomentumInitial-mu_mass*mu_mass)/mu_mass;
    G4double mu_init_energy=gamma_mu_init*mu_mass;
    G4ThreeVector beta_mu_init_vec=ParentMomentumInitial/mu_init_energy;
    
    G4int muparentID=NuParentTrack->GetParentID();
    NumiTrajectory* MuParentTrack=GetParentTrajectory(muparentID);
    G4int point_no=MuParentTrack->GetPointEntries();
    G4ThreeVector MuParentMomentumFinal=MuParentTrack->GetMomentum(point_no-1);
    partial=gamma_mu_init*(beta_mu_init_vec*MuParentMomentumFinal);
    G4double MuParentMass=MuParentTrack->GetMass();
    G4double MuParentEnergy=sqrt(MuParentMass*MuParentMass+MuParentMomentumFinal*MuParentMomentumFinal);
    
    partial=MuParentEnergy-partial/(gamma_mu_init+1);
    G4ThreeVector P_pcm_mp=MuParentMomentumFinal-beta_mu_init_vec*gamma_mu_init*partial;
    G4double costh=P_dcm_nu*P_pcm_mp/(P_dcm_nu.mag()*P_pcm_mp.mag());
    G4double wt_ratio=0.;
    if ((nutrack.GetDefinition()==G4NeutrinoE::NeutrinoEDefinition()) ||
	(nutrack.GetDefinition()==G4AntiNeutrinoE::AntiNeutrinoEDefinition())) 
      wt_ratio=1.-costh;
    else if ((nutrack.GetDefinition()==G4NeutrinoMu::NeutrinoMuDefinition()) ||
	     (nutrack.GetDefinition()==G4AntiNeutrinoMu::AntiNeutrinoMuDefinition())) {
      G4double xnu=2.*enuzr/mu_mass;
      wt_ratio=((3.-2.*xnu)-(1.-2*xnu)*costh)/(3.-2.*xnu);}
    
    wgt=wgt*wt_ratio;
  }
  return wgt;
}

G4double NumiAnalysis::GetTheta(G4ThreeVector vertex_r,G4ThreeVector momentum,G4double x_det,G4double y_det,G4double z_det)
{
  G4double rad=sqrt((vertex_r[0]/m-x_det/m)*(vertex_r[0]/m-x_det/m)+(vertex_r[1]/m-y_det/m)*(vertex_r[1]/m-y_det/m)+(vertex_r[2]/m-z_det/m)*(vertex_r[2]/m-z_det/m))*m;
  G4double parentp=sqrt(momentum*momentum);
  G4double theta_pardet=((momentum[0])*(-vertex_r[0]/m+x_det/m)+(momentum[1])*(-vertex_r[1]/m+y_det/m)+(momentum[2])*(-vertex_r[2]/m+z_det/m))/(parentp*(rad/m));
  if (theta_pardet>1.) theta_pardet=1.;
  if (theta_pardet<-1.) theta_pardet=-1.;
			 
  return acos(theta_pardet);
}
G4double NumiAnalysis::GetNuEnergy(G4double enuzr,G4double gamma, G4double beta, G4double theta_pardet)
{
  //returns energy of the neutrino passing through some point X; depends on the nu energy in parent rest frame (enuzr) and the angle between parent beta and r(vertex->X)
  
  G4double emrat= 1/(gamma*(1-beta*cos(theta_pardet)));
  return enuzr*emrat;  
}

G4int NumiAnalysis::GetParticleCode(G4String particle_name)
{
  //returns particle code
  G4int ptype=0;
  if (particle_name=="nu_tau") ptype=0; 
  if (particle_name=="anti_nu_tau") ptype=0;

  if (particle_name=="eta_prime") ptype=0; //?

  if (particle_name=="mu+") ptype=5;
  if (particle_name=="mu-") ptype=6;
  if (particle_name=="pi0") ptype=7;
  if (particle_name=="pi+") ptype=8;
  if (particle_name=="pi-") ptype=9;
  if (particle_name=="kaon0L") ptype=10;
  if (particle_name=="kaon+") ptype=11;
  if (particle_name=="kaon-") ptype=12;
  if (particle_name=="neutron") ptype=13;
  if (particle_name=="proton") ptype=14;
  if (particle_name=="anti_proton") ptype=15; 
  if (particle_name=="kaon0S") ptype=16;
  if (particle_name=="eta") ptype=17; 
  if (particle_name=="lambda") ptype=18;
  if (particle_name=="sigma+") ptype=19; 
  if (particle_name=="sigma0") ptype=20;
  if (particle_name=="sigma-") ptype=21;
  if (particle_name=="xi0") ptype=22; 
  if (particle_name=="xi-") ptype=23;
  if (particle_name=="omega-") ptype=24; 
  if (particle_name=="anti_neutron") ptype=25; 
  if (particle_name=="anti_lambda") ptype=26;
  if (particle_name=="anti_sigma-") ptype=27; 
  if (particle_name=="anti_sigma0") ptype=28;
  if (particle_name=="anti_sigma+") ptype=29;
  if (particle_name=="anti_xi0") ptype=30; 
  if (particle_name=="anti_xi-") ptype=31;//?
  if (particle_name=="anti_nu_e") ptype=52;
  if (particle_name=="nu_e") ptype=53;
  if (particle_name=="anti_nu_mu") ptype=55;
  if (particle_name=="nu_mu") ptype=56;

  if (ptype==0) G4cout<<"NumiAnalysis: "<<particle_name<<" code not found"<<G4endl;

  return ptype;
}
