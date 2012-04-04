
#include <vector>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <stdlib.h>

//Root 
#include <TSystem.h>        
#include <TStopwatch.h>    
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
#include "G4Proton.hh"

//g4na49
#include "NA49Analysis.hh"
//#include "ProdTuple_t.hh"

using namespace std;

NA49Analysis* NA49Analysis::instance = 0;

//------------------------------------------------------------------------------------
NA49Analysis::NA49Analysis()
{
#ifdef G4ANALYSIS_USE
#endif
 }
//------------------------------------------------------------------------------------
NA49Analysis::~NA49Analysis()
{ 


#ifdef G4ANALYSIS_USE
#endif
}
//------------------------------------------------------------------------------------
NA49Analysis* NA49Analysis::getInstance()
{
  if (instance == 0) instance = new NA49Analysis;
  return instance;
}
//------------------------------------------------------------------------------------
void NA49Analysis::book()
{
  G4RunManager* pRunManager = G4RunManager::GetRunManager();
    sprintf(NtupleFileName,"QGSP_NA49_%04d.root",pRunManager->GetCurrentRun()->GetRunID());
    FileNtuple = new TFile(NtupleFileName, "RECREATE","pion from p+C ntuple");
   
    //   FileNtuple->cd();
    
    ProdTree = new TTree("pCinfo","g4NA49 info from p+C");
    ProdTree->Branch("npart",&g4Proddata.NPart,"NPart/I");
    ProdTree->Branch("pdg", &g4Proddata.PDG,"PDG[NPart]/I");
    ProdTree->Branch("intType", &g4Proddata.InterType,"IntType[NPart]/I");
    ProdTree->Branch("part_types", &g4Proddata.PartTypes,"PartTypes[12]/I");
    ProdTree->Branch("x",  &g4Proddata.X,"X[NPart][3]/D");
    ProdTree->Branch("p",  &g4Proddata.P,"P[NPart][4]/D");
    ProdTree->Branch("xf", &g4Proddata.XF,"XF[NPart]/D");
    ProdTree->Branch("pt", &g4Proddata.PT,"PT[NPart]/D");
   
}


//------------------------------------------------------------------------------------
void NA49Analysis::finish()
{
    FileNtuple->cd();
     
    ProdTree->Write();
    FileNtuple->Close();
    delete FileNtuple;
}

//------------------------------------------------------------------------------------

void NA49Analysis::FillNtuple(std::vector<TrackInfo_t> trackInfoVec)
{
 G4RunManager* pRunManager = G4RunManager::GetRunManager();
 g4Proddata.NPart= trackInfoVec.size();

 Int_t pdg_t;
 Int_t partNum   = 0;
 Int_t NPiPlus   = 0;
 Int_t NPiMinus  = 0;
 Int_t NPi0      = 0;
 Int_t NKPlus    = 0;
 Int_t NKMinus   = 0;
 Int_t NK0S      = 0;
 Int_t NK0L      = 0;
 Int_t NProtons  = 0; 
 Int_t NAProtons = 0; 
 Int_t NNeutrons = 0; 
 Int_t NANeutrons= 0;
 Int_t NOthers = 0;

 //Variables for PT, XF;
 Double_t XF,PT,Ecm,PL,beta,gamma,massPart,Pxx,Pyy,Pzz,PartE;
 Double_t BeamEnergy = 158.003*1000.;//Mev
 //Double_t massProton = 0.938242046*1000.;//MeV
 Double_t massProton = 0.938*1000.;//MeV
 // Double_t massCarbon = 11.17802*1000.;//MeV

  std::vector<TrackInfo_t>::iterator iteTrackInfo = trackInfoVec.begin();
   for(;iteTrackInfo != trackInfoVec.end();iteTrackInfo++){  

     massPart = (*iteTrackInfo).massPart;
     pdg_t = (*iteTrackInfo).PDGcode; 
     Pxx   = (*iteTrackInfo).Mom.X();
     Pyy   = (*iteTrackInfo).Mom.Y();
     Pzz   = (*iteTrackInfo).Mom.Z();
     PartE = (*iteTrackInfo).Mom.E();

     PT    = sqrt(Pxx*Pxx+Pyy*Pyy);

     Ecm   = sqrt(massProton*massProton+massProton*massProton+2.*BeamEnergy*massProton);

     beta  = sqrt(BeamEnergy*BeamEnergy-massProton*massProton)/(BeamEnergy+massProton);
     gamma = sqrt(1.-beta*beta);
     PL    = gamma*(Pzz-beta*PartE);    
     XF    = 2.*PL/Ecm;
     
     g4Proddata.PDG[partNum] = pdg_t;
     g4Proddata.InterType[partNum] = (*iteTrackInfo).interType;
     g4Proddata.X[partNum][0]= (*iteTrackInfo).Pos.X();
     g4Proddata.X[partNum][1]= (*iteTrackInfo).Pos.Y();
     g4Proddata.X[partNum][2]= (*iteTrackInfo).Pos.Z();
     g4Proddata.P[partNum][0]= Pxx;
     g4Proddata.P[partNum][1]= Pyy;
     g4Proddata.P[partNum][2]= Pzz;
     g4Proddata.P[partNum][3]= PartE;
     g4Proddata.PT[partNum]  = PT;
     g4Proddata.XF[partNum]  = XF;
     

     //Types of particles counter:
     if(pdg_t ==  211)NPiPlus++;
     if(pdg_t == -211)NPiMinus++; 
     if(pdg_t ==  111)NPi0++; 
     if(pdg_t ==  321)NKPlus++;
     if(pdg_t == -321)NPiMinus++;  
     if(pdg_t ==  310)NK0S++; 
     if(pdg_t ==  130)NK0L++;  
     if(pdg_t == 2212)NProtons++; 
     if(pdg_t ==-2212)NAProtons++;
     if(pdg_t == 2112)NNeutrons++;
     if(pdg_t ==-2112)NANeutrons++;
     if((3100<pdg_t && pdg_t<3400)||(-3400<pdg_t && pdg_t<-3100)||(pdg_t==221))NOthers++;
     partNum++; 
 }
    
   g4Proddata.PartTypes[0]=NPiPlus;
   g4Proddata.PartTypes[1]=NPiMinus;
   g4Proddata.PartTypes[2]=NPi0;
   g4Proddata.PartTypes[3]=NKPlus;
   g4Proddata.PartTypes[4]=NKMinus;
   g4Proddata.PartTypes[5]=NK0S;
   g4Proddata.PartTypes[6]=NK0L;
   g4Proddata.PartTypes[7]=NProtons;
   g4Proddata.PartTypes[8]=NAProtons;
   g4Proddata.PartTypes[9]=NNeutrons;
   g4Proddata.PartTypes[10]=NANeutrons;
   g4Proddata.PartTypes[11]=NOthers;

 if (g4Proddata.NPart>0)WriteNtuple();
}

void NA49Analysis::WriteNtuple(){
    
    ProdTree->Fill();
    
}
