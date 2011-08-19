
#include <vector>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <stdlib.h>
#include <map.h>

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

//g4na49
#include "ProdTuple_t.hh"
#include "NA49Analysis.hh"

using namespace std;

NA49Analysis* NA49Analysis::instance = 0;

//------------------------------------------------------------------------------------
NA49Analysis::NA49Analysis()
{
   G4cout << "NA49Analysis" << G4endl;
   g4Proddata = new ProdTuple_t();

#ifdef G4ANALYSIS_USE
#endif
 }
//------------------------------------------------------------------------------------
NA49Analysis::~NA49Analysis()
{ 


#ifdef G4ANALYSIS_USE
  // delete things
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
    sprintf(NtupleFileName,"NA49_%04d.root",pRunManager->GetCurrentRun()->GetRunID());
    FileNtuple = new TFile(NtupleFileName, "RECREATE","pion from p+C ntuple");
    FileNtuple->cd();
    ProdTree = new TTree("PionInfo","g4na49 Pion from p+C");
    ProdTree->Branch("piondata", "ProdTuple_t",&g4Proddata,32000,1);
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

void NA49Analysis::FillNtuple(const G4Track& track)
{
 G4RunManager* pRunManager = G4RunManager::GetRunManager();
 G4ParticleDefinition* pd = track.GetDefinition();
 // if(pd != G4Electron::Electron())G4cout<<"Particle "<<pd->GetParticleName()<<G4endl;
 if((pd == G4PionPlus::PionPlus())||(pd == G4PionMinus::PionMinus())){
 g4Proddata->xpos    = track.GetPosition()[0]/cm;
 g4Proddata->ypos    = track.GetPosition()[1]/cm;
 g4Proddata->zpos    = track.GetPosition()[2]/cm;
 g4Proddata->xmom    = track.GetMomentum()[0]/GeV;
 g4Proddata->ymom    = track.GetMomentum()[1]/GeV;
 g4Proddata->zmom    = track.GetMomentum()[2]/GeV;
 g4Proddata->ener    = track.GetTotalEnergy()/GeV;
 g4Proddata->pdgcode = pd->GetPDGEncoding();
 WriteNtuple();
 }
}

void NA49Analysis::WriteNtuple(){

    ProdTree->Fill();

    g4Proddata->xpos    = -10000;
    g4Proddata->ypos    = -10000;
    g4Proddata->zpos    = -10000;
    g4Proddata->xmom    = -10000;
    g4Proddata->ymom    = -10000;
    g4Proddata->zmom    = -10000;
    g4Proddata->ener    = -10000;
    g4Proddata->pdgcode = -10000;
}
