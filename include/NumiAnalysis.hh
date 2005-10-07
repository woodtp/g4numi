//
// NumiAnalysis.hh
//
// Modified Jul 2005 by A. MArino to make data_t and hadmmtuple_t classes

#ifndef NUMIANALYSIS_HH
#define NUMIANALYSIS_HH

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

#include "NumiDataInput.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "NumiTrajectory.hh"
#include "G4TrajectoryContainer.hh"

class G4Track;
class data_t;
class hadmmtuple_t;

class NumiAnalysis
{
public:

  ~NumiAnalysis();

  void book();
  
  void finish();
  
  void FillNeutrinoNtuple(const G4Track& track);
  void FillHadmmNtuple(const G4Track& track);
  G4double GetWeight(const G4Track& nutrack, 
		     G4double enuzr, 
		     G4ThreeVector vertex_r,
		     G4double gamma,
		     G4ThreeVector beta_vec,
		     G4double theta_pardet, 
		     G4double x_det, 
		     G4double y_det, 
		     G4double z_det);
  G4double GetNuEnergy(G4double Parent_mass, G4double gamma, G4double beta, G4double theta);
  G4double GetTheta(G4ThreeVector vertex_r,G4ThreeVector momentum,G4double x_det,G4double y_det,G4double z_det);
  
  G4int GetParticleCode(G4String);
  G4String GetParticleName(G4int);
  NumiTrajectory* GetParentTrajectory(G4int parentID);
  
  static NumiAnalysis* getInstance();
private:

  NumiAnalysis();

  static NumiAnalysis* instance;

  G4double x;
  G4double y;
  G4double z;

  G4double noProtons;
  char asciiFileName[30], nuNtupleFileName[30], hadmmNtupleFileName[30];

  NumiDataInput* NumiData;

  TFile* hadmmNtuple;
  TFile* nuNtuple;
 
  TTree* tree;
  TTree* hadmmtree;

  data_t *g4data;
  hadmmtuple_t *g4hmmdata;
  /*
  typedef struct {
    Int_t run;        //
    Int_t evtno;
    Double_t beamHWidth;
    Double_t beamVWidth;
    Double_t beamX;
    Double_t beamY;
    Double_t protonX;
    Double_t protonY;
    Double_t nuTarZ;
    Double_t hornCurrent;
    Double_t Ndxdz;
    Double_t Ndydz;
    Double_t Npz;
    Double_t Nenergy;
    Double_t NdxdzNea;
    Double_t NdydzNea;
    Double_t NenergyN[10];
    Double_t NWtNear[10];
    Double_t NdxdzFar;
    Double_t NdydzFar;
    Double_t NenergyF[10];
    Double_t NWtFar[10];
    Int_t    Norig;
    Int_t    Ndecay;
    Int_t    Ntype;
    Double_t Vx;
    Double_t Vy;
    Double_t Vz;
    Double_t pdPx;
    Double_t pdPy;
    Double_t pdPz;
    Double_t ppdxdz;
    Double_t ppdydz;
    Double_t pppz;
    Double_t ppenergy;
    Double_t ppmedium;
    Int_t    ptype;
    Double_t ppvx;
    Double_t ppvy;
    Double_t ppvz;
    Double_t muparpx;
    Double_t muparpy;
    Double_t muparpz;
    Double_t mupare;
    Double_t Necm;
    Double_t Nimpwt;
    Double_t xpoint;
    Double_t ypoint;
    Double_t zpoint;
    Double_t tvx;
    Double_t tvy;
    Double_t tvz;
    Double_t tpx;
    Double_t tpy;
    Double_t tpz;
    Int_t tptype;
    Int_t tgen;
    Double_t trkx[10];
    Double_t trky[10];
    Double_t trkz[10];
    Double_t trkpx[10];
    Double_t trkpy[10];
    Double_t trkpz[10];
  } data_t;
  data_t g4data;
  //  std::ofstream asciiFile;

  typedef struct{
    Int_t run;        //
    Int_t evtno; 
    Double_t mtgthpos;
    Double_t mtgtvpos;
    Double_t mtgthsig;
    Double_t mtgtvsig;
    Int_t ptype;
    Double_t hmmenergy;
    Double_t hmmxpos;
    Double_t hmmypos;
    Double_t hmmzpos;
    Double_t hmmpx;
    Double_t hmmpy;
    Double_t hmmpz;
  } hadmmtuple_t;
  hadmmtuple_t g4hmmdata;
*/  
};

#endif 
