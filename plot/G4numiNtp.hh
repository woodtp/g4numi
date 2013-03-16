#ifndef G4NUMINTP_H
#define G4NUMINTP_H 


#include <TH1.h>
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <math.h>
#include <TLeaf.h>

class G4numiNtp {

public:
  G4numiNtp(void);
  ~G4numiNtp(void);

      Int_t run;
      Int_t evtno;
      Double_t Ndxdz;
      Double_t Ndydz;
      Double_t Npz;
      Double_t Nenergy;

      Int_t Norig;
      Int_t Ndecay;
      Int_t Ntype;
      Double_t Vx;
      Double_t Vy;
      Double_t Vz;
      Double_t pdpx;
      Double_t pdpy;
      Double_t pdpz;
      Double_t ppdxdz;
      Double_t ppdydz;
      Double_t pppz;
      Double_t ppenergy;
      Int_t ppmedium;
      Int_t ptype;
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
      Int_t tgptype;
      Double_t tgppx;
      Double_t tgppy;
      Double_t tgppz;
      Double_t tprivx;
      Double_t tprivy;
      Double_t tprivz;
      Double_t beamx;
      Double_t beamy;
      Double_t beamz;
      Double_t beampx;
      Double_t beampy;
      Double_t beampz;

      Double_t NdxdzNear;
      Double_t NdydzNear;
      Double_t NenergyN;
      Double_t NWtNear;
      Double_t NWtNear_one;
      Double_t NdxdzFar;
      Double_t NdydzFar;
      Double_t NenergyF;
      Double_t NWtFar;
  
  //void SetTree(TTree &t, Int_t entry);
  void SetTree(TTree &t);
  
private:
  Int_t entry;
  Double_t NdxdzNear_array[11];
  Double_t NdydzNear_array[11];
  Double_t NenergyN_array[11];
  Double_t NWtNear_array[11];
  Double_t NdxdzFar_array[2];
  Double_t NdydzFar_array[2];
  Double_t NenergyF_array[2];
  Double_t NWtFar_array[2];
  TLeaf *l_NdxdzNear;
  TLeaf *l_NdydzNear;
  TLeaf *l_NenergyN;
  TLeaf *l_NWtNear;
  TLeaf *l_NdxdzFar;
  TLeaf *l_NdydzFar;
  TLeaf *l_NenergyF;
  TLeaf *l_NWtFar;

      TLeaf *l_run;
      TLeaf *l_evtno;
      TLeaf *l_Ndxdz;
      TLeaf *l_Ndydz;
      TLeaf *l_Npz;
      TLeaf *l_Nenergy;
      TLeaf *l_Norig;
      TLeaf *l_Ndecay;
      TLeaf *l_Ntype;
      TLeaf *l_Vx;
      TLeaf *l_Vy;
      TLeaf *l_Vz;
      TLeaf *l_pdpx;
      TLeaf *l_pdpy;
      TLeaf *l_pdpz;
      TLeaf *l_ppdxdz;
      TLeaf *l_ppdydz;
      TLeaf *l_pppz;
      TLeaf *l_ppenergy;
      TLeaf *l_ppmedium;
      TLeaf *l_ptype;
      TLeaf *l_ppvx;
      TLeaf *l_ppvy;
      TLeaf *l_ppvz;
      TLeaf *l_muparpx;
      TLeaf *l_muparpy;
      TLeaf *l_muparpz;
      TLeaf *l_mupare;
      TLeaf *l_Necm;
      TLeaf *l_Nimpwt;
      TLeaf *l_xpoint;
      TLeaf *l_ypoint;
      TLeaf *l_zpoint;
      TLeaf *l_tvx;
      TLeaf *l_tvy;
      TLeaf *l_tvz;
      TLeaf *l_tpx;
      TLeaf *l_tpy;
      TLeaf *l_tpz;
      TLeaf *l_tptype;
      TLeaf *l_tgen;
      TLeaf *l_tgptype;
      TLeaf *l_tgppx;
      TLeaf *l_tgppy;
      TLeaf *l_tgppz;
      TLeaf *l_tprivx;
      TLeaf *l_tprivy;
      TLeaf *l_tprivz;
      TLeaf *l_beamx;
      TLeaf *l_beamy;
      TLeaf *l_beamz;
      TLeaf *l_beampx;
      TLeaf *l_beampy;
      TLeaf *l_beampz;


};

#endif
