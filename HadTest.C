#include <TH1.h>
#include <TBranch.h>
#include <iostream>
#include "TMath.h" 
#include <TString.h>
#include <iostream>

void HadTest(TString fname){
 
  TFile* f = new TFile(fname);

  TTree* tRun  = (TTree*)f->Get("RunInfo");
  TTree* tEvts = (TTree*)f->Get("nTarget");

    // Declaration of leaf types
   Double_t        inelasticXS;
   Double_t        elasticXS;
   Double_t        tickness;
   Double_t        radius;
   Double_t        density;
   string          *material;
   Double_t        enerPrimGen;
   string          *PartName;
   Int_t           numberEvts;
   Int_t           npart;
   Int_t           pdg[150];   //[NPart]
   Int_t           intType[150];   //[NPart]
   Int_t           part_types[12];
   Double_t        x[150][3];   //[NPart]
   Double_t        p[150][4];   //[NPart]
   Double_t        xf[150];   //[NPart]
   Double_t        pt[150];   //[NPart]

   // List of branches
   TBranch        *b_inelXS;   //!
   TBranch        *b_elXS;   //!
   TBranch        *b_tickness;   //!
   TBranch        *b_radius;   //!
   TBranch        *b_density;   //!
   TBranch        *b_material;   //!
   TBranch        *b_enerPrimGen;   //!
   TBranch        *b_PartName;   //!
   TBranch        *b_numberEvts;   //!
   TBranch        *b_NPart;   //!
   TBranch        *b_pdg;   //!
   TBranch        *b_intType;   //!
   TBranch        *b_PartTypes;   //!
   TBranch        *b_x;   //!
   TBranch        *b_p;   //!
   TBranch        *b_xf;   //!
   TBranch        *b_pt;   //!


   //set branches
   tRun->SetBranchAddress("inelasticXS", &inelasticXS, &b_inelXS);
   tRun->SetBranchAddress("elasticXS", &elasticXS, &b_elXS);
   tRun->SetBranchAddress("tickness", &tickness, &b_tickness);
   tRun->SetBranchAddress("radius", &radius, &b_radius);
   tRun->SetBranchAddress("density", &density, &b_density);
   tRun->SetBranchAddress("material", &material, &b_material);
   tRun->SetBranchAddress("enerPrimGen", &enerPrimGen, &b_enerPrimGen);
   tRun->SetBranchAddress("PartName", &PartName, &b_PartName);
   tRun->SetBranchAddress("numberEvts", &numberEvts, &b_numberEvts);
   tEvts->SetBranchAddress("npart", &npart, &b_NPart);
   tEvts->SetBranchAddress("pdg", pdg, &b_pdg);
   tEvts->SetBranchAddress("intType", intType, &b_intType);
   tEvts->SetBranchAddress("part_types", part_types, &b_PartTypes);
   tEvts->SetBranchAddress("x", x, &b_x);
   tEvts->SetBranchAddress("p", p, &b_p);
   tEvts->SetBranchAddress("xf", xf, &b_xf);
   tEvts->SetBranchAddress("pt", pt, &b_pt);


//Get info:
  Int_t nentriesRun  = (Int_t)tRun->GetEntries();
  Int_t nentriesEvts = (Int_t)tEvts->GetEntries();
 
  tRun->GetEntry(0);
  tRun->Show(0);  
 
  tEvts->GetEntry(1);
  tEvts->Show(1);

}
