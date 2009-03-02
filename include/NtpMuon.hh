#ifndef NTPMUON_H
#define NTPMUON_H


//C++
#include <map>

// ROOT
#include "Rtypes.h" // include Float_t, Double_t, etc definitions

#include "TTree.h"
class TTree;
class TBranch;
//class Rtypes;

//Local
class NtpMuon
{
public:
   
   NtpMuon();
   ~NtpMuon();

   void Clear();

   bool SetTree(TTree* tree);

   int GetNEntries();
   
   void GetEntry(unsigned int entry);

   void GetBranchEntry(std::string name, unsigned int entry);
   
public:

   Int_t evtno;
   Double_t muvx;
   Double_t muvy;
   Double_t muvz;
   Double_t mupx;
   Double_t mupy;
   Double_t mupz;
   Double_t muweight;
   Double_t tvx;
   Double_t tvy;
   Double_t tvz;
   Double_t tpx;
   Double_t tpy;
   Double_t tpz;
   Int_t tptype;
   Double_t nimpwt;
   Double_t ppvx;
   Double_t ppvy;
   Double_t ppvz;
   Double_t pdvx;
   Double_t pdvy;
   Double_t pdvz;
   Double_t pdpx;
   Double_t pdpy;
   Double_t pdpz;
   Int_t ptype;
   Int_t ppmedium;
   Int_t pgen;

   UShort_t evtnoS;
   Float_t muvxF;
   Float_t muvyF;
   Float_t muvzF;
   Float_t mupxF;
   Float_t mupyF;
   Float_t mupzF;
   Double_t muweightD;
   Float_t tvxF;
   Float_t tvyF;
   Float_t tvzF;
   Float_t tpxF;
   Float_t tpyF;
   Float_t tpzF;
   Short_t tptypeS;
   Double_t nimpwtD;
   Float_t ppvxF;
   Float_t ppvyF;
   Float_t ppvzF;
   Float_t pdvxF;
   Float_t pdvyF;
   Float_t pdvzF;
   Float_t pdpxF;
   Float_t pdpyF;
   Float_t pdpzF;
   Short_t ptypeS;
   Short_t ppmediumS;
   Short_t pgenS;
 
        
private:

   bool SetBranches();

   void SetBranchesStatus();

   bool SetBranch(void *var, const std::string &name);



   
   std::map<std::string, TBranch*> fBranch;
   
   TTree *ftree;
};

#endif
