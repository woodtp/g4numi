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

   Int_t    evtno;
   Double_t muvx;
   Double_t muvy;
   Double_t muvz;
   Double_t mupx;
   Double_t mupy;
   Double_t mupz;
   Double_t muweight;
   Double_t tpx;
   Double_t tpy;
   Double_t tpz;
   Int_t    tptype;
   Double_t nimpwt;
   Int_t ptype;

   UShort_t evtnoS;
   Float_t muvxF;
   Float_t muvyF;
   Float_t muvzF;
   Float_t mupxF;
   Float_t mupyF;
   Float_t mupzF;
   Float_t tpxF;
   Float_t tpyF;
   Float_t tpzF;
   Short_t tptypeS;
   Short_t ptypeS;


        
private:

   bool SetBranches();

   void SetBranchesStatus();

   bool SetBranch(void *var, const std::string &name);



   
   std::map<std::string, TBranch*> fBranch;
   
   TTree *ftree;
};

#endif
