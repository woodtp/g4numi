
#ifndef PRODTUPLE_T_HH
#define PRODTUPLE_T_HH

#include "TROOT.h"          // Top level (or root) structure for all classes
#include "TObject.h"

class ProdTuple_t 
{
  
 public:
  // a constructor and a destructor
  ProdTuple_t();
  virtual ~ProdTuple_t();
  
  
  // the following variables are placed in the root tree
  Double_t xpos;
  Double_t ypos;
  Double_t zpos;
  Double_t xmom;
  Double_t ymom;
  Double_t zmom;
  Double_t ener;
  Int_t    pdgcode;

 private:
  ClassDef(ProdTuple_t , 1)
    
    };
#endif 

