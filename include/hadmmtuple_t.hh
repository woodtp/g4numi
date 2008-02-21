//----------------------------------------------------------------------
// hadmmtuple_t.hh
//
// ADM, July 2005
// This is a class that defines the hadmmtuple_t object that is used to store 
// the g4numi output in a root tree.
//
// $Id: hadmmtuple_t.hh,v 1.3 2008/02/21 15:43:03 koskinen Exp $
//----------------------------------------------------------------------

#ifndef HADMMTUPLE_T_HH
#define HADMMTUPLE_T_HH

#include "TROOT.h"          // Top level (or root) structure for all classes
#include "TObject.h"

class  hadmmtuple_t 
{
  
 public:
  // a constructor and a destructor
  hadmmtuple_t();
  virtual ~hadmmtuple_t();
  
  
  // the following variables are placed in the root tree
  Int_t    run;       
  Int_t    evtno; 
  Double_t mtgthpos;
  Double_t mtgtvpos;
  Double_t mtgthsig;
  Double_t mtgtvsig;
  Int_t    ptype;
  Double_t hmmenergy;
  Double_t hmmxpos;
  Double_t hmmypos;
  Double_t hmmzpos;
  Double_t hmmpx;
  Double_t hmmpy;
  Double_t hmmpz;
  Double_t mmxpos[4];
  Double_t mmypos[4];
  Double_t mmzpos[4];
  Double_t mmpx[4];
  Double_t mmpy[4];
  Double_t mmpz[4];
  Int_t    cell[4];
  
 private:
  ClassDef(hadmmtuple_t , 1)
    
    };
#endif // admmtuple_t_h

