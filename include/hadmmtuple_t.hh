//----------------------------------------------------------------------
// hadmmtuple_t.hh
//
// ADM, July 2005
// This is a class that defines the hadmmtuple_t object that is used to store 
// the g4numi output in a root tree.
//
// $Id: hadmmtuple_t.hh,v 1.5 2008/09/17 17:36:23 loiacono Exp $
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
  Int_t run;       
  Int_t evtno; 
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
  Int_t    tpptype;
  Double_t nimpwt;
  Double_t mtgthpos;
  Double_t mtgtvpos;
  Double_t mtgthsig;
  Double_t mtgtvsig;
  Int_t pptype;
  Int_t ptype;
  Double_t hmmenergy;
  Double_t hmmxpos;
  Double_t hmmypos;
  Double_t hmmzpos;
  Double_t hmmpx;
  Double_t hmmpy;
  Double_t hmmpz;
  Double_t mmxpos[3];
  Double_t mmypos[3];
  Double_t mmzpos[3];
  Double_t mmpx[3];
  Double_t mmpy[3];
  Double_t mmpz[3];
  Double_t cell[3];
  Double_t mmxpos_Edep[3];
  Double_t mmypos_Edep[3];
  Double_t mmzpos_Edep[3];
  Double_t mmpx_Edep[3];
  Double_t mmpy_Edep[3];
  Double_t mmpz_Edep[3];
  
 private:
  ClassDef(hadmmtuple_t , 1)
    
    };
#endif // admmtuple_t_h

