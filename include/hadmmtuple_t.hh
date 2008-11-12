//----------------------------------------------------------------------
// hadmmtuple_t.hh
//
// ADM, July 2005
// This is a class that defines the hadmmtuple_t object that is used to store 
// the g4numi output in a root tree.
//
// $Id: hadmmtuple_t.hh,v 1.6 2008/11/12 00:21:40 loiacono Exp $
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

  void Clear();

   // the following variables are placed in the root tree
  UShort_t evtno; 
  Float_t muvx;
  Float_t muvy;
  //Float_t muvz;
  Float_t mupx;
  Float_t mupy;
  Float_t mupz;
  Double_t muweight;
  Float_t tpx;
  Float_t tpy;
  Float_t tpz;
  Short_t tpptype;
  Double_t nimpwt;
  Short_t pptype;
  Short_t ptype;
  Float_t mmxpos[3];
  Float_t mmypos[3];
  //Float_t mmzpos[3];
  Float_t mmpx[3];
  Float_t mmpy[3];
  Float_t mmpz[3];
  Short_t cell[3];
  //Float_t mmxpos_Edep[3];
  //Float_t mmypos_Edep[3];
  //Float_t mmzpos_Edep[3];
  Float_t mmpx_Edep[3];
  Float_t mmpy_Edep[3];
  Float_t mmpz_Edep[3];

  /*Float_t hmmenergy;
  Float_t hmmxpos;
  Float_t hmmypos;
  Float_t hmmzpos;
  Float_t hmmpx;
  Float_t hmmpy;
  Float_t hmmpz; 
  Short_t run;
  Float_t mtgthsig;
  Float_t mtgtvsig;
  Float_t mtgthpos;
  Float_t mtgtvpos;
  */  
  
 private:
  ClassDef(hadmmtuple_t , 1)
    
    };
#endif // admmtuple_t_h

