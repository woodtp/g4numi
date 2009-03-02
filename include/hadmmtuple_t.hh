//----------------------------------------------------------------------
// hadmmtuple_t.hh
//
// ADM, July 2005
// This is a class that defines the hadmmtuple_t object that is used to store 
// the g4numi output in a root tree.
//
// $Id: hadmmtuple_t.hh,v 1.7 2009/03/02 03:35:12 loiacono Exp $
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
  Float_t tvx;
  Float_t tvy;
  Float_t tvz;
  Float_t tpx;
  Float_t tpy;
  Float_t tpz;
  Short_t tpptype;
  Double_t nimpwt;
   Float_t ppvx;
   Float_t ppvy;
   Float_t ppvz;
   Float_t pdvx;
   Float_t pdvy;
   Float_t pdvz;
   Float_t pdpx;
  Float_t pdpy;
   Float_t pdpz;
   Short_t pptype;
   Short_t ptype;
    Short_t ppmedium;
    Short_t pgen;
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

