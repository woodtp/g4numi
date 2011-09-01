
#ifndef PRODTUPLE_T_HH
#define PRODTUPLE_T_HH

#include "TROOT.h"        
#include "TObject.h"
#include "Rtypes.h"

const Int_t maxPart = 50;
//Particle types:
// 0: pi+    1: pi-   2: pi0
// 3: K+     4: K-    5: K0S   6: K0L
// 7: p      8: antip 9: n    10: antin
//11: lambda, sigma, omega, eta, xi //are these all particles for here?

class ProdTuple_t 
{
  
 public:
  // a constructor and a destructor
  ProdTuple_t();
  virtual ~ProdTuple_t();
  
  // the following variables are placed in the root tree
  Int_t PDG[maxPart];
  Double_t X[maxPart][3];//Initial position of the track
  Double_t P[maxPart][4];//4-momentum. P[maxPart][3] is the energy
  Double_t XF[maxPart];
  Double_t PT[maxPart];
 
  Int_t NPart;
  Int_t PartTypes[12];
 
 private:
  ClassDef(ProdTuple_t , 1)
    
    };
#endif 

