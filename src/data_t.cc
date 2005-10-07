//
// data_t.cc
//
//  ADM, July 2005
//  This is a class that defines the data_t object that is used to store the g4numi output
//   in a root tree.
//------------------

#include "data_t.hh"
ClassImp(data_t)

data_t::data_t():
     run(-1),       
     evtno(-1),
     beamHWidth(-1.e4),
     beamVWidth(-1.e4),
     beamX(-1.e4),
     beamY(-1.e4),
     protonX(-1.e4),
     protonY(-1.e4),
     protonZ(-1.e4),
     protonPx(-1.e4),
     protonPy(-1.e4),
     protonPz(-1.e4),
     nuTarZ(-1.e4),
     hornCurrent(-1.e4),
     Ndxdz(-1.e4),
     Ndydz(-1.e4),
     Npz(-1.e4),
     Nenergy(-1.e4),
     Norig(-1),
     Ndecay(-1),
     Ntype(-1),
     Vx(-1.e4),
     Vy(-1.e4),
     Vz(-1.e4),
     pdPx(-1.e4),
     pdPy(-1.e4),
     pdPz(-1.e4),
     ppdxdz(-1.e4),
     ppdydz(-1.e4),
     pppz(-1.e4),
     ppenergy(-1.e4),
     ppmedium(-1.e4),
     ptype(-1),
     ppvx(-1.e4),
     ppvy(-1.e4),
     ppvz(-1.e4),
     muparpx(-1.e4),
     muparpy(-1.e4),
     muparpz(-1.e4),
     mupare(-1.e4),
     Necm(-1.e4),
     Nimpwt(-1.e4),
     xpoint(-1.e4),
     ypoint(-1.e4),
     zpoint(-1.e4),
     tvx(-1.e4),
     tvy(-1.e4),
     tvz(-1.e4),
     tpx(-1.e4),
     tpy(-1.e4),
     tpz(-1.e4),
     tptype(-1),
     tgen(-1)
{
  
  // data_t record constructor
  for ( Int_t index=0; index<9; ++index){
    NdxdzNear[index]= -1.e4;
    NdydzNear[index]= -1.e4;
    NenergyN[index]= -1.e4;
    NWtNear[index]= -1.e4 ;
  }
  for ( Int_t index=0; index<2; ++index){
    NdxdzFar[index]= -1.e4;
    NdydzFar[index]= -1.e4;
    NenergyF[index]= -1.e4;
    NWtFar[index]= -1.e4;
  }
  for ( Int_t index=0; index<10; ++index){
    trkx[index]=-1.e4;
    trky[index]=-1.e4;
    trkz[index]=-1.e4;
    trkpx[index]=-1.e4;
    trkpy[index]=-1.e4;
    trkpz[index]=-1.e4;
  }
}


//----------------------------------------------------------------------
data_t::~data_t()
{
  // data_t destructor
}

