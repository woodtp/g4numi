//----------------------------------------------------------------------
// Sets the relevant memebers of the data class for storing the
// MC data for the Hadron and Muon Monitors.
//
// $Id: hadmmtuple_t.cc,v 1.6 2008/09/17 17:36:23 loiacono Exp $
//----------------------------------------------------------------------

#include "hadmmtuple_t.hh"

ClassImp(hadmmtuple_t)

hadmmtuple_t::hadmmtuple_t():
  run(-81000),       
  evtno(-81000), 
  muvx(-81000.),
  muvy(-81000.),
  muvz(-81000.),
  mupx(-81000.),
  mupy(-81000.),
  mupz(-81000.),
  muweight(-81000.),
  tpx(-81000.),
  tpy(-81000.),
  tpz(-81000.),
  tpptype(-81000),
  nimpwt(-81000.),
  mtgthpos(-81000.),
  mtgtvpos(-81000.),
  mtgthsig(-81000.),
  mtgtvsig(-81000.),
  pptype(-81000),
  ptype(-81000),
  hmmenergy(-81000.),
  hmmxpos(-81000.),
  hmmypos(-81000.),
  hmmzpos(-81000.),
  hmmpx(-81000.),
  hmmpy(-81000.),
  hmmpz(-81000.)
{
  for(Int_t index=0;index<3;++index)
    {
      mmxpos[index]=-81000.;
      mmypos[index]=-81000.;
      mmzpos[index]=-81000.;
      mmpx[index]=-81000.;
      mmpy[index]=-81000.;
      mmpz[index]=-81000.;
      cell[index]=-81000;

      mmxpos_Edep[index]=-81000.;
      mmypos_Edep[index]=-81000.;
      mmzpos_Edep[index]=-81000.;
      mmpx_Edep[index]=-81000.;
      mmpy_Edep[index]=-81000.;
      mmpz_Edep[index]=-81000.;
    }

}
  //-----------------------------------------------------------------------------------
hadmmtuple_t::~hadmmtuple_t()
{
  // hadmmtuple_t destructor
}

