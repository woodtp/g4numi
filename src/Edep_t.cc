//----------------------------------------------------------------------
// Sets the relevant memebers of the data class for storing the
// MC data for the Hadron and Muon Monitors.
//
// $Id: Edep_t.cc,v 1.1.2.1 2010/08/19 19:50:54 minervacvs Exp $
//----------------------------------------------------------------------


#include "globals.hh"
#include "G4ios.hh"

#include "Edep_t.hh"

ClassImp(Edep_t)

//-----------------------------------------------------------------------------------
   Edep_t::Edep_t()
      :nTracks(0),
       wghtedNTracks(0.0),
       sumEdepWghts(0.0),
       sumWghtdEdep(0.0),
       sumWghtdEdep2(0.0),
       trackVec()
{
   
}

//-----------------------------------------------------------------------------------
Edep_t::~Edep_t()
{
  // Edep_t destructor
}

//------------------------------------------------------------------------------------------
void Edep_t::ClearVector()
{
   trackVec.clear();
}
