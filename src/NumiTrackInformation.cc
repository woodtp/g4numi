//
// NumiTrackInformation.cc
//

#include "NumiTrackInformation.hh"
#include "G4ios.hh"

G4Allocator<NumiTrackInformation> aTrackInformationAllocator;

NumiTrackInformation::NumiTrackInformation()
{
  decay_code = 0;
  tgen = 0;
}


NumiTrackInformation::NumiTrackInformation(const NumiTrackInformation* aTrackInfo)
{
  decay_code = aTrackInfo->decay_code;
  tgen = aTrackInfo->tgen;
}

NumiTrackInformation::~NumiTrackInformation(){}

void NumiTrackInformation::Print() const
{
    G4cout 
     << "Decay code = " << decay_code << G4endl;
    G4cout 
     << "tgen = " << tgen << G4endl;
}

