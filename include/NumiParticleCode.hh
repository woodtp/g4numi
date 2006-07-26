#ifndef NUMIPARTICLECODE_H
#define NUMIPARTICLECODE_H

#include "globals.hh"
#include "Rtypes.h"
#include <string>

namespace NumiParticleCode
{
  typedef enum ENumiParticleCode
    { 
      kMuonPlus             = 5,
      kMuonMinus            = 6,
      kPion0                = 7,
      kPionPlus             = 8,
      kPionMinus            = 9,
      kKaon0L               = 10,
      kKaonPlus             = 11,
      kKaonMinus            = 12,
      kNeutron              = 13,
      kProton               = 14,
      kAntiProton           = 15,
      kKaon0S               = 16,
      kEta                  = 17,
      kLambda               = 18,
      kSigmaPlus            = 19,
      kSigma0               = 20,
      kSigmaMinus           = 21,
      kXi0                  = 22,
      kXiMinus              = 23,
      kOmegaMinus           = 24,
      kAntiNeutron          = 25,
      kAntiLambda           = 26,
      kAntiSigmaMinus       = 27,
      kAntiSigma0           = 28,
      kAntiSigmaPlus        = 29,
      kAntiXi0              = 30,
      kAntiXiMinus          = 31,
      kElectronAntiNeutrino = 52,
      kElectronNeutrino     = 53,
      kMuonAntiNeutrino     = 55,
      kMuonNeutrino         = 56,
      kOther                = 99 
    } NumiParticleCode_t;

  G4String NumiParticleCode::AsString(NumiParticleCode_t pCode);
  G4int NumiParticleCode::AsInt(NumiParticleCode_t pCode);
  NumiParticleCode_t NumiParticleCode::IntToEnum(G4int particleInt);
  NumiParticleCode_t NumiParticleCode::StringToEnum(G4String particleName);

}

#endif
