
#ifndef NumiPrimaryGeneratorAction_h
#define NumiPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class NumiPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    NumiPrimaryGeneratorAction();
    ~NumiPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4ParticleGun* particleGun;
};

#endif


