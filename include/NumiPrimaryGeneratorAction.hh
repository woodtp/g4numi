
#ifndef NumiPrimaryGeneratorAction_h
#define NumiPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class NumiDataInput;

class NumiPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    NumiPrimaryGeneratorAction();
    ~NumiPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    NumiDataInput* ND;
    G4ParticleGun* particleGun;
    G4double meanx;
    G4double meany;
    G4double meanz;
};

#endif


