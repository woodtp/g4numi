#ifndef NumiPrimaryMessenger_h
#define NumiPrimaryMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ios.hh"
#include "NumiPrimaryGeneratorAction.hh"

class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithoutParameter;

class NumiPrimaryMessenger: public G4UImessenger
{
public:
  NumiPrimaryMessenger(NumiPrimaryGeneratorAction* );
  ~NumiPrimaryMessenger();
  
  void SetNewValue(G4UIcommand* ,G4String );
 
private:
  NumiPrimaryGeneratorAction*   PrimaryAction;
  G4ParticleTable*          particleTable;


  G4UIdirectory*            BeamDir;

  
  G4UIcmdWithADouble* setCosX;
  G4UIcmdWithADouble* setCosY;
  G4UIcmdWithADoubleAndUnit* setZ;
  G4UIcmdWithADoubleAndUnit* setMomentum;
  G4UIcmdWithADoubleAndUnit* setSpread;
  G4UIcmdWithADoubleAndUnit* setDivergence;
  G4UIcmdWithADoubleAndUnit* setInnerR;
  G4UIcmdWithADoubleAndUnit* setOuterR;
  G4UIcmdWithAString*        setParticle;  


};

#endif
