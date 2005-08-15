#ifndef NumiRunActionMessenger_h
#define NumiRunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ios.hh"

class NumiRunAction;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

class NumiRunActionMessenger: public G4UImessenger
{
public:
  NumiRunActionMessenger(NumiRunAction* );
  ~NumiRunActionMessenger();
  
  void SetNewValue(G4UIcommand* ,G4String );
 
private:
  NumiRunAction*            runAction;
  G4UIdirectory*            RndmDir;
  G4UIcmdWithAString*       readRndmCmd;  
  G4UIcmdWithoutParameter*  showRndmCmd;
  G4UIcmdWithAnInteger*     setRndmSeedCmd;
  G4UIdirectory*            NumiRunDir;
  G4UIcmdWithAnInteger*     setRunID;

};

#endif
