
#ifndef NA49StackingMessenger_h
#define NA49StackingMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class NA49StackingAction;
class G4UIcmdWithABool;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class NA49StackingMessenger: public G4UImessenger
{
public:

  NA49StackingMessenger(NA49StackingAction*);
  virtual ~NA49StackingMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:
    
  NA49StackingAction*     stackAction;
    
  G4UIcmdWithABool*   killCmd;
  G4UIcmdWithAString* kCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
