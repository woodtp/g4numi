
#ifndef NA49StackingAction_h
#define NA49StackingAction_h 1

#include "G4UserStackingAction.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "NA49Analysis.hh"

class NA49StackingMessenger;
class G4Track;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class NA49StackingAction : public G4UserStackingAction
{
public:

  NA49StackingAction();
  virtual ~NA49StackingAction();
   
  void SetKillStatus(G4bool value)    {killSecondary = value;};
  void SetKill(const G4String& name)  {pname = name;};
     
  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*);
    
private:

  NA49StackingMessenger*  stackMessenger;

  G4String                    pname;
  G4bool                      killSecondary;
  const G4ParticleDefinition* primaryDef;
  G4double                    primaryTotalEnergy;
  G4Element*                  elm;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

