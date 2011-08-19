
#ifndef NA49TargetSD_h
#define NA49TargetSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Step;
class G4TouchableHistory;
class G4HCofThisEvent;

class NA49TargetSD : public G4VSensitiveDetector
{
public: // Without description

  NA49TargetSD(const G4String&);
  virtual ~NA49TargetSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*,G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);
  void clear();
  void PrintAll();

private:

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

