//
// NumiSteppingAction.hh
//

#ifndef NumiSteppingAction_H
#define NumiSteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"
#include "NumiDataInput.hh"

class G4EventManager;
class NumiEventAction;

class NumiSteppingAction : public G4UserSteppingAction
{
  
 public:
  NumiSteppingAction();
  virtual ~NumiSteppingAction();
  
  virtual void UserSteppingAction(const G4Step*);

private:
   
  NumiDataInput *NDI;
   G4EventManager *EvtManager;
   NumiEventAction *NumiEvtAct;

   G4bool fPrintAllSteps;
   G4bool fPrintSplitting;

   G4bool fPrintMuAlcove1;
   G4bool fPrintMuAlcove2;
   G4bool fPrintMuAlcove3;
   G4bool fPrintDeltaAlcove1;
   G4bool fPrintDeltaAlcove2;
   G4bool fPrintDeltaAlcove3;

   G4bool fPrintProcesses;
   G4bool fPrintTouchableHistory;
  

};

#endif

