//
// NumiSteppingAction.hh
//

#ifndef NumiSteppingAction_H
#define NumiSteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class NumiSteppingAction : public G4UserSteppingAction
{

  public:
    NumiSteppingAction();
    virtual ~NumiSteppingAction();

    virtual void UserSteppingAction(const G4Step*);

};

#endif

