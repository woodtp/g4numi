//
// NumiTrackingAction.hh
//
#ifndef NumiTrackingAction_h
#define NumiTrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "NumiTrajectory.hh"

class NumiTrackingAction : public G4UserTrackingAction
{
  public:
    NumiTrackingAction();
    virtual ~NumiTrackingAction();

    virtual void PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);
};

#endif
