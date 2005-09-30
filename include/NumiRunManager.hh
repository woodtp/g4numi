//
// NumiRunManager.hh
//

#ifndef NumiRunManager_h
#define NumiRunManager_h 1

#include "G4RunManager.hh"

class NumiPrimaryGeneratorAction;

class NumiRunManager : public G4RunManager
{
 public:
  NumiRunManager();
  virtual ~NumiRunManager();

  virtual void BeamOn(G4int n_event,const char* macroFile=0,G4int n_select=-1);
  inline G4int GetNumberOfEvents(){
    return numberOfEventToBeProcessed;
  }

  G4int nEvents;

 protected:
  NumiPrimaryGeneratorAction * primaryGeneratorAction;

};

#endif
