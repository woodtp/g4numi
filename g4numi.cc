
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "NumiDetectorConstruction.hh"
#include "QGSP.hh"
#include "NumiPrimaryGeneratorAction.hh"
#include "NumiStackingAction.hh"
#include "NumiSteppingAction.hh"
#include "NumiTrackingAction.hh"
#include "NumiRunAction.hh"

#ifdef G4VIS_USE
#include "NumiVisManager.hh"
#endif


int main(int argc,char** argv)
{
  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // set mandatory initialization classes
  runManager->SetUserInitialization(new NumiDetectorConstruction);
  QGSP * theQGSP = new QGSP;
  runManager->SetUserInitialization(theQGSP);

#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new NumiVisManager;
  visManager->Initialize();
#endif

  // set mandatory user action class
  runManager->SetUserAction(new NumiPrimaryGeneratorAction);

  // set user action classes
  runManager->SetUserAction(new NumiSteppingAction);
  runManager->SetUserAction(new NumiStackingAction);
  runManager->SetUserAction(new NumiTrackingAction);
  runManager->SetUserAction(new NumiRunAction);

  // Initialize G4 kernel
  runManager->Initialize();

  // get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();

  if(argc==1)
  // Define (G)UI terminal for interactive mode  
  { 
    G4UIsession * session = 0;

    // G4UIterminal is a (dumb) terminal.
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif    

    session->SessionStart();
    delete session;
  }
  else
  // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }

#ifdef G4VIS_USE
  delete visManager;
#endif

#  // job termination
  delete runManager;
  return 0;
}


