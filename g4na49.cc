
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh"

#include "NA49DetectorConstruction.hh"
#include "NA49PhysicsList.hh"
//#include "QBBC.hh"
//#include "QGSP.hh"
#include "FTFP_BERT.hh"

#include "NA49PrimaryGeneratorAction.hh"

#include "NA49RunAction.hh"
#include "NA49EventAction.hh"
#include "NA49StackingAction.hh"
#include "NA49TrackingAction.hh"

#include "G4VisExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());

  //Construct the default run manager
  G4RunManager * runManager = new G4RunManager();

  //set mandatory initialization classes
  runManager->SetUserInitialization(new NA49DetectorConstruction());

  FTFP_BERT* physicsList = new FTFP_BERT;

  runManager->SetUserInitialization(physicsList);
  //runManager->SetUserInitialization(new QBBC(1,"QBEC_HP"));
  //runManager->SetUserInitialization(new QGSP);

  runManager->SetUserAction(new NA49PrimaryGeneratorAction());

  //set user action classes
  runManager->SetUserAction(new NA49RunAction());
  runManager->SetUserAction(new NA49EventAction());
  runManager->SetUserAction(new NA49StackingAction());
  runManager->SetUserAction(new NA49TrackingAction());

  //get the pointer to the User Interface manager
  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4VisManager* visManager = 0;

  if (argc==1)   // Define UI terminal for interactive mode
    {
#ifdef G4VIS_USE
      //visualization manager
      visManager = new G4VisExecutive;
      visManager->Initialize();
#endif
      G4UIsession* session = 0;
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);
#else
      session = new G4UIterminal();
#endif
      session->SessionStart();
      delete session;
    }
  else           // Batch mode
    {
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UI->ApplyCommand(command+fileName);
    }

  //job termination
  if(visManager) delete visManager;
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
