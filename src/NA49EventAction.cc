
#include "NA49EventAction.hh"
#include "G4Event.hh"
#include "NA49EventActionMessenger.hh"

#include "G4UImanager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

NA49EventAction::NA49EventAction():
  printModulo(100),
  nSelected(0),
  drawFlag("all"),
  debugStarted(false)
{
  eventMessenger = new NA49EventActionMessenger(this);
  UI = G4UImanager::GetUIpointer();
  selectedEvents.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

NA49EventAction::~NA49EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void NA49EventAction::BeginOfEventAction(const G4Event* evt)
{
  // New event
  G4int nEvt = evt->GetEventID();
  if(nEvt%1000 == 0)G4cout<<"EventID " <<nEvt<<G4endl;
 if(nSelected>0) {
    for(G4int i=0; i<nSelected; i++) {
      if(nEvt == selectedEvents[i]) {
        UI->ApplyCommand("/random/saveThisEvent");
        UI->ApplyCommand("/tracking/verbose  2");
        debugStarted = true;
        break;
      }
    }
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void NA49EventAction::EndOfEventAction(const G4Event*)
{
  if(debugStarted) {
    UI->ApplyCommand("/tracking/verbose  0");
    debugStarted = false;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
