#include "NumiRunActionMessenger.hh"

#include "NumiRunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4RunManager.hh"


NumiRunActionMessenger::NumiRunActionMessenger(NumiRunAction* RA)
  :runAction (RA)
{
  RndmDir = new G4UIdirectory("/NuMI/rndm/");
  RndmDir->SetGuidance("Rndm status control.");
  
  readRndmCmd = new G4UIcmdWithAString("/NuMI/rndm/read",this);
  readRndmCmd->SetGuidance("get rndm status from an external file.");
  readRndmCmd->SetParameterName("fileName",true);
  readRndmCmd->SetDefaultValue ("");
  readRndmCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  showRndmCmd = new G4UIcmdWithoutParameter("/NuMI/rndm/show",this);
  showRndmCmd->SetGuidance("show rndm status.");
  showRndmCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  setRndmSeedCmd = new G4UIcmdWithAnInteger("/NuMI/rndm/setRndmSeed",this);
  setRndmSeedCmd->SetGuidance("set rndm seed.");
  setRndmSeedCmd->SetParameterName("rndmSeed",true);
  setRndmSeedCmd->SetDefaultValue (0);
  setRndmSeedCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NumiRunDir = new G4UIdirectory("/NuMI/run/");
  NumiRunDir->SetGuidance("NuMI run management");

  setRunID = new G4UIcmdWithAnInteger("/NuMI/run/setRunID",this);
  setRunID->SetGuidance("set run ID.");
  setRunID->SetParameterName("run ID number",true);
  setRunID->SetDefaultValue (0);
  setRunID->AvailableForStates(G4State_PreInit,G4State_Idle);  
}

NumiRunActionMessenger::~NumiRunActionMessenger()
{
  delete readRndmCmd;
  delete showRndmCmd;
  delete setRndmSeedCmd;
  delete RndmDir;

  delete setRunID;
  delete NumiRunDir;

}

void NumiRunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValues)
{
  if (command == readRndmCmd)
    { 
      G4cout << "\n---> rndm status restored from file: " << newValues << G4endl;
      G4String rndmFile="rndm/";
      rndmFile.append(newValues);
      HepRandom::restoreEngineStatus(rndmFile);
    }   

  if (command == showRndmCmd)
    { 
      HepRandom::showEngineStatus();
    }  

  if (command == setRndmSeedCmd)
    { 
      HepRandom::setTheSeed(setRndmSeedCmd->GetNewIntValue(newValues));
    }  

  if (command == setRunID){
    G4RunManager* runManager = G4RunManager::GetRunManager();
    runManager->SetRunIDCounter(setRunID->GetNewIntValue(newValues));
  }
}

