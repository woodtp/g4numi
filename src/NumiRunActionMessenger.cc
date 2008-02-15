//----------------------------------------------------------------------
// $Id
//----------------------------------------------------------------------

#include "NumiRunActionMessenger.hh"

#include "NumiRunAction.hh"
#include "NumiDataInput.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4RunManager.hh"

NumiRunActionMessenger::NumiRunActionMessenger(NumiRunAction* RA)
  :runAction (RA)
{
  NumiData=NumiDataInput::GetNumiDataInput();

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

  debugOn = new G4UIcmdWithABool("/NuMI/run/debugOn",this);
  debugOn->SetGuidance("Output some debugging info");
  debugOn->SetParameterName("debugOn",true);
  debugOn->SetDefaultValue (NumiData->IsDebugOn());
  debugOn->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  setRunID = new G4UIcmdWithAnInteger("/NuMI/run/setRunID",this);
  setRunID->SetGuidance("set run ID.");
  setRunID->SetParameterName("run ID number",true);
  setRunID->SetDefaultValue (0);
  setRunID->AvailableForStates(G4State_PreInit,G4State_Idle); 

  setMaterialSigma = new G4UIcmdWithADouble("/material/materialSigma",this);
  setMaterialSigma->SetGuidance("set Sigma");
  setMaterialSigma->SetParameterName("Sigma value of density/material",true);
  setMaterialSigma->SetDefaultValue (0);
  setMaterialSigma->AvailableForStates(G4State_PreInit,G4State_Idle); 

  useNImpWeight = new G4UIcmdWithABool("/NuMI/run/useNImpWeight",this);
  useNImpWeight->SetGuidance("use importance weighting (true/false)");
  useNImpWeight->SetParameterName("NImpWeight",true);
  useNImpWeight->SetDefaultValue (NumiData->NImpWeightOn);
  useNImpWeight->AvailableForStates(G4State_PreInit,G4State_Idle);

  useFlukaInput = new G4UIcmdWithABool("/NuMI/run/useFlukaInput",this);
  useFlukaInput->SetGuidance("use fluka input ntuple");
  useFlukaInput->SetParameterName("useFlukaInput",true);
  useFlukaInput->SetDefaultValue (NumiData->useFlukaInput);
  useFlukaInput->AvailableForStates(G4State_PreInit,G4State_Idle);

  useMarsInput = new G4UIcmdWithABool("/NuMI/run/useMarsInput",this);
  useMarsInput->SetGuidance("use mars input ntuple");
  useMarsInput->SetParameterName("useMarsInput",true);
  useMarsInput->SetDefaultValue (NumiData->useMarsInput);
  useMarsInput->AvailableForStates(G4State_PreInit,G4State_Idle);

  extNtupleFileName = new G4UIcmdWithAString("/NuMI/run/extNtupleFileName",this);
  extNtupleFileName->SetGuidance("set external (fluka/mars) ntuple file name");
  extNtupleFileName->SetParameterName("extNtupleFileName",true);
  extNtupleFileName->SetDefaultValue (NumiData->GetExtNtupleFileName());
  extNtupleFileName->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  NumiOutputDir = new G4UIdirectory("/NuMI/output/");
  NumiOutputDir->SetGuidance("NuMI output management");

  setNuNtupleFile = new G4UIcmdWithAString("/NuMI/output/setNuNtupleFile",this);
  setNuNtupleFile->SetGuidance("set neutrino ntuple file name");
  setNuNtupleFile->SetParameterName("fileName",true);
  setNuNtupleFile->SetDefaultValue (NumiData->nuNtupleName);
  setNuNtupleFile->AvailableForStates(G4State_PreInit,G4State_Idle);

  outputNuNtuple = new G4UIcmdWithABool("/NuMI/output/outputNuNtuple",this);
  outputNuNtuple->SetGuidance("output neutrino ntuple (true/false)");
  outputNuNtuple->SetParameterName("outputNuNtuple",true);
  outputNuNtuple->SetDefaultValue (NumiData->createNuNtuple);
  outputNuNtuple->AvailableForStates(G4State_PreInit,G4State_Idle);

  setHadmmNtupleFile = new G4UIcmdWithAString("/NuMI/output/setHadmmNtupleFile",this);
  setHadmmNtupleFile->SetGuidance("set hadron and muon monitor ntuple file name");
  setHadmmNtupleFile->SetParameterName("fileName",true);
  setHadmmNtupleFile->SetDefaultValue (NumiData->hadmmNtupleName);
  setHadmmNtupleFile->AvailableForStates(G4State_PreInit,G4State_Idle);

  outputHadmmNtuple = new G4UIcmdWithABool("/NuMI/output/outputHadmmNtuple",this);
  outputHadmmNtuple->SetGuidance("output hadron and muon monitor ntuple (true/false)");
  outputHadmmNtuple->SetParameterName("outputHadmmNtuple",true);
  outputHadmmNtuple->SetDefaultValue (NumiData->createHadmmNtuple);
  outputHadmmNtuple->AvailableForStates(G4State_PreInit,G4State_Idle);

  setASCIIFile = new G4UIcmdWithAString("/NuMI/output/setASCIIFile",this);
  setASCIIFile->SetGuidance("set ASCII file name");
  setASCIIFile->SetParameterName("fileName",true);
  setASCIIFile->SetDefaultValue (NumiData->asciiName);
  setASCIIFile->AvailableForStates(G4State_PreInit,G4State_Idle);

  outputASCIIFile = new G4UIcmdWithABool("/NuMI/output/outputASCIIFile",this);
  outputASCIIFile->SetGuidance("output ASCII file (true/false)");
  outputASCIIFile->SetParameterName("outputASCIIFile",true);
  outputASCIIFile->SetDefaultValue (NumiData->createASCII);
  outputASCIIFile->AvailableForStates(G4State_PreInit,G4State_Idle);
}

NumiRunActionMessenger::~NumiRunActionMessenger()
{
  delete readRndmCmd;
  delete showRndmCmd;
  delete setRndmSeedCmd;
  delete RndmDir;

  delete setRunID;
  delete setMaterialSigma;
  delete useNImpWeight;
  delete NumiRunDir;
  delete useFlukaInput;
  delete useMarsInput;

  delete setNuNtupleFile;
  delete outputNuNtuple;
  delete setHadmmNtupleFile;
  delete outputHadmmNtuple;
  delete setASCIIFile;
  delete outputASCIIFile;
  delete NumiOutputDir;
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

  if (command == setMaterialSigma){
    NumiData->materialSigma = setMaterialSigma->GetNewDoubleValue(newValues);
  }

  if (command == debugOn){
    NumiData->SetDebugOn(debugOn->GetNewBoolValue(newValues));
  }

  if (command == useNImpWeight){
    NumiData->SetNImpWeight(useNImpWeight->GetNewBoolValue(newValues));
  }

  if (command == useFlukaInput){
    NumiData->SetFlukaInput(useFlukaInput->GetNewBoolValue(newValues));
  }

  if (command == useMarsInput){
    NumiData->SetMarsInput(useMarsInput->GetNewBoolValue(newValues));
  }
  if (command == extNtupleFileName){
    NumiData->SetExtNtupleFileName(newValues);
  }
  if (command == setNuNtupleFile){
      NumiData->SetNuNtupleName(newValues);
  }

  if (command == outputNuNtuple){
      NumiData->OutputNuNtuple(outputNuNtuple->GetNewBoolValue(newValues));
  }

  if (command == setHadmmNtupleFile){
      NumiData->SetHadmmNtupleName(newValues);
  }

  if (command == outputHadmmNtuple){
      NumiData->OutputHadmmNtuple(outputNuNtuple->GetNewBoolValue(newValues));
  }

 if (command == setASCIIFile){
      NumiData->SetASCIIName(newValues);
  }

  if (command == outputASCIIFile){
      NumiData->OutputASCII(outputASCIIFile->GetNewBoolValue(newValues));
  }
}

