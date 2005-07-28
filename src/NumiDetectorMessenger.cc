
#include "NumiDetectorMessenger.hh"
#include "NumiDetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"


NumiDetectorMessenger::NumiDetectorMessenger( NumiDetectorConstruction* NumiDet):NumiDetector(NumiDet) {

	NumiDir = new G4UIdirectory("/NuMi/");
	NumiDir->SetGuidance("UI commands for detector geometry");

	detDir = new G4UIdirectory("/NuMi/det/");
	detDir->SetGuidance("detector control");

	TargetGasCmd = new G4UIcmdWithAString("/NuMi/det/setTarGas",this);
	TargetGasCmd->SetGuidance("Select gas inside the target.");
	TargetGasCmd->SetParameterName("choice",false);
	TargetGasCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
	TargetZ0Cmd = new G4UIcmdWithADoubleAndUnit("/NuMi/det/setTargetZ0",this);
	TargetZ0Cmd->SetGuidance("Set Z0 position of target");
	TargetZ0Cmd->SetParameterName("Size",false);
	TargetZ0Cmd->SetRange("Size<=0.");
	TargetZ0Cmd->SetUnitCategory("Length");
	TargetZ0Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	HornCurrentCmd = new G4UIcmdWithADoubleAndUnit("/NuMi/det/setHornCurrent",this);
	HornCurrentCmd->SetGuidance("Set horn current");
	HornCurrentCmd->SetParameterName("current",false);
	HornCurrentCmd->SetRange("current>=0.");
	HornCurrentCmd->SetUnitCategory("ampere");
	HornCurrentCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	UpdateCmd = new G4UIcmdWithoutParameter("/NuMi/det/update",this);
	UpdateCmd->SetGuidance("Update NuMi geometry.");
	UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
	UpdateCmd->SetGuidance("if you changed geometrical value(s).");
	UpdateCmd->AvailableForStates(G4State_Idle);

}

NumiDetectorMessenger::~NumiDetectorMessenger() {

	delete TargetGasCmd;
	delete TargetZ0Cmd;
	delete HornCurrentCmd;
	delete UpdateCmd;

	delete detDir;
	delete NumiDir;
}


void NumiDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue){
	if ( command == TargetGasCmd ) {
		//NumiDetector->SetTargetGas(newValue);
	}
	if ( command == TargetZ0Cmd ) {
		NumiDetector->SetTargetZ0(TargetZ0Cmd->GetNewDoubleValue(newValue));
	}
	if ( command == HornCurrentCmd ) {
		NumiDetector->SetHornCurrent(HornCurrentCmd->GetNewDoubleValue(newValue));
	}

	if ( command == UpdateCmd ) {
		NumiDetector->UpdateGeometry();
	}
	
}
	
	
