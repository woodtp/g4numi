//----------------------------------------------------------------------
// $Id
//----------------------------------------------------------------------

#include "NumiDetectorMessenger.hh"
#include "NumiDetectorConstruction.hh"
#include "NumiDataInput.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UnitsTable.hh"

NumiDetectorMessenger::NumiDetectorMessenger( NumiDetectorConstruction* NumiDet):NumiDetector(NumiDet)
{

        NumiDataInput *ND = NumiDataInput::GetNumiDataInput();

        if(ND->fPrintInfo > 0 || ND->IsDebugOn())
        {
           G4cout << "NumiDetectorMessenger Constructor Called." << G4endl;
        }

	NumiDir = new G4UIdirectory("/NuMI/");
	NumiDir->SetGuidance("UI commands for detector geometry");

	detDir = new G4UIdirectory("/NuMI/det/");
	detDir->SetGuidance("detector control");

	TargetGasCmd = new G4UIcmdWithAString("/NuMI/det/setTarGas",this);
	TargetGasCmd->SetGuidance("Select gas inside the target.");
	TargetGasCmd->SetParameterName("choice",false);
	TargetGasCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	ConstructTarget = new G4UIcmdWithABool("/NuMI/det/constructTarget",this); 
	ConstructTarget->SetGuidance("Target construction on/off"); 
	ConstructTarget->SetParameterName("constructTarget",true); 
	ConstructTarget->SetDefaultValue (ND->constructTarget); 
	ConstructTarget->AvailableForStates(G4State_PreInit,G4State_Idle);

	HeInDecayPipe = new G4UIcmdWithABool("/NuMI/det/heInDecayPipe",this); 
	HeInDecayPipe->SetGuidance("Insert 0.9atm, 300K He in the decay pipe - on/off"); 
	HeInDecayPipe->SetParameterName("heInDecayPipe",true); 
	HeInDecayPipe->SetDefaultValue (ND->HeInDecayPipe); 
	HeInDecayPipe->AvailableForStates(G4State_PreInit,G4State_Idle);

        
        LengthOfWaterInTgt = new G4UIcmdWithADoubleAndUnit("/NuMI/det/LengthOfWaterInTgt",this);
        LengthOfWaterInTgt->SetGuidance("Set length of water in target for Water In the Target simulation.");
	LengthOfWaterInTgt->SetParameterName("LengthOfWaterInTgt",true);
	LengthOfWaterInTgt->SetRange("LengthOfWaterInTgt>=3.");
	LengthOfWaterInTgt->SetUnitCategory("Length");
        LengthOfWaterInTgt->SetDefaultValue(ND->GetLengthOfWaterInTgt()); 
	LengthOfWaterInTgt->AvailableForStates(G4State_PreInit,G4State_Idle);
  
	//TargetZ0Cmd = new G4UIcmdWithADoubleAndUnit("/NuMI/det/setTargetZ0",this);
	//TargetZ0Cmd->SetGuidance("Set Z0 position of target");
	//TargetZ0Cmd->SetParameterName("Size",false);
	//TargetZ0Cmd->SetRange("Size<=0.");
	//TargetZ0Cmd->SetUnitCategory("Length");
	//TargetZ0Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	new G4UnitDefinition("kiloampere" , "kA", "Electric current", 1000.*ampere);
        //HornCurrentCmd = new G4UIcmdWithADoubleAndUnit("/NuMI/det/setHornCurrent",this);
	//HornCurrentCmd->SetGuidance("Set horn current");
	//HornCurrentCmd->SetParameterName("current",false);
	//HornCurrentCmd->SetDefaultValue(ND->HornCurrent);
	//HornCurrentCmd->SetRange("current>=0.");
	//HornCurrentCmd->SetUnitCategory("Electric current");
	//HornCurrentCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

        ConstructSolidMuMons = new G4UIcmdWithABool("/NuMI/det/constructSolidMuMons",this); 
	ConstructSolidMuMons->SetGuidance("Construct Full Solid Muon Monitors yes/no"); 
	ConstructSolidMuMons->SetParameterName("solidMuMons",true); 
	ConstructSolidMuMons->SetDefaultValue (ND->GetSolidMuMons()); 
	ConstructSolidMuMons->AvailableForStates(G4State_PreInit,G4State_Idle);


        RunPeriod = new G4UIcmdWithAnInteger("/NuMI/det/RunPeriod",this);
        RunPeriod -> SetGuidance("Set the Run Period, default is -999 which is no run period at all. MUST SET IT. Run Period 0 is like a \"default\" run period, everything is aligned, it is not specific to any ACTUAL NuMI beam running period.");
        RunPeriod -> SetParameterName("RunPeriod",false);
        RunPeriod -> SetDefaultValue (ND->GetRunPeriod()); 
        RunPeriod -> AvailableForStates(G4State_PreInit,G4State_Idle);
        
        BeamConfig = new G4UIcmdWithAString("/NuMI/det/BeamConfig",this);
        BeamConfig -> SetGuidance("Set the Beam Configuration. Must be in the form \"LE#z#i\" or \"le#z#i\". Note cannot handle ME config right now.");
        BeamConfig -> SetParameterName("BeamConfig",false);
        BeamConfig -> SetDefaultValue (ND->GetBeamConfig()); 
        BeamConfig -> AvailableForStates(G4State_PreInit,G4State_Idle);

        UseCorrHornCurrent = new G4UIcmdWithABool("/NuMI/det/UseCorrHornCurrent",this);
        UseCorrHornCurrent -> SetGuidance("Use the Calibration Corrected Horn Current value of the horn current configuration defined by the BeamConfiguration, true/false.");
        UseCorrHornCurrent -> SetParameterName("UseCorrHornCurrent",false);
        UseCorrHornCurrent -> AvailableForStates(G4State_PreInit,G4State_Idle);
        

        AbsorberConfig = new G4UIcmdWithAString("/NuMI/det/absorberConfig",this);
        AbsorberConfig -> SetGuidance("Set Absorber Configuration. Full Converage, Center Coverage, None, Configured.");
        AbsorberConfig -> SetParameterName("AbsorberConfig",false);
        AbsorberConfig -> SetDefaultValue (ND->GetAbsorberConfig()); 
        AbsorberConfig -> AvailableForStates(G4State_PreInit,G4State_Idle);

        Mon0AbsMatCmd = new G4UIcmdWithAString("/NuMI/det/setMon0AbsMater",this);
        Mon0AbsMatCmd->SetGuidance("Select Material of the Absorber in front of Mon 0.");
        Mon0AbsMatCmd->SetParameterName("choice",false);
        Mon0AbsMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

        Mon1AbsMatCmd = new G4UIcmdWithAString("/NuMI/det/setMon1AbsMater",this);
        Mon1AbsMatCmd->SetGuidance("Select Material of the Absorber in front of Mon 1.");
        Mon1AbsMatCmd->SetParameterName("choice",false);
        Mon1AbsMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

        Mon2AbsMatCmd = new G4UIcmdWithAString("/NuMI/det/setMon2AbsMater",this);
        Mon2AbsMatCmd->SetGuidance("Select Material of the Absorber in front of Mon 2.");
        Mon2AbsMatCmd->SetParameterName("choice",false);
        Mon2AbsMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

        Mon0AbsThickCmd = new G4UIcmdWithADoubleAndUnit("/NuMI/det/setMon0AbsThick",this);
        Mon0AbsThickCmd->SetGuidance("Set Thickness of the Absorber in front of Mon 0.");
        Mon0AbsThickCmd->SetParameterName("Abs0Thick",false);
        Mon0AbsThickCmd->SetRange("Abs0Thick>=0.");
        Mon0AbsThickCmd->SetUnitCategory("Length");
        Mon0AbsThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

        Mon1AbsThickCmd = new G4UIcmdWithADoubleAndUnit("/NuMI/det/setMon1AbsThick",this);
        Mon1AbsThickCmd->SetGuidance("Set Thickness of the Absorber in front of Mon 1.");
        Mon1AbsThickCmd->SetParameterName("Abs1Thick",false);
        Mon1AbsThickCmd->SetRange("Abs1Thick>=0.");
        Mon1AbsThickCmd->SetUnitCategory("Length");
        Mon1AbsThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

        Mon2AbsThickCmd = new G4UIcmdWithADoubleAndUnit("/NuMI/det/setMon2AbsThick",this);
        Mon2AbsThickCmd->SetGuidance("Set Thickness of the Absorber in front of Mon 2.");
        Mon2AbsThickCmd->SetParameterName("Abs2Thick",false);
        Mon2AbsThickCmd->SetRange("Abs2Thick>=0.");
        Mon2AbsThickCmd->SetUnitCategory("Length");
        Mon2AbsThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

        Mon0AbsDistCmd = new G4UIcmdWithADoubleAndUnit("/NuMI/det/setMon0AbsDist",this);
        Mon0AbsDistCmd->SetGuidance("Set Distance between the Absorber and Mon 0.");
        Mon0AbsDistCmd->SetParameterName("Abs0Dist",false);
        Mon0AbsDistCmd->SetRange("Abs0Dist>=0.");
        Mon0AbsDistCmd->SetUnitCategory("Length");
        Mon0AbsDistCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

        Mon1AbsDistCmd = new G4UIcmdWithADoubleAndUnit("/NuMI/det/setMon1AbsDist",this);
        Mon1AbsDistCmd->SetGuidance("Set Distance between the Absorber and Mon 1.");
        Mon1AbsDistCmd->SetParameterName("Abs1Dist",false);
        Mon1AbsDistCmd->SetRange("Abs1Dist>=0.");
        Mon1AbsDistCmd->SetUnitCategory("Length");
        Mon1AbsDistCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

        Mon2AbsDistCmd = new G4UIcmdWithADoubleAndUnit("/NuMI/det/setMon2AbsDist",this);
        Mon2AbsDistCmd->SetGuidance("Set Distance between the Absorber and Mon 2.");
        Mon2AbsDistCmd->SetParameterName("Abs2Dist",false);
        Mon2AbsDistCmd->SetRange("Abs2Dist>=0.");
        Mon2AbsDistCmd->SetUnitCategory("Length");
        Mon2AbsDistCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	UpdateCmd = new G4UIcmdWithoutParameter("/NuMI/det/update",this);
	UpdateCmd->SetGuidance("Update NuMi geometry.");
	UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
	UpdateCmd->SetGuidance("if you changed geometrical value(s).");
	UpdateCmd->AvailableForStates(G4State_Idle);

}

NumiDetectorMessenger::~NumiDetectorMessenger() {

	delete TargetGasCmd;
	//delete TargetZ0Cmd;
	//delete HornCurrentCmd;
        delete ConstructSolidMuMons;
        delete RunPeriod;
        delete BeamConfig;        
        delete AbsorberConfig;
        delete Mon0AbsMatCmd;
        delete Mon1AbsMatCmd;
        delete Mon2AbsMatCmd;
        delete Mon0AbsThickCmd;
        delete Mon1AbsThickCmd;
        delete Mon2AbsThickCmd;
        delete Mon0AbsDistCmd;
        delete Mon1AbsDistCmd;
        delete Mon2AbsDistCmd;
        delete UpdateCmd;
        delete LengthOfWaterInTgt;

	delete detDir;
	delete NumiDir;
}


void NumiDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{

   NumiDataInput *NumiData = NumiDataInput::GetNumiDataInput();

   if(NumiData->fPrintInfo > 0 || NumiData->IsDebugOn())
   {
      G4cout << "NumiDetectorMessenger::SetNewValue - Setting Parameter value from input macro." << G4endl;
   }

   
   if( command == RunPeriod )          { NumiData->SetRunPeriod(RunPeriod->GetNewIntValue(newValue)); }
   if( command == BeamConfig )         { NumiData->SetBeamConfig(newValue); }
   if( command == UseCorrHornCurrent)  { NumiData->SetUseCorrHornCurrent(UseCorrHornCurrent->GetNewBoolValue(newValue)); }
   if( command == LengthOfWaterInTgt)  { NumiData->SetLengthOfWaterInTgt(LengthOfWaterInTgt->GetNewDoubleValue(newValue));}

   
   if ( command == TargetGasCmd ) {
      //NumiDetector->SetTargetGas(newValue);
   }
   //if ( command == TargetZ0Cmd ) {
   //   NumiDetector->SetTargetZ0(TargetZ0Cmd->GetNewDoubleValue(newValue));
   //}
   //if ( command == HornCurrentCmd ) {
   //   NumiDetector->SetHornCurrent(HornCurrentCmd->GetNewDoubleValue(newValue));
   //}
   if ( command == ConstructTarget ) {
      NumiDataInput *NumiData=NumiDataInput::GetNumiDataInput();
      NumiData->constructTarget=ConstructTarget->GetNewBoolValue(newValue);
   }
   if ( command == HeInDecayPipe ) {
      NumiDataInput *NumiData=NumiDataInput::GetNumiDataInput();
      NumiData->HeInDecayPipe=HeInDecayPipe->GetNewBoolValue(newValue);
   }
   if ( command == ConstructSolidMuMons ) {
      NumiDataInput *NumiData=NumiDataInput::GetNumiDataInput();
      NumiData->SetSolidMuMons(ConstructSolidMuMons->GetNewBoolValue(newValue));
   }
   if( command == AbsorberConfig )
   {
      NumiDetector->SetAbsorberConfig(newValue);
   }
   
   if( command == Mon0AbsMatCmd )
   {
      NumiDetector->SetMonAbsorberMaterial(newValue, 0);
   }
   if( command == Mon1AbsMatCmd )
   {
           NumiDetector->SetMonAbsorberMaterial(newValue, 1);
   }
   if( command == Mon2AbsMatCmd )
   {
      NumiDetector->SetMonAbsorberMaterial(newValue, 2);
   }
   if( command == Mon0AbsThickCmd )
   {
      NumiDetector->SetMonAbsorberThickness(Mon0AbsThickCmd->GetNewDoubleValue(newValue), 0);
   }
   if( command == Mon1AbsThickCmd )
   {
      NumiDetector->SetMonAbsorberThickness(Mon1AbsThickCmd->GetNewDoubleValue(newValue), 1);
   }
   if( command == Mon2AbsThickCmd )
   {
      NumiDetector->SetMonAbsorberThickness(Mon2AbsThickCmd->GetNewDoubleValue(newValue), 2);
   }
   if( command == Mon0AbsDistCmd )
   {
      NumiDetector->SetAbsorberDistFromMon(Mon0AbsDistCmd->GetNewDoubleValue(newValue), 0);
   }
   if( command == Mon1AbsDistCmd )
   {
      NumiDetector->SetAbsorberDistFromMon(Mon1AbsDistCmd->GetNewDoubleValue(newValue), 1);
   }
   if( command == Mon2AbsDistCmd )
   {
      NumiDetector->SetAbsorberDistFromMon(Mon2AbsDistCmd->GetNewDoubleValue(newValue), 2);
   }
   
   if ( command == UpdateCmd ) {
#ifndef FLUGG
      NumiDetector->UpdateGeometry();
#endif
      return;
   }
   
}
	
	
