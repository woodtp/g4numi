#include "NumiPrimaryMessenger.hh"

#include "NumiPrimaryGeneratorAction.hh"
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

NumiPrimaryMessenger::NumiPrimaryMessenger(NumiPrimaryGeneratorAction* RA)
  :PrimaryAction (RA)
{
  particleTable = G4ParticleTable::GetParticleTable();

  BeamDir = new G4UIdirectory("/NuMI/Beam/");
  BeamDir->SetGuidance("Beam Paramaters a la Fluka");
  
  setCosX = new G4UIcmdWithADouble("/NuMI/Beam/cosx",this);
  setCosX->SetGuidance("x direction cosine of the beam");
  setCosX->SetParameterName("cosx",true);
  setCosX->SetDefaultValue (0.0);
  
  setCosY = new G4UIcmdWithADouble("/NuMI/Beam/cosy",this);
  setCosY->SetGuidance("y direction cosine of the beam");
  setCosY->SetParameterName("cosy",true);
  setCosY->SetDefaultValue (0.0);

  setZ = new G4UIcmdWithADoubleAndUnit("/NuMI/Beam/zpos",this);
  setZ->SetGuidance("Z position of the beam (cm by default)");
  setZ->SetParameterName("zpos",true,true);
  setZ->SetDefaultValue (-70.0);
  setZ->SetDefaultUnit("cm");
  setZ->SetUnitCandidates("micron mm cm m km");

  setMomentum = new G4UIcmdWithADoubleAndUnit("/NuMI/Beam/momentum",this);
  setMomentum->SetGuidance("Average beam momentum in GeV");
  setMomentum->SetParameterName("momentum",true);
  setMomentum->SetDefaultValue (10.);
  setMomentum->SetDefaultUnit("GeV");
  setMomentum->SetUnitCandidates("eV keV MeV GeV TeV");

  setSpread = new G4UIcmdWithADoubleAndUnit("/NuMI/Beam/spread",this);
  setSpread->SetGuidance("Momentum Spread in GeV");
  setSpread->SetParameterName("spread",true);
  setSpread->SetDefaultValue (0.);
  setSpread->SetRange("spread >= 0");


  setDivergence = new G4UIcmdWithADoubleAndUnit("/NuMI/Beam/divergence",this);
  setDivergence->SetGuidance("Beam Divergence in mrad");
  setDivergence->SetParameterName("divergence",true);
  setDivergence->SetDefaultValue (20.);


  setInnerR = new G4UIcmdWithADoubleAndUnit("/NuMI/Beam/innerR",this);
  setInnerR->SetGuidance("Inner Radius of annular beam spot (cm by default)");
  setInnerR->SetParameterName("innerR",true,true);
  setInnerR->SetDefaultValue (0.);
  setInnerR->SetDefaultUnit("cm");
  setInnerR->SetUnitCandidates("micron mm cm m km");

  setOuterR = new G4UIcmdWithADoubleAndUnit("/NuMI/Beam/outerR",this);
  setOuterR->SetGuidance("Outer Radius of annular beam spot (cm by default)");
  setOuterR->SetParameterName("outerR",true, true);
  setOuterR->SetDefaultValue (0.);
  setOuterR->SetDefaultUnit("cm");
  setOuterR->SetUnitCandidates("micron mm cm m km");

 
  setParticle = new G4UIcmdWithAString("/NuMI/Beam/particle",this);
  setParticle->SetGuidance("Name of beam particle");
  setParticle->SetParameterName("particle",true);
  setParticle->SetDefaultValue ("pi+");
  G4String candidateList; 
  G4int nPtcl = particleTable->entries();
  for(G4int i=0;i<nPtcl;i++)
  {
    candidateList += particleTable->GetParticleName(i);
    candidateList += " ";
  }
  setParticle->SetCandidates(candidateList);


  //PrimaryAction->PG_Z(-36079.0*cm);
  PrimaryAction->PG_Z(-70*cm);
  PrimaryAction->PG_Momentum(10.0*GeV);
  PrimaryAction->PG_Spread(0.0);
  PrimaryAction->PG_Divergence(20*mrad);
  PrimaryAction->PG_InnerR(0.0);
  PrimaryAction->PG_OuterR(0.0);
  PrimaryAction->PG_Particle(particleTable->FindParticle("pi+"));
  
}

NumiPrimaryMessenger::~NumiPrimaryMessenger()
{
  delete setCosX;
  delete setCosY;
  delete setZ;
  delete setMomentum;
  delete setSpread;
  delete setDivergence;
  delete setInnerR;
  delete setOuterR;
  delete setParticle;
}

void NumiPrimaryMessenger::SetNewValue(G4UIcommand* command,G4String newValues)
{
  G4cout << "NumiPrimaryMessenger: " << newValues << G4endl;
  if ( command==setCosX ) {
	PrimaryAction->PG_CosX(setCosX->GetNewDoubleValue(newValues));
  }
  else if ( command==setCosY ) {
	PrimaryAction->PG_CosY(setCosY->GetNewDoubleValue(newValues));
  }
  else if ( command==setZ ) {
	PrimaryAction->PG_Z(setZ->GetNewDoubleValue(newValues));
  }
  else if( command==setMomentum ) { 
	PrimaryAction->PG_Momentum(setMomentum->GetNewDoubleValue(newValues));
  }
  else if( command==setSpread ) { 
	PrimaryAction->PG_Spread(setSpread->GetNewDoubleValue(newValues));
  }
  else if( command==setDivergence ) { 
	PrimaryAction->PG_Divergence(setDivergence->GetNewDoubleValue(newValues));
  }
  else if( command==setInnerR ) { 
	PrimaryAction->PG_InnerR(setInnerR->GetNewDoubleValue(newValues));
  }
  else if( command==setOuterR ) { 
	PrimaryAction->PG_OuterR(setOuterR->GetNewDoubleValue(newValues));
  }
  else if( command==setParticle ) {
	G4ParticleDefinition* pd = particleTable->FindParticle(newValues);
    if(pd != NULL) {
	  PrimaryAction->PG_Particle(pd);
	}
  }
  G4cout << "NumiPrimaryMessenger complete" << G4endl;
}

