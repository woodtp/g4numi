
#include "NumiDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"
#include "G4Cons.hh"
#include "G4Torus.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4PVReplica.hh"
#include "G4AssemblyVolume.hh"
#include "NumiMagneticField.hh"
#include "G4Transform3D.hh"

NumiDetectorConstruction::NumiDetectorConstruction()
{
    //Scan the input file     
    NumiData=new NumiDataInput();

// Pointers for magnetic fields ***    
    numiMagField = new NumiMagneticField(); 
    numiMagFieldIC = new NumiMagneticFieldIC();
    numiMagFieldOC = new NumiMagneticFieldOC();
}

NumiDetectorConstruction::~NumiDetectorConstruction()
{
  delete numiMagField; 
  delete numiMagFieldIC;
  delete numiMagFieldOC;
  delete NumiData;
  DestroyMaterials();
}

G4VPhysicalVolume* NumiDetectorConstruction::Construct()
{
  DefineMaterials();

  //Define world volume
  G4Tubs* ROCK_solid = new G4Tubs("ROCK_solid",0.,NumiData->RockRadius,NumiData->RockHalfLen,0,360.*deg);
  ROCK_log = new G4LogicalVolume(ROCK_solid,Concrete,"ROCK_log",0,0,0); 
  ROCK_log->SetVisAttributes(G4VisAttributes::Invisible);
  ROCK = new G4PVPlacement(0,G4ThreeVector(),ROCK_log,"ROCK",0,false,0);
 
  ConstructTargetHall();
  ConstructDecayPipe();
  ConstructTarget();
  ConstructHorns();
  ConstructHadronAbsorber();

  return ROCK;
}
