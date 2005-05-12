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
#include "G4AssemblyVolume.hh"
#include "NumiDataInput.hh"

void NumiDetectorConstruction::ConstructBaffle()

{  G4ThreeVector target_hall_position=G4ThreeVector(0,0,NumiData->TargetAreaLength/2.+NumiData->TargetAreaZ0);

 //HPBV
  G4double HPBRin=NumiData->HPBaffleRin;
  G4double HPBRout=NumiData->HPBaffleRout;
  G4double HPBlength=NumiData->HPBaffleLength;
  G4Tubs* SBaffle=new G4Tubs("SBaffle",HPBRin,HPBRout,HPBlength/2.,0.,360.*deg);
  G4LogicalVolume* LVBaffle=new G4LogicalVolume(SBaffle,Target,"LVBaffle",0,0,0);
  G4ThreeVector Baffle_pos=G4ThreeVector(NumiData->HPBaffleX0,NumiData->HPBaffleY0,NumiData->HPBaffleZ0+HPBlength/2.)-target_hall_position;
  new G4PVPlacement(0,Baffle_pos,"Baffle",LVBaffle,TGAR,false,0);
/*  
//HPBV 
  G4double HPBheight=NumiData->HPBaffleHeight;
  G4double HPBwidth=NumiData->HPBaffleWidth;
  G4double HPBlength=NumiData->HPBaffleLength;
  G4Box* HPBV_solid=new G4Box("HPBV_solid",HPBwidth/2.,HPBheight/2.,HPBlength/2.);
  G4LogicalVolume* HPBV_log=new G4LogicalVolume(HPBV_solid,Air,"HPBV_log",0,0,0);
  G4ThreeVector HPB_position=G4ThreeVector(NumiData->HPBaffleX0,NumiData->HPBaffleY0,NumiData->HPBaffleZ0+HPBlength/2.)-target_hall_position;
  G4VPhysicalVolume* HPBV = new G4PVPlacement(0,HPB_position,"HPBV",HPBV_log,TGAR,false,0);
  
  G4Box* HPBT_solid=new G4Box("HPBT_solid",5.*cm,2.2*cm,HPBlength/2.);
  G4LogicalVolume* HPBT_log=new G4LogicalVolume(HPBT_solid,Target,"HPBT_log",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,2.8*cm,0),"HPBT",HPBT_log,HPBV,false,0);
  
  G4Box* HPBB_solid=new G4Box("HPBB_solid",5.*cm,2.2*cm,HPBlength/2.);
  G4LogicalVolume* HPBB_log=new G4LogicalVolume(HPBB_solid,Target,"HPBB_log",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,-2.8*cm,0),"HPBB",HPBB_log,HPBV,false,0);

  G4Box* HPBL_solid=new G4Box("HPBL_solid",2.365*cm,0.59995002*cm,HPBlength/2.);
  G4LogicalVolume* HPBL_log=new G4LogicalVolume(HPBL_solid,Target,"HPBL_log",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(2.635*cm,0,0),"HPBL",HPBL_log,HPBV,false,0);

  G4Box* HPBR_solid=new G4Box("HPBR_solid",2.365*cm,0.59995002*cm,HPBlength/2.);
  G4LogicalVolume* HPBR_log=new G4LogicalVolume(HPBR_solid,Target,"HPBR_log",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(-2.635*cm,0,0),"HPBR",HPBR_log,HPBV,false,0);
*/
}
