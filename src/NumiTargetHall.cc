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

void NumiDetectorConstruction::ConstructTargetHall()
{
  //TGAR
  G4double TGAR_width=NumiData->TargetAreaWidth/2.;
  G4double TGAR_height=NumiData->TargetAreaHeight/2.;
  G4double TGAR_length=NumiData->TargetAreaLength/2.;
  G4ThreeVector target_hall_position=G4ThreeVector(0,0,NumiData->TargetAreaLength/2.+NumiData->TargetAreaZ0);

  G4Box* TGAR_solid=new G4Box("TGAR_solid",TGAR_width,TGAR_height,TGAR_length);
  TGAR_log=new G4LogicalVolume(TGAR_solid,Air,"TGAR_log",0,0,0);
  TGAR = new G4PVPlacement(0,target_hall_position,"TGAR",TGAR_log,ROCK,false,0);
  
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

  //Blocks
  for (G4int ii=0;ii<NumiData->BlockNblock;ii++){
    if ((NumiData->BlockZ0[ii]>=NumiData->TargetAreaZ0) && 
	(NumiData->BlockZ0[ii]<=(NumiData->TargetAreaZ0 + NumiData->TargetAreaLength))){

      char no[3];
      sprintf(no,"%d",ii+1);
      G4String vol_name="BLK_solid";
      vol_name.append(no);
      BLK_solid[ii] = new G4Box(vol_name,NumiData->BlockHdx[ii],NumiData->BlockHdy[ii],NumiData->BlockLength[ii]/2.);
      vol_name="BLK_log";
      vol_name.append(no);
      G4Material* material=GetMaterial(NumiData->BlockGeantmat[ii]);
      BLK_log[ii] = new G4LogicalVolume(BLK_solid[ii],material,vol_name,0,0,0);
      vol_name="BLK";
      vol_name.append(no);
      G4ThreeVector block_position=G4ThreeVector(NumiData->BlockX0[ii],NumiData->BlockY0[ii],NumiData->BlockZ0[ii]+NumiData->BlockLength[ii]/2.)-target_hall_position;
      new G4PVPlacement(0,block_position,vol_name,BLK_log[ii],TGAR,false,0);
    }
  }
  G4cout << "Target Hall Constructed" << G4endl;
 }

