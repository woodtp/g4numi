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
      G4Material* material=GetMaterial(NumiData->BlockGeantMaterial[ii]);
      BLK_log[ii] = new G4LogicalVolume(BLK_solid[ii],material,vol_name,0,0,0);
      vol_name="BLK";
      vol_name.append(no);
      G4ThreeVector block_position=G4ThreeVector(NumiData->BlockX0[ii],NumiData->BlockY0[ii],NumiData->BlockZ0[ii]+NumiData->BlockLength[ii]/2.)-target_hall_position;
      new G4PVPlacement(0,block_position,vol_name,BLK_log[ii],GetPhysicalVolume(NumiData->BlockMotherVolume[ii]),false,0);
    }
  }
  G4cout << "Blocks Constructed" << G4endl;
}
 
