#include "NumiDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "NumiDataInput.hh" 

void NumiDetectorConstruction::ConstructTargetHall()
{

//TGAR
  G4double TGAR_width=NumiData->TargetAreaWidth/2.;
  G4double TGAR_height=NumiData->TargetAreaHeight/2.;
  G4double TGAR_length=NumiData->TargetAreaLength/2.;
  G4ThreeVector targetHallPosition=G4ThreeVector(0,0,NumiData->TargetAreaLength/2.+NumiData->TargetAreaZ0);

  G4Box* sTGAR=new G4Box("sTGAR",TGAR_width,TGAR_height,TGAR_length);
  G4LogicalVolume *lvTGAR=new G4LogicalVolume(sTGAR,Air,"lvTGAR",0,0,0);
  TGAR = new G4PVPlacement(0,targetHallPosition,"TGAR",lvTGAR,ROCK,false,0);
   
  // G4VisAttributes invisible=G4VisAttributes(false);
  //Blocks
  for (G4int ii=0;ii<NumiData->THBlockNblock;ii++){
    char no[3];
    sprintf(no,"%d",ii+1);
    G4String volName="sTHBLK";
    volName.append(no);
    BLK_solid[ii] = new G4Box(volName,NumiData->THBlockHdx[ii],NumiData->THBlockHdy[ii],NumiData->THBlockLength[ii]/2.);
    volName="lvTHBLK"; volName.append(no);
    G4Material* material=GetMaterial(NumiData->THBlockGeantMaterial[ii]);
    G4LogicalVolume *lvTHBLK = new G4LogicalVolume(BLK_solid[ii],material,volName,0,0,0);
    //      BLK_log[ii]->SetVisAttributes(invisible);
    volName="THBLK"; volName.append(no);
    G4ThreeVector blockPosition=G4ThreeVector(NumiData->THBlockX0[ii],NumiData->THBlockY0[ii],NumiData->THBlockZ0[ii]+NumiData->THBlockLength[ii]/2.)-targetHallPosition;
    new G4PVPlacement(0,blockPosition,volName,lvTHBLK,TGAR,false,0);
  }
  
  G4cout << "Target Hall Constructed" << G4endl;
}
 
