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
#include "NumiDataInput.hh"

void NumiDetectorConstruction::ConstructHadronAbsorber()
{

  G4ThreeVector tunnelPos=G4ThreeVector(0,0,NumiData->TunnelLength/2.+NumiData->TunnelZ0);
  for (G4int ii=0;ii<NumiData->HABlockNblock;ii++){
    //if ((NumiData->BlockZ0[ii]>=NumiData->TunnelZ0) && 
    //	(NumiData->BlockZ0[ii]<=(NumiData->TunnelZ0 + NumiData->TunnelLength))){
      char no[3];
      sprintf(no,"%d",ii+1);
      G4String volName="sBLK"; volName.append(no);
      BLK_solid[ii] = new G4Box(volName,NumiData->HABlockHdx[ii],NumiData->HABlockHdy[ii],NumiData->HABlockLength[ii]/2.);
      volName="lvBLK"; volName.append(no);
      G4Material* material=GetMaterial(NumiData->HABlockGeantMaterial[ii]);
      BLK_log[ii] = new G4LogicalVolume(BLK_solid[ii],material,volName,0,0,0);
      volName="BLK"; volName.append(no);
      G4ThreeVector blockPos=G4ThreeVector(NumiData->HABlockX0[ii],NumiData->HABlockY0[ii],NumiData->HABlockZ0[ii]+NumiData->HABlockLength[ii]/2.)-tunnelPos;
      new G4PVPlacement(0,blockPos,volName,BLK_log[ii],pvTUNE,false,0);
      //}
  }
  G4cout << "Hadron Absorber Constructed" << G4endl;

}


