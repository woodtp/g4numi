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

  G4ThreeVector tunnel_position=G4ThreeVector(0,0,NumiData->TunnelLength/2.+NumiData->TunnelZ0);
  for (G4int ii=0;ii<NumiData->BlockNblock;ii++){
    if ((NumiData->BlockZ0[ii]>=NumiData->TunnelZ0) && 
	(NumiData->BlockZ0[ii]<=(NumiData->TunnelZ0 + NumiData->TunnelLength))){
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
      G4ThreeVector block_position=G4ThreeVector(NumiData->BlockX0[ii],NumiData->BlockY0[ii],NumiData->BlockZ0[ii]+NumiData->BlockLength[ii]/2.)-tunnel_position;
      new G4PVPlacement(0,block_position,vol_name,BLK_log[ii],TUNE,false,0);
    }
  }
  G4cout << "Hadron Absorber Constructed" << G4endl;
}


