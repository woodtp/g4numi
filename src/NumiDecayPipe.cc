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

void NumiDetectorConstruction::ConstructDecayPipe()
{
  // Tunnel
  G4double l=NumiData->TunnelLength/2;
  G4double r=NumiData->TunnelRadius;
  G4ThreeVector tunnel_position=G4ThreeVector(0,0,l+NumiData->TunnelZ0);
 
  G4Tubs* TUNE_solid = new G4Tubs("TUNE_solid",0.,r,l,0,360.*deg);
  TUNE_log = new G4LogicalVolume(TUNE_solid,GetMaterial(NumiData->TunnelGEANTmat),"TUNE_log",0,0,0); 
  TUNE_log->SetVisAttributes(G4VisAttributes::Invisible);
  TUNE = new G4PVPlacement(0,tunnel_position,"TUNE",TUNE_log,ROCK,false,0);

  // Shielding
  l=NumiData->ShieldLength/2.;
  G4double r_in=NumiData->ShieldRin;
  G4double r_out=NumiData->ShieldRout;
  G4ThreeVector sc01_position=G4ThreeVector(NumiData->ShieldX0,NumiData->ShieldY0,NumiData->ShieldZ0+l)-tunnel_position;

  G4Tubs* SC01_solid = new G4Tubs("SC01_solid",r_in,r_out,l,0,360.*deg);
  G4LogicalVolume* SC01_log = new G4LogicalVolume(SC01_solid,GetMaterial(NumiData->ShieldGEANTmat),"SC01_log",0,0,0); 
  new G4PVPlacement(0,sc01_position,"SC01",SC01_log,TUNE,false,0);

  // Decay Pipe
  l=NumiData->DecayPipeLength/2.;
  r=NumiData->DecayPipeRadius;
  G4ThreeVector decay_pipe_position=G4ThreeVector(0,0,NumiData->DecayPipeZ0+l)-tunnel_position;

  G4Tubs* DPIP_solid = new G4Tubs("DPIP_solid",0.,r,l,0,360.*deg);
  G4LogicalVolume* DPIP_log = new G4LogicalVolume(DPIP_solid,GetMaterial(NumiData->DecayPipeGEANTmat),"DPIP_log",0,0,0); 
  G4VPhysicalVolume* DPIP = new G4PVPlacement(0,decay_pipe_position,"DPIP",DPIP_log,TUNE,false,0);
 
  // Decay Volume
  l=(NumiData->DecayPipeLength-NumiData->DecayPipeFWinThick-NumiData->DecayPipeEWinThick)/2.;
  r=NumiData->DecayPipeRadius-NumiData->DecayPipeWallThick;
  G4ThreeVector decay_volume_position=G4ThreeVector(0,0,(NumiData->DecayPipeFWinThick-NumiData->DecayPipeEWinThick)/2.);

  G4Tubs* DVOL_solid = new G4Tubs("DVOL_solid",0.,r,l,0,360.*deg);
  G4LogicalVolume* DVOL_log = new G4LogicalVolume(DVOL_solid,Vacuum,"DVOL_log",0,0,0); 
  new G4PVPlacement(0,decay_volume_position,"DVOL",DVOL_log,DPIP,false,0);

  // Upstream window
  l=NumiData->DecayPipeFWinThick/2.;
  r=NumiData->DecayPipeRadius;
  G4ThreeVector upwn_position=G4ThreeVector(0,0,-NumiData->DecayPipeLength/2.+l);

  G4Tubs* UPWN_solid = new G4Tubs("UPWN_solid",0.,r,l,0,360.*deg);
  G4LogicalVolume* UPWN_log = new G4LogicalVolume(UPWN_solid,GetMaterial(NumiData->DecayPipeFWinmat),"UPWN_log",0,0,0); 
  new G4PVPlacement(0,upwn_position,"UPWN",UPWN_log,DPIP,false,0);

  // Downstream window
  l=NumiData->DecayPipeEWinThick/2.;
  r=NumiData->DecayPipeRadius;
  G4ThreeVector dnwn_position=G4ThreeVector(0,0,NumiData->DecayPipeLength/2.-l);

  G4Tubs* DNWN_solid = new G4Tubs("DNWN_solid",0.,r,l,0,360.*deg);
  G4LogicalVolume* DNWN_log = new G4LogicalVolume(DNWN_solid,GetMaterial(NumiData->DecayPipeEWinmat),"DNWN_log",0,0,0); 
  new G4PVPlacement(0,dnwn_position,"DNWN",DNWN_log,DPIP,false,0);

  G4cout << "Decay Pipe Constructed" <<G4endl;
}
