#include "NumiDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "NumiDataInput.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include <math.h>

static const G4float in = 2.54*cm;
static const G4float ft = 12*in;
void NumiDetectorConstruction::ConstructSecMonitors()
{
  // Only want to retrieve some information once, i.e. z position of the end of the decay pipe.
  G4ThreeVector tunnelPos = G4ThreeVector(0,0,NumiData->TunnelLength/2.+NumiData->TunnelZ0);
  G4float DP_end = NumiData->DecayPipeLength + NumiData->TunnelZ0;
  G4float  MuAlcv_x = 0;
  G4float  MuAlcv_y = 0;
  G4float  MuAlcv_z = 0;

  G4float  MuAlcv_width = 12*ft;
  G4float  MuAlcv_height =12*ft;
  G4float  MuAlcv_length = 9*ft+10*in;
  // the last bit is the distance from the end of the decay pipe to the 
  // beginning of the Hadron Absorber
  G4float D_MA1 = NumiData->HadrBox_length+(9*ft+9*in)/2.0+(4*ft+2.25*in);
  G4float E_MA1 = 11.85*ft-(14*ft+8*in)/2.0;
  G4float H_MA1 = sqrt(pow(E_MA1,2)+pow(D_MA1,2));
  G4float a_MA1 = atan(E_MA1/D_MA1);
  G4float b_MA1 = 3.3*pi/180-a_MA1;
  G4float z_MA1 = H_MA1*cos(b_MA1)+DP_end;
  G4float y_MA1 = H_MA1*sin(b_MA1);

  G4float D_MA2 = D_MA1+(9*ft+9*in)/2.0+39*ft+6*in+(9*ft+10*in)/2.0;
  G4float E_MA2 = 11.85*ft-MuAlcv_height/2.0;
  G4float H_MA2 = sqrt(pow(E_MA2,2)+pow(D_MA2,2));
  G4float a_MA2 = atan(E_MA2/D_MA2);
  G4float b_MA2 = 3.3*pi/180-a_MA2;
  G4float z_MA2 = H_MA2*cos(b_MA2)+DP_end;
  G4float y_MA2 = H_MA2*sin(b_MA2);

  G4float D_MA3 = D_MA2+59*ft+3*in+9*ft+10*in;
  G4float E_MA3 = 11.85*ft+(464.50-460.48)*ft-MuAlcv_height/2.0;
  G4float H_MA3 = sqrt(pow(E_MA3,2)+pow(D_MA3,2));
  G4float a_MA3 = atan(E_MA3/D_MA3);
  G4float b_MA3 = 3.3*pi/180-a_MA3;
  G4float z_MA3 = H_MA3*cos(b_MA3)+DP_end;
  G4float y_MA3 = H_MA3*sin(b_MA3);

  //rotation to match vertical position of detector equipment to downslope angle of tunnel
  G4RotationMatrix *zrot = new G4RotationMatrix();
  zrot->rotateX(3.3*M_PI/180);

  // Gotta put the Hadron Monitor inside the Concrete Shielding G4Box
  // otherwise the particle doesn't exit to the pvTUNE volume and
  // tunnels past the monitor.

  G4Box *Shadmon = new G4Box("SHadMon",1.1*m,1.1*m,5*mm);
  G4LogicalVolume *LVhadmon = new G4LogicalVolume(Shadmon,Air,"LVHadMon",0,0,0);
  G4ThreeVector hadmon_pos = G4ThreeVector(0.,0.2*m,0.2*m);
  
  new G4PVPlacement(0,hadmon_pos,"PVHadMon",LVhadmon,ShldBox,false,0);
  
  G4cout << "Hadron Monitor Constructed" << G4endl;
  G4VPhysicalVolume* MuMonAlcv_0;
  G4VPhysicalVolume* MuMonAlcv_1;
  G4VPhysicalVolume* MuMonAlcv_2;
  
  // Muon monitor 1 has slightly different dimensions than the other three alcoves.
  MuAlcv_z = z_MA1;
  MuAlcv_y = y_MA1;
  G4Box *Smumonalcv_0 = new G4Box("SMuMonAlcv_0",(12*ft)/2.0,(14*ft+8*in)/2.0,(9*ft+9*in)/2.0);
  G4LogicalVolume *LVmumonalcv_0 = new G4LogicalVolume(Smumonalcv_0,Air,"LVMuMonAlcv_0",0,0,0);
  G4ThreeVector mumonalcv_0_pos = G4ThreeVector(0.,MuAlcv_y,MuAlcv_z);
  MuMonAlcv_0 = new G4PVPlacement(zrot,mumonalcv_0_pos,"MuMonAlcv_0",LVmumonalcv_0,ROCK,false,0);

  MuAlcv_z = z_MA2;
  MuAlcv_y = y_MA2;
  G4Box *Smumonalcv = new G4Box("SMuMonAlcv",MuAlcv_width/2.0,MuAlcv_height/2.0,MuAlcv_length/2.0);
  G4LogicalVolume *LVmumonalcv_1 = new G4LogicalVolume(Smumonalcv,Air,"LVMuMonAlcv_1",0,0,0);
  G4ThreeVector mumonalcv_pos = G4ThreeVector(0.,MuAlcv_y,MuAlcv_z);
  MuMonAlcv_1 = new G4PVPlacement(zrot,mumonalcv_pos,"MuMonAlcv_1",LVmumonalcv_1,ROCK,false,0);

  MuAlcv_z = z_MA3;
  MuAlcv_y = y_MA3;
  mumonalcv_pos = G4ThreeVector(0.,MuAlcv_y,MuAlcv_z);
  G4LogicalVolume *LVmumonalcv_2 = new G4LogicalVolume(Smumonalcv,Air,"LVMuMonAlcv_2",0,0,0);
  MuMonAlcv_2 = new G4PVPlacement(zrot,mumonalcv_pos,"MuMonAlcv_2",LVmumonalcv_2,ROCK,false,0);

  G4cout << "Muon Monitor Alcoves(1,2,3) Constructed" << G4endl;

  // The first iteration of the Muon Monitors will be 2m x 2m x 2.25in boxes which will become 
  // perfectly active elements via the NumiSteppingAction.cc tracking a particle (hopefully muon) 
  // through the monitor. 
  G4float MuMon_x = 2*m;
  G4float MuMon_y = 2*m;
  G4float MuMon_z = 2.25*in;

  G4Box *SMuMon = new G4Box("SMuMon",MuMon_x/2.0,MuMon_y/2.0,MuMon_z/2.0);
  G4LogicalVolume *LVMuMon = new G4LogicalVolume(SMuMon,Air,"LVMuMon",0,0,0);

  G4ThreeVector MuMon_pos0 = G4ThreeVector(0.,-(14*ft+8*in)/2.0+105*in,0.);
  G4ThreeVector MuMon_pos1_2 = G4ThreeVector(0.,0.,0.);
  new G4PVPlacement(0,MuMon_pos0,"MuMon_0",LVMuMon,MuMonAlcv_0,false,0);
  new G4PVPlacement(0,MuMon_pos1_2,"MuMon_1",LVMuMon,MuMonAlcv_1,false,0);
  new G4PVPlacement(0,MuMon_pos1_2,"MuMon_2",LVMuMon,MuMonAlcv_2,false,0);
  
}
