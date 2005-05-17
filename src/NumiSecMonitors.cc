#include "NumiDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "NumiDataInput.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"


void NumiDetectorConstruction::ConstructSecMonitors()
{
 G4ThreeVector tunnel_position=G4ThreeVector(0,0,NumiData->TunnelLength/2.+NumiData->TunnelZ0);

  G4Box *Shadmon=new G4Box("SHadMon",2.*m,2.*m,0.5*mm);
  G4LogicalVolume *LVhadmon=new G4LogicalVolume(Shadmon,Vacuum,"LVHadMon",0,0,0);
  G4ThreeVector hadmon_pos=G4ThreeVector(0.,0.,723.25*m)-tunnel_position;
  
  new G4PVPlacement(0,hadmon_pos,"PVHadMon",LVhadmon,TUNE,false,0);
  
  
  G4cout << "Hadron Monitor Constructed" << G4endl;
  
  G4Box *Smumona1=new G4Box("SMuMonA1",2.*m,2.*m,0.5*mm);
  G4LogicalVolume *LVmumona1=new G4LogicalVolume(Smumona1,Vacuum,"LVMuMonA1",0,0,0);
  G4ThreeVector mumona1_pos=G4ThreeVector(0.,0.,739.*m)-tunnel_position;
  
  new G4PVPlacement(0,mumona1_pos,"PVMuMonA1",LVmumona1,TUNE,false,0);
  
  G4cout << "Muon Monitor Alcove 1 Constructed" << G4endl;
  
}
