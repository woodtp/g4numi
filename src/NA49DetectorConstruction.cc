
#include "NA49DetectorConstruction.hh"
#include "NA49DetectorMessenger.hh"

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4RunManager.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UnitsTable.hh"
#include "G4ios.hh"

#include "NA49TargetSD.hh"
#include "NA49CheckVolumeSD.hh"
#include "G4SDManager.hh"
#include "G4NistManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NA49DetectorConstruction::NA49DetectorConstruction()
{
  logicTarget = 0;
  logicCheck  = 0;
  logicWorld  = 0;
  detectorMessenger = new NA49DetectorMessenger(this);

  radius = 0.3*cm;

  targetMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_GRAPHITE");
  worldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");

  // Prepare sensitive detectors
  checkSD = new NA49CheckVolumeSD("checkSD");
  (G4SDManager::GetSDMpointer())->AddNewDetector( checkSD );
  targetSD = new NA49TargetSD("targetSD");
  (G4SDManager::GetSDMpointer())->AddNewDetector( targetSD );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NA49DetectorConstruction::~NA49DetectorConstruction()
{ 
  delete detectorMessenger;
}

G4VPhysicalVolume* NA49DetectorConstruction::Construct()
{
  // Cleanup old geometry

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // Sizes
  G4double checkR  = radius + mm;
  G4double worldR  = radius + cm;
  G4double targetZ = 0.7*cm*0.5 ; 
  G4double checkZ  = targetZ + mm;
  G4double worldZ  = targetZ + cm;

  G4int nSlices    = 30;
  G4double sliceZ  = targetZ/G4double(nSlices);

  //
  // World
  //
  G4Tubs* solidW = new G4Tubs("World",0.,worldR,worldZ,0.,twopi);
  logicWorld = new G4LogicalVolume( solidW,worldMaterial,"World");
  G4VPhysicalVolume* world = new G4PVPlacement(0,G4ThreeVector(),
                                       logicWorld,"World",0,false,0);
  //
  // Check volume
  //
  G4Tubs* solidC = new G4Tubs("Check",0.,checkR,checkZ,0.,twopi);
  logicCheck = new G4LogicalVolume( solidC,worldMaterial,"World");
  //  G4VPhysicalVolume* physC = 
  new G4PVPlacement(0,G4ThreeVector(),logicCheck,"World",logicWorld,false,0);
  logicCheck->SetSensitiveDetector(checkSD);

  //
  // Target volume
  //
  G4Tubs* solidA = new G4Tubs("Target",0.,radius,sliceZ,0.,twopi);
  logicTarget = new G4LogicalVolume( solidA,targetMaterial,"Target");
  logicTarget->SetSensitiveDetector(targetSD);

  G4double z = sliceZ - targetZ;

  for(G4int i=0; i<nSlices; i++) {
    // physC = 
    new G4PVPlacement(0,G4ThreeVector(0.0,0.0,z),logicTarget,"Target",logicCheck,false,i);
    z += 2.0*sliceZ;
  }
 
  // colors
  G4VisAttributes zero = G4VisAttributes::Invisible;
  logicWorld->SetVisAttributes(&zero);

  G4VisAttributes regWcolor(G4Colour(0.3, 0.3, 0.3));
  logicCheck->SetVisAttributes(&regWcolor);

  G4VisAttributes regCcolor(G4Colour(0., 0.3, 0.7));
  logicTarget->SetVisAttributes(&regCcolor);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  return world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NA49DetectorConstruction::SetTargetMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(mat);

  if (material && material != targetMaterial) {
    targetMaterial = material;
    if(logicTarget) logicTarget->SetMaterial(targetMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NA49DetectorConstruction::SetWorldMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(mat);

  if (material && material != worldMaterial) {
    worldMaterial = material;
    if(logicWorld) logicWorld->SetMaterial(worldMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NA49DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NA49DetectorConstruction::SetTargetRadius(G4double val)  
{
  if(val > 0.0) {
    radius = val;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
