
#ifndef NumiDetectorConstruction_H
#define NumiDetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4BREPSolidPCone.hh"
#include "G4Transform3D.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class NumiDataInput;
class NumiMagneticField;
class NumiMagneticFieldIC;
class NumiMagneticFieldOC;

class NumiDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    NumiDetectorConstruction();
    ~NumiDetectorConstruction();

    G4VPhysicalVolume* Construct();

  private:

  NumiDataInput* NumiData;
  NumiMagneticField* numiMagField;
  NumiMagneticFieldIC* numiMagFieldIC;
  NumiMagneticFieldOC* numiMagFieldOC;
  
  void ConstructTargetHall();
  void ConstructTarget();
  void ConstructDecayPipe();
  void ConstructHadronAbsorber();  
  void ConstructHorns();  
  void DefineMaterials();
  void DestroyMaterials();
  G4Material* GetMaterial(G4int matcode);
  G4double phornRgivenZ(G4double a, G4double b, G4double c, G4double z);
 
    // Materials
  G4Material* Vacuum;
  G4Material* Air;
  G4Material* Water;
  G4Material* Be;
  G4Material* C;
  G4Material* Al;
  G4Material* Pb;
  G4Material* Fe;
  G4Material* CT852;
  G4Material* Concrete;
  G4Material* Target;

    // Logical volumes
    //
    G4LogicalVolume* ROCK_log;
    G4LogicalVolume* TRGT_lv;
    G4LogicalVolume* TUNE_log;
    G4LogicalVolume* BLK_log[20]; 
    G4LogicalVolume* TGAR_log;
    G4LogicalVolume* Horn_PM_lv[8];

    // Physical volumes
    //
    G4VPhysicalVolume* ROCK;
    G4VPhysicalVolume* TUNE;
    G4VPhysicalVolume* TGAR;
    G4VPhysicalVolume* TRGT;
    G4VPhysicalVolume* PHORN[8];
    G4VPhysicalVolume* CPIP[20];
    G4VPhysicalVolume* CNT[20];

    //Solids
    //
    G4Box* BLK_solid[20];
    G4BREPSolidPCone* Horn_PM[8];
};

#endif

