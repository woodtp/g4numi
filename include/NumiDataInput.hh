//----------------------------------------------------------------------
//
// $Id: NumiDataInput.hh,v 1.14 2008/05/01 18:09:46 loiacono Exp $
//----------------------------------------------------------------------

#ifndef NumiDataInput_h
#define NumiDataInput_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include "NumiHornSpiderSupport.hh"

typedef std::vector<G4double> vdouble_t;
typedef std::vector<G4int> vint_t;
typedef std::vector<G4String> vstring_t;
typedef std::vector<G4bool> vbool_t;
typedef std::vector<NumiHornSpiderSupport> vNumiHornSpiderSupport_t;

class NumiDataInput
{
public:
  NumiDataInput();
  ~NumiDataInput();

  static NumiDataInput* GetNumiDataInput();
private:
  static NumiDataInput* fNumiDataInput;

public:
  void SetTargetZ0(G4double val) {
    TargetZ0 = val;
    HPBaffleZ0 += (val + 0.45*m);
  }
  void SetHornCurrent(G4double val) {
    HornCurrent = val;
  }
  void SetDebugOn(G4bool val) {
    debugOn = val;
  }
  G4bool IsDebugOn(){
    return debugOn;
  }
  void SetNImpWeight(G4bool val) {
    NImpWeightOn = val;
  }
  void SetFlukaInput(G4bool val) {
    useFlukaInput = val;
  }
  void SetMarsInput(G4bool val) {
    useMarsInput = val;
  }
  void SetMuonInput(G4bool val) {
    useMuonInput = val;
  }
  void SetRunNumber(G4String runNum){
    RunNumber=runNum;
  }
  G4String GetRunNumber(){
    return RunNumber;
  }
  void SetGeometryTag(G4String geoName){
    geometry = geoName;
  }
  G4String GetGeometryTag(){
    return geometry;
  }
  G4double GetMaterialSigma(){
    return materialSigma;
  }
  void SetExtNtupleFileName(G4String fileName){
    extNtupleFileName=fileName;
  }
  G4String GetExtNtupleFileName(){
    return extNtupleFileName;
  }
  void SetNuNtupleName(G4String fileName){
    nuNtupleName=fileName;
  }
  void SetHadmmNtupleName(G4String fileName){
    hadmmNtupleName=fileName;
  }
  void SetHadmmNtupleDir(G4String fileDir){
    hadmmNtupleDir=fileDir;
  }
  G4String GetHadmmNtupleDir(){
    return hadmmNtupleDir;
  }
  void SetASCIIName(G4String fileName){
    asciiName=fileName;
  }
  void OutputNuNtuple(G4bool output){
    createNuNtuple=output;
  }
  void OutputHadmmNtuple(G4bool output){
    createHadmmNtuple=output;
  }
  void OutputASCII(G4bool output){
    createASCII=output;
  }
  G4bool GetKillTracking(){
    return KillTracking;
  }
  G4double GetKillTrackingThreshold(){
    return KillTrackingThreshold;
  }
  


 private:
  G4bool debugOn;
  G4String extNtupleFileName;

 public:
  G4bool NImpWeightOn, createNuNtuple,createHadmmNtuple, createASCII;
  G4bool useFlukaInput, useMarsInput, useMuonBeam, useMuonInput;
  G4bool KillTracking;
  G4String nuNtupleName, hadmmNtupleName, hadmmNtupleDir, asciiName, RunNumber, geometry;

  G4double protonMomentum, beamSigmaX, beamSigmaY, protonKineticEnergy, materialSigma;
  G4double KillTrackingThreshold;
  G4ThreeVector beamPosition, beamDirection;

  // World Volume
  G4double RockRadius, RockHalfLen, RockDensity, RockRadLen;

  // Target Area
  G4double TargetAreaZ0, TargetAreaLength, TargetAreaHeight,TargetAreaWidth;
  G4int TargetAreaGEANTmat;

 // added by Zachary Barnett - Target Area ConcretePit
  G4int NTHConcreteSectionsN;
  vdouble_t THConcreteX0, THConcreteY0, THConcreteZ0;
  vdouble_t THConcreteLength, THConcreteHdx, THConcreteHdy;
  vdouble_t THConcreteDxdz, THConcreteDydz;
  vint_t THConcreteGeantMaterial;
  vstring_t THConcreteName;
  //-------------------------------------
  
  // Target
  G4bool constructTarget;
  G4double TargetX0, TargetY0, TargetZ0, TargetDxdz, TargetDydz;
  G4double TargetSLength, TargetSWidth, TargetSHeight, TargetA, TargetDensity;
  G4double  TargetZ, TargetRL;
  G4int TargetGEANTmat,TargetSegmentNo;
  G4double TargetSegmentPitch,TargetCPGRadius,TargetCPGPosition;
  G4bool TargetEndRounded;
	
  //Rings holding target and cooling pipes
  G4int NTgtRingN;
  vdouble_t TgtRingZ0, TgtRingLength, TgtRingRin,TgtRingRout;
  vint_t TgtRingGeantMaterial;
  vstring_t TgtRingVolName;

  //Budal
  G4double BudalX0,BudalY0,BudalZ0,BudalDxdz,BudalDydz;

  //HPBaffle 
  G4int HPBaffle, HPBaffleGEANTMat;
  G4double HPBaffleX0,HPBaffleY0,HPBaffleZ0,HPBaffleDXDZ,HPBaffleDYDZ;
  G4double HPBaffleLength,HPBaffleRin,HPBaffleRout;

  //Cooling pipes
  G4int NCPipeN;
  vbool_t CPipeFilledWater;
  vdouble_t CPipeX0,CPipeY0,CPipeZ0,CPipeLength,CPipeRadiusIn,CPipeRadiusOut,CPipeWallThick;
  vdouble_t CPipeDXDZ,CPipeDYDZ;
  vint_t CPGeantMat;
  vdouble_t CPipeCurvRad,CPipeOpenAng,CPipeCloseAng;
  vstring_t CPipeVolName;
  
  //Container
  G4int NContainerN;
  vdouble_t CTubeZ0,CTubeLength,CTubeRin,CTubeRout;
  vint_t CTubeGeantMat;
  vstring_t CTubeVolName;

  // Tunnel
  G4double TunnelZ0, TunnelRadius, TunnelLength, TunnelA, TunnelZ;
  G4int TunnelGEANTmat;
  G4double BeamAngle;//nominally .058 radians

  // Decay Pipe
  G4double DecayPipeZ0,DecayPipeRadius,DecayPipeLength,DecayPipeFWinThick,DecayPipeEWinThick,DecayPipeWallThick;
  G4double DecayPipeA,DecayPipeZ;
  G4int DecayPipeGEANTmat,DecayPipeFWinmat, DecayPipeEWinmat;
  
  // Decay Pipe Shield
  G4double ShieldX0,ShieldY0,ShieldZ0,ShieldDxdz, ShieldDydz,ShieldLength,ShieldRout,ShieldRin;
  G4int ShieldGEANTmat;
  
  //Target Hall Blocks
  G4double THBlockNblock;
  vdouble_t THBlockX0,THBlockY0,THBlockZ0, THBlockDxdz,THBlockDydz;
  vdouble_t THBlockLength,THBlockHdx,THBlockHdy;
  vint_t THBlockGeantMaterial;
  vstring_t THBlockName;

  //Hadron Absorber Box
  G4double HadrBox_width, HadrBox_height, HadrBox_length;

  // Horn 1 & 2
  G4int PhornNphorn;

  vdouble_t PhornZ1, PhornZ2;  
  vint_t  PhornNpoint;  
  vdouble_t PhornAin, PhornBin, PhornCin, PhornAout, PhornBout, PhornCout;
  vdouble_t  PhornROCin, PhornROCout;
  vdouble_t PhornThickFront, PhornThickEnd;  
  vdouble_t PhornX0, PhornY0, PhornZ0, PhornDXDZ, PhornDYDZ, PhornCurrent;
  vint_t PhornGEANTmat; 
  vstring_t PhornName;
  
  // Flux Area
  G4int nNear,nFar;
  vdouble_t xdet_near,ydet_near,zdet_near;
  vdouble_t xdet_far,ydet_far,zdet_far;

  G4double HornCurrent; //old
	
 // Horn 1
  G4int NPHorn1OCN,NPHorn1ICN,NPHorn1EndN;
  
  vdouble_t PHorn1OCRout,PHorn1OCRin,PHorn1OCZ0;
  vdouble_t PHorn1ICZ0;
  vint_t PHorn1ICNpoint;

  vdouble_t PHorn1EndZ0,PHorn1EndLength,PHorn1EndRout,PHorn1EndRin;
  vint_t PHorn1EndGeantMat;
  vstring_t PHorn1EndVolName;

  G4int NHorn1SpiderSupportPlanesN,NHorn1SpidersPerPlaneN;
  vdouble_t Horn1SpiderSupportZ0;
  vNumiHornSpiderSupport_t Horn1SS;
  
  // Horn 2
  G4int NPHorn2OCN,NPHorn2ICN,NPHorn2EndN;
  
  vdouble_t PHorn2OCRout,PHorn2OCRin,PHorn2OCZ0;
  vdouble_t PHorn2ICZ0;
  vint_t PHorn2ICNpoint;

  vdouble_t PHorn2EndZ0,PHorn2EndLength,PHorn2EndRout,PHorn2EndRin;
  vint_t PHorn2EndGeantMat;
  vstring_t PHorn2EndVolName;

  G4int NHorn2SpiderSupportPlanesN,NHorn2SpidersPerPlaneN;
  vdouble_t Horn2SpiderSupportZ0;
  vNumiHornSpiderSupport_t Horn2SS;
};
#endif
