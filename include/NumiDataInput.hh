#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>

#ifndef NumiDataInput_h
#define NumiDataInput_h 1

typedef std::vector<G4double> vdouble_t;
typedef std::vector<G4int> vint_t;
typedef std::vector<G4String> vstring_t;
typedef std::vector<G4bool> vbool_t;

class NumiDataInput
{
public:
  NumiDataInput();
  ~NumiDataInput();

  static NumiDataInput* GetNumiDataInput();
private:
  static NumiDataInput* fNumiDataInput;

public:
 
  G4bool NImpWeightOn, CreateNuNtuple,CreateHadmmNtuple;

  G4double protonMomentum, beamSigmaX, beamSigmaY, protonKineticEnergy;
  G4ThreeVector beamPosition, beamDirection;

  // World Volume
  G4double RockRadius, RockHalfLen, RockDensity, RockRadLen;

  // Target Area
  G4double TargetAreaZ0, TargetAreaLength, TargetAreaHeight,TargetAreaWidth;
  G4int TargetAreaGEANTmat;
 
  // Target
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

  // Decay Pipe
  G4double DecayPipeZ0,DecayPipeRadius,DecayPipeLength,DecayPipeFWinThick,DecayPipeEWinThick,DecayPipeWallThick;
  G4double DecayPipeA,DecayPipeZ;
  G4int DecayPipeGEANTmat,DecayPipeFWinmat, DecayPipeEWinmat;
  
  // Decay Pipe Shield
  G4double ShieldX0,ShieldY0,ShieldZ0,ShieldDxdz, ShieldDydz,ShieldLength,ShieldRout,ShieldRin;
  G4int ShieldGEANTmat;
  
  //Blocks
  G4double BlockNblock;
  vdouble_t BlockX0,BlockY0,BlockZ0, BlockDxdz,BlockDydz;
  vdouble_t BlockLength,BlockHdx,BlockHdy;
  vint_t BlockGeantMaterial;
  vstring_t BlockName;
  vstring_t BlockMotherVolume;

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
  vdouble_t xdet_near,ydet_near,zdet_near;
  vdouble_t xdet_far,ydet_far,zdet_far;

  G4double HornCurrent; //old
  

};
#endif
