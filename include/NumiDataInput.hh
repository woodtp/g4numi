#include "globals.hh"
#include "G4ThreeVector.hh"

#ifndef NumiDataInput_h
#define NumiDataInput_h 1

class NumiDataInput
{
public:
  NumiDataInput();
  ~NumiDataInput();

  static NumiDataInput* GetNumiDataInput();
private:
  static NumiDataInput* fNumiDataInput;

public:

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

  //HPBaffle 
  G4int HPBaffle, HPBaffleGEANTMat;
  G4double HPBaffleX0,HPBaffleY0,HPBaffleZ0,HPBaffleDXDZ,HPBaffleDYDZ;
  G4double HPBaffleHeight,HPBaffleWidth,HPBaffleLength;
  G4double HPBaffleHoleHeight,HPBaffleHoleWidth;

  //Cooling pipes
  G4int NCPipeN;
  G4bool CPipeFilledWater[20];
  G4double CPipeX0[20],CPipeY0[20],CPipeZ0[20],CPipeLength[20],CPipeRadiusIn[20],CPipeRadiusOut[20],CPipeWallThick[20];
  G4double CPipeDXDZ[20],CPipeDYDZ[20];
  G4int CPGeantMat[20];
  G4double CPipeCurvRad[20],CPipeOpenAng[20],CPipeCloseAng[20];
  
  //Container
  G4int NContainerN;
  G4double CTubeZ0[20],CTubeLength[20],CTubeRin[20],CTubeRout[20];
  G4int CTubeGeantMat[20],CTubeHoleIndex[20];
  G4String CTubeVolName[20];

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

  // Blocks
  G4double BlockNblock;
  G4double BlockX0[20],BlockY0[20],BlockZ0[20], BlockDxdz[20],BlockDydz[20];
  G4double BlockLength[20],BlockHdx[20],BlockHdy[20];
  G4int BlockGeantmat[20];

  // Horn 1 & 2
  G4int PhornNphorn;

  G4double PhornZ1[8], PhornZ2[8];  
  G4int  PhornNpoint[8];  
  G4double PhornAin[8], PhornBin[8], PhornCin[8], PhornAout[8], PhornBout[8], PhornCout[8];
  G4double  PhornROCin[8], PhornROCout[8];
  G4double PhornThickFront[8], PhornThickEnd[8];
  
  G4double PhornX0[8], PhornY0[8], PhornZ0[8], PhornDXDZ[8], PhornDYDZ[8], PhornCurrent[8];
  G4int PhornGEANTmat[8]; 
  
  // Flux Area
  G4double xdet_near[10],ydet_near[10],zdet_near[10];
  G4double xdet_far[10],ydet_far[10],zdet_far[10];

  G4double HornCurrent; //old
  

};
#endif
