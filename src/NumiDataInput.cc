#include "NumiDataInput.hh"
#include "G4ThreeVector.hh"

NumiDataInput* NumiDataInput::fNumiDataInput = 0;

NumiDataInput* NumiDataInput::GetNumiDataInput()
{return fNumiDataInput;}

NumiDataInput::NumiDataInput()
{
  if (fNumiDataInput)
    { G4Exception("NumiDataInput constructed twice.");}
  fNumiDataInput = this;

  //Rock        1
  //=======================================================================
  RockRadius  = 10.0*m;
  RockHalfLen = 1200.0*m;
  RockDensity = 0.0; // not
  RockRadLen  = 0.0; // used

  //TargetArea          1
  //=======================================================================
  TargetAreaZ0       = -4.0*m;
  TargetAreaLength   = 49.28*m;
  TargetAreaHeight   = 6.0*m;
  TargetAreaWidth    = 6.0*m;
  TargetAreaGEANTmat = 15;

  // Target   1
  //=======================================================================
  TargetX0           = 0.0;
  TargetY0           = 0.0;
  TargetZ0           = -0.35*m;
  TargetDxdz         = 0.0; // doesn't
  TargetDydz         = 0.0; // work properly yet
  TargetSLength      = 20.*mm;
  TargetSWidth       = 6.40E-03*m;
  TargetSHeight      = 18.0E-03*m;
  TargetCPGRadius    = 3.2*mm; // Cooling pipe groove
  TargetCPGPosition  = 10.7*mm;
  TargetEndRounded   = true;
  TargetSegmentNo    = 47;
  TargetSegmentPitch = 0.3*mm;
  TargetA            = 12.01*g/mole;
  TargetZ            = 6.;
  TargetDensity      = 1.754*g/cm3;
  TargetRL           = 25.692;
  TargetGEANTmat     = 18;

  //HPBaffle           1 // Only the length and position can be changed currently
  //=======================================================================
  HPBaffleGEANTMat   =  18;
  HPBaffleX0         =  0.00;
  HPBaffleY0         =  0.00;
  HPBaffleZ0         = -3.04*m;
  HPBaffleDXDZ       =  0.0;
  HPBaffleDYDZ       =  0.0;
  HPBaffleHeight     =  0.10*m;
  HPBaffleWidth      =  0.10*m;
  HPBaffleLength     =  1.20*m;
  HPBaffleHoleHeight = 12.0E-3*m;
  HPBaffleHoleWidth  =  5.4E-3*m;

  //Cooling pipes
  NCPipeN = 9;
  //=======================================================================
  // CPipeDXDZ=-99999 and CPipeDYDZ=0 puts curved part in z-y plane
  //=======================================================================
  G4int CPGeantMat_[]        = {31     ,   31   ,10     , 10     , 10    , 10    ,  10   , 10    , 10 };
  G4bool CPipeFilledWater_[] = { true  , true   ,true   , true   , true  ,true   ,false  ,true   , false};
  G4double CPipeX0_[]        = {0.     ,0.      ,0.     , 0.     ,   0.  ,0      , 0.    , 0.    ,  0.};
  G4double CPipeY0_[]        = {1.05e-2,-1.05e-2,1.05e-2,-1.05e-2,3.5e-2 ,-3.5e-2, 0.    , 0.    ,  0.};
  G4double CPipeZ0_[]        = {-0.625 , -0.625 ,-0.65  ,-0.65   , -.6501,-.6501 ,.619   , .6195 , .6235};
  G4double CPipeDXDZ_[]      = {    0  ,  0     , 0     ,  0     ,-99999 ,-99999 , 0     ,  0    , 0}; 
  G4double CPipeDYDZ_[]      = {    0  ,  0     , 0     ,  0     , 0     ,0      , 0     ,  0    , 0};
  G4double CPipeLength_[]    = {1.230  , 1.230  , .025  , .025   , 0     ,0      , 5e-4  ,  4e-3 , 3e-3};
  G4double CPipeRadiusOut_[] = {3e-3   , 3e-3   ,3e-3   ,3e-3    , 3e-3  ,3e-3   ,13.52e-3,14.1e-3,14.1e-3};
  G4double CPipeRadiusIn_[]  = {   0.  ,    0.  ,   0.  ,   0.   , 0.    ,   0.  ,7.52e-3,7.3e-3,7.3e-3}; //This is not the inner pipe radius - use radiusOut and wallthickness for that
  G4double CPipeWallThick_[] = {4e-4   ,4e-4    ,4e-4   ,4e-4    ,4e-4   ,4e-4   , 0     ,1.6e-3 , 0}; 
  G4double CPipeCurvRad_[]   = {0      , 0      , 0     , 0      ,2.45e-2,2.45e-2, 0     ,   0   , 0}; // straight tube if this is 0;
  G4double CPipeOpenAng_[]   = {0      , 0      , 0     , 0      ,90.    ,90.    , 0     ,   0   , 0};
  G4double CPipeCloseAng_[]  = {0      , 0      , 0     , 0      ,270.   ,270.   , 0     ,   0   , 0};

  for (G4int ii=0;ii<NCPipeN;ii++){
    CPGeantMat[ii]= CPGeantMat_[ii];
    CPipeFilledWater[ii]=CPipeFilledWater_[ii];
    CPipeX0[ii] = CPipeX0_[ii]*m;
    CPipeY0[ii] = CPipeY0_[ii]*m;   
    CPipeZ0[ii] = CPipeZ0_[ii]*m;
    CPipeDXDZ[ii] = CPipeDXDZ_[ii];
    CPipeDYDZ[ii] = CPipeDYDZ_[ii];
    CPipeLength[ii] = CPipeLength_[ii]*m;
    CPipeRadiusOut[ii] = CPipeRadiusOut_[ii]*m;
    CPipeRadiusIn[ii] = CPipeRadiusIn_[ii]*m;
    CPipeWallThick[ii] = CPipeWallThick_[ii]*m;
    CPipeCurvRad[ii] = CPipeCurvRad_[ii]*m;
    CPipeOpenAng[ii] = CPipeOpenAng_[ii]*deg;
    CPipeCloseAng[ii] = CPipeCloseAng_[ii]*deg;
  }
  
  //Container
  //=======================================================================
  NContainerN=16;
  
  G4double CTubeZ0_[]     ={0.6285  , -0.3525 , -0.338 , -0.38   , -0.138   , 0.105    , 0.347   , 0.589  ,-0.3992   ,-0.4072, -.71481,-.714556,-.70855,-0.7022,-.6832  ,-0.4302  };
  G4double CTubeLength_[] ={5e-4    , 14.5e-3 , 0.9785 , 8e-3    ,  8e-3    , 8e-3     , 8e-3    , 8e-3   , 17e-3    ,8e-3   , .254e-3,6.006e-3,6.35e-3 , 19e-3, 253e-3 , 23e-3 };
  G4double CTubeRin_[]    ={0.0     , 16e-3   ,14.6e-3 , 9.5e-3  , 9.5e-3   , 9.5e-3   , 9.5e-3  , 9.5e-3 ,24e-3     , 74e-3 , 0.     ,1.27e-2,22.22e-3 ,17.5e-3,77e-3  , 77e-3};
  G4double CTubeRout_[]   ={14.15e-3,16.4e-3  , 15e-3  , 14.55e-3, 14.55e-3 , 14.55e-3 , 14.55e-3, 14.55e-3,106e-3   ,106e-3 , 22e-3  ,34.67e-3,34.67e-3,107e-3 , 83e-3 , 120e-3  };
  G4int CTubeHoleIndex_[] ={0       ,  0      , 0       , 1       , 1        , 1        , 1       , 1       , 0      , 0     ,  0    ,   0    ,  0     ,  0    ,  0    , 0};
  G4int CTubeGeantMat_[]  ={5       ,  9      , 9       , 9       , 9        , 9        , 9       , 9       , 10     , 10    ,  5    ,  10    , 10     , 10    , 10    , 10};
  G4String CTubeVolName_[]={"BeDW"  ,"AlTube1","AlTube2", "Ring1" ,"Ring2"   ,"Ring3"   ,"Ring4"  , "Ring5" ,"CLid2" ,"CLid1","BeUp1", "BeUp2","BeUp3" ,"BFront","Body","BEnd" };

 for (G4int ii=0;ii<NContainerN;ii++){
    CTubeZ0[ii]=CTubeZ0_[ii]*m;
    CTubeLength[ii]=CTubeLength_[ii]*m;
    CTubeRin[ii]=CTubeRin_[ii]*m;
    CTubeRout[ii]=CTubeRout_[ii]*m;
    CTubeHoleIndex[ii]=CTubeHoleIndex_[ii];
    CTubeGeantMat[ii]= CTubeGeantMat_[ii];
    CTubeVolName[ii]= CTubeVolName_[ii];
 }

  //Tunnel         1
  //=======================================================================
  TunnelZ0       = 45.28*m;
  TunnelRadius   = 3.3*m;
  TunnelLength   = 693.86*m;
  TunnelA        = 0.0;
  TunnelZ        = 0.0;
  TunnelGEANTmat = 15;

  //ShieldNshield  5 
  //======================================================================= 
  ShieldX0       = 0.0; 
  ShieldY0       = 0.0; 
  ShieldZ0       = 45.28*m;   
  ShieldDxdz     = 0.0; // not 
  ShieldDydz     = 0.0; // used
  ShieldLength   = 677.1*m;  
  ShieldRout     = 2.23*m; 
  ShieldRin      = 1.0097*m; 
  ShieldGEANTmat = 17;
  
  //DecayPipe          1    
  //=======================================================================
  DecayPipeZ0        = 45.28*m;
  DecayPipeRadius    = 0.9906*m;
  DecayPipeLength    = 677.1*m;
  DecayPipeFWinThick = 1.60E-3*m;
  DecayPipeEWinThick = 4.76E-3*m;
  DecayPipeWallThick = 1.905E-2*m;
  DecayPipeA         = 55.85;
  DecayPipeZ         = 26.0;
  DecayPipeGEANTmat  = 10;
  DecayPipeFWinmat   = 9;
  DecayPipeEWinmat   = 10;

  BlockNblock = 14;
  //======================================================================= 
  // Target area shielding (1-6) ; Absorber (7-13) ; Budal monitor (14) 
 
  G4double BlockX0_[]    = { 1.3970, -1.3970,  0.    ,  0.    , 0     , 0     , 1.651 , -1.3208, 0.0     , 0.3302  , 0.0   ,  0.0   , 0.0   , -7.5E-3 };
  G4double BlockY0_[]    = { 0.8707,  0.8707, -1.4153,  1.3665, 1.3792, 1.3665, 0.0   ,  0.0   , 0.0     , 0.0     , 0.0   , -1.651 , 1.651 ,  0.0    };
  G4double BlockZ0_[]    = {-4.0   , -4.0   , -4.0   , -4.0   , 9.0   , 15.0  , 723.48,  723.48, 726.4772, 728.7632, 723.48,  723.48, 723.48,  -51.72E-2};
  G4double BlockDxdz_[]  = { 0.0   ,  0.0   ,  0.0   ,  0.0   , 0     , 0     , 0.0   ,  0.0   , 0.0     , 0.0     , 0.0   ,  0.0   , 0.0   ,  0.0 };
  G4double BlockDydz_[]  = { 0.0   ,  0.0   ,  0.0   ,  0.0   , 0     , 0     , 0.0   ,  0.0   , 0.0     , 0.0     , 0.0   ,  0.0   , 0.0   ,  0.0 };
  G4double BlockLength_[]= { 49.28 ,  49.28 ,  49.28 ,  13.0  , 6.0   , 30.279, 5.2832,  5.2832, 2.2860  , 0.9144  , 2.4384,  2.9972, 2.9972,  20.0E-3 };
  G4double BlockHdx_[]   = { 0.8128,  0.8128,  2.2098,  0.5842, 0.5842, 0.5842, 0.9906,  0.6604, 0.6604  , 2.3114  , 0.6477,  0.6604, 0.6604,  17.5E-3 };
  G4double BlockHdy_[]   = { 1.8232,  1.8232,  0.4628,  1.0109, 0.9982, 1.0109, 2.6416,  2.6416, 2.6416  , 2.6416  , 0.6477,  0.9906, 0.9906,  3.2E-3 };
  G4int BlockGeantmat_[] = { 10    ,  10    ,  10    ,  10    , 10    , 10    , 10    ,  10    , 10      , 17      , 9     ,  10    , 10    ,  18 };

  for (G4int ii=0;ii<BlockNblock;ii++){
    BlockX0[ii]=BlockX0_[ii]*m;
    BlockY0[ii]=BlockY0_[ii]*m;
    BlockZ0[ii]=BlockZ0_[ii]*m;
    BlockDxdz[ii]=BlockDxdz_[ii];
    BlockDydz[ii]=BlockDydz_[ii];
    BlockLength[ii]= BlockLength_[ii]*m;
    BlockHdx[ii]= BlockHdx_[ii]*m;
    BlockHdy[ii]=BlockHdy_[ii]*m;
    BlockGeantmat[ii]=BlockGeantmat_[ii];
  }

  HornCurrent=205; //kA
  PhornNphorn=8;
  //=======================================================================
  // Horn 1 & 2
  //=======================================================================
  G4double PhornZ1_[]         = {   0.00000 ,  0.44047  ,  0.80000 ,  0.83982  ,  0.95128 , 0.00000  , 0.97617 ,   1.04803};
  G4double PhornZ2_[]         = {   0.44047 ,  0.80000  ,  0.83982 ,  0.95128  ,  3.00000 , 0.97617  , 1.04803 ,   3.00000};

  G4int PhornNpoint_[]        = {       10  ,       10  ,       1  ,        5  ,       10 ,       10 ,       1 ,        10};
  G4double PhornAin_[]        = {  92.8454  ,  85.7091  ,  0.0000  , -82.2123  , -80.0000 , 100.0000 , 0.0000  , -100.0000};
  G4double PhornBin_[]        = {   7.0483  ,   7.0483  ,  0.0000  ,  2.1850   ,  2.1850  ,  0.1351  , 0.0000  ,    0.2723};
  G4double PhornCin_[]        = {  -0.2000  ,   0.0000  ,  0.9000  ,  0.0000   , -0.2000  , -0.3000  , 3.9000  ,   -0.3000};
  G4double PhornAout_[]       = {  92.8454  ,  92.8454  ,  0.0000  , -80.0000  , -80.0000 , 100.0000 , 0.0000  , -100.0000};
  G4double PhornBout_[]       = {   7.0483  ,   7.0483  ,  0.0000  ,  2.1850   ,  2.1850  ,  0.1351  , 0.0000  ,    0.2723};
  G4double PhornCout_[]       = {   0.0000  ,   0.0000  ,  1.3500  ,  0.0000   , 0.0000   ,  0.0000  , 4.4000  ,    0.0000};
  
  G4double PhornROCin_[]      = {   0.1533  ,   0.1533  ,  0.1533  ,  0.1533   , 0.1533   ,  0.3700  , 0.3700  ,    0.3700};
  G4double PhornROCout_[]     = {   0.1620  ,   0.1620  ,  0.1620  ,  0.1620   , 0.1620   ,  0.3787  , 0.3787  ,    0.3787};
  
  G4double PhornThickFront_[] = {   0.0030  ,   0.0000  ,  0.0000  ,  0.0000   , 0.0000   ,  0.0030  , 0.0000  ,    0.0000};
  G4double PhornThickEnd_[]   = {   0.0000  ,   0.0000  ,  0.0000  ,  0.0000   , 0.0030   ,  0.0000  , 0.0000  ,    0.0030};

  G4double PhornX0_[]         = {   0.0000  ,   0.0000  ,  0.0000  ,  0.0000   , 0.0000   ,  0.0000  , 0.0000  ,  0.0000  };
  G4double PhornY0_[]         = {   0.0000  ,   0.0000  ,  0.0000  ,  0.0000   , 0.0000   ,  0.0000  , 0.0000  ,  0.0000  };
  G4double PhornZ0_[]         = {   0.0000  ,   0.0000  ,  0.0000  ,  0.0000   , 0.0000   , 10.0000  , 10.0000 , 10.0000  };
  G4double PhornDXDZ_[]       = {   0.0000  ,   0.0000  ,  0.0000  ,  0.0000   , 0.0000   ,  0.0000  , 0.0000  ,  0.0000  };
  G4double PhornDYDZ_[]       = {   0.0000  ,   0.0000  ,  0.0000  ,  0.0000   , 0.0000   ,  0.0000  , 0.0000  ,  0.0000  };
  G4double PhornCurrent_[]    = {   0.2000  ,   0.2000  ,  0.2000  ,  0.2000   , 0.2000   ,  0.2000  , 0.2000  ,  0.2000  };
  G4int PhornGEANTmat_[]      = {       20  ,       20  ,      20  ,      20   ,     20   ,  20      , 20      ,  20      };
  
  for (G4int ii=0;ii<PhornNphorn;ii++){
    PhornZ1[ii]         = PhornZ1_[ii]*m;
    PhornZ2[ii]         = PhornZ2_[ii]*m;  
    
    PhornNpoint[ii]     = PhornNpoint_[ii];
    PhornAin[ii]        = PhornAin_[ii]*cm; 
    PhornBin[ii]        = PhornBin_[ii]/cm;
    PhornCin[ii]        = PhornCin_[ii]*cm;   
    PhornAout[ii]       = PhornAout_[ii]*cm;   
    PhornBout[ii]       = PhornBout_[ii]/cm;
    PhornCout[ii]       = PhornCout_[ii]*cm;
    
    PhornROCin[ii]      = PhornROCin_[ii]*m;  
    PhornROCout[ii]     = PhornROCout_[ii]*m;
    
    PhornThickFront[ii] = PhornThickFront_[ii]*m;
    PhornThickEnd[ii]   = PhornThickEnd_[ii]*m;
    
    PhornX0[ii]         = PhornX0_[ii]*m;    
    PhornY0[ii]         = PhornY0_[ii]*m;  
    PhornZ0[ii]         = PhornZ0_[ii]*m;  
    PhornDXDZ[ii]       = PhornDXDZ_[ii];   
    PhornDYDZ[ii]       = PhornDYDZ_[ii]; 
    PhornCurrent[ii]    = PhornCurrent_[ii];
    PhornGEANTmat[ii]   = PhornGEANTmat_[ii];
  }

 //Near & Far Detector location
  xdet_near[0]=0;
  ydet_near[0]=0;
  zdet_near[0]=1040.*m;
  xdet_far[0]=0;
  ydet_far[0]=0;
  zdet_far[0]=735000.*m;

  for(G4int ii=1;ii<10;ii++){
  xdet_near[ii]=0;
  ydet_near[ii]=0;
  zdet_near[ii]=0;
  xdet_far[ii]=0;
  ydet_far[ii]=0;
  zdet_far[ii]=0;
  }
  
}
NumiDataInput::~NumiDataInput()
{
}
