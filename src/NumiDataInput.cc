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


  NImpWeightOn = false; CreateNuNtuple=true;CreateHadmmNtuple=false;
  protonMomentum = 120.*GeV;
  beamSigmaY = 1.*mm;
  beamSigmaX = 0.9*mm;
  beamDirection = G4ThreeVector(0.,0.,1.);
  beamPosition = G4ThreeVector(0.,0.,-4.*m);

  protonKineticEnergy = sqrt(pow((.938*GeV),2)+pow(protonMomentum,2))-0.938*GeV;

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

  //Budal Monitor
  BudalX0 = 0.0;
  BudalY0 = 0.0;
  BudalZ0 = -16.72*cm;
  BudalDxdz = 0.0;
  BudalDydz = 0.0;
  
  //HPBaffle           1 // Only the length and position can be changed currently
  //=======================================================================
  HPBaffleGEANTMat   =  18;
  HPBaffleX0         =  0.00;
  HPBaffleY0         =  0.00;
  HPBaffleZ0         = -3.04*m;
  HPBaffleDXDZ       =  0.0;
  HPBaffleDYDZ       =  0.0;
  HPBaffleLength     =  1.20*m;
  HPBaffleRin        =  5.5*mm;
  HPBaffleRout       =  3.*cm;
 
  /*
  //Cooling pipes
  NCPipeN = 9;
  //=======================================================================
  // CPipeDXDZ=-99999 and CPipeDYDZ=0 puts curved part in z-y plane
  // Z0 is with respect to the first target fin 
  //=======================================================================
  G4int CPGeantMat_[]        = {31     ,   31   ,10     , 10     , 10    , 10    ,  10   , 10    , 10 };
  G4bool CPipeFilledWater_[] = { true  , true   ,true   , true   , true  ,true   ,false  ,true   , false};
  G4double CPipeX0_[]        = {0.     ,0.      ,0.     , 0.     ,   0.  ,0.     , 0.    , 0.    ,  0.};
  G4double CPipeY0_[]        = {1.05e-2,-1.05e-2,1.05e-2,-1.05e-2,3.5e-2 ,-3.5e-2, 0.    , 0.    ,  0.};
  G4double CPipeZ0_[]        = {-0.275 , -0.275 ,-0.3  ,-0.3     , -.301 ,-.301  ,.969   , .9695 , .9735};
  G4double CPipeDXDZ_[]      = {    0  ,  0     , 0     ,  0     ,-99999 ,-99999 , 0     ,  0    , 0}; 
  G4double CPipeDYDZ_[]      = {    0  ,  0     , 0     ,  0     , 0     ,0      , 0     ,  0    , 0};
  G4double CPipeLength_[]    = {1.230  , 1.230  , .025  , .025   , 0     ,0      , 5e-4  ,  4e-3 , 3e-3};
  G4double CPipeRadiusOut_[] = {3e-3   , 3e-3   ,3e-3   ,3e-3    , 3e-3  ,3e-3   ,13.52e-3,14.1e-3,14.1e-3};
  G4double CPipeRadiusIn_[]  = {   0.  ,    0.  ,   0.  ,   0.   , 0.    ,   0.  ,7.52e-3,7.3e-3,7.3e-3}; //This is not the inner pipe radius - use radiusOut and wallthickness for that
  G4double CPipeWallThick_[] = {4e-4   ,4e-4    ,4e-4   ,4e-4    ,4e-4   ,4e-4   , 0     ,1.6e-3 , 0}; 
  G4double CPipeCurvRad_[]   = {0      , 0      , 0     , 0      ,2.45e-2,2.45e-2, 0     ,   0   , 0}; // straight tube if this is 0;
  G4double CPipeOpenAng_[]   = {0      , 0      , 0     , 0      ,90.    ,90.    , 0     ,   0   , 0};
  G4double CPipeCloseAng_[]  = {0      , 0      , 0     , 0      ,270.   ,270.   , 0     ,   0   , 0};
  G4String CPipeVolName_[]   = {"Pipe1","Pipe2" ,"Pipe3","Pipe4" ,"Pipe5","Pipe6","Pipe7","Pipe8","Pipe9"};
  */
  //Cooling pipes
  NCPipeN = 17;
  //=======================================================================
  // CPipeDXDZ=-99999 and CPipeDYDZ=0 puts curved part in z-y plane
  //=======================================================================
  G4int CPGeantMat_[]        = {10     ,   10   ,10     , 10     , 10    , 10    ,  10   , 10    , 10     , 10           ,10            , 31     , 31      ,  10     ,  10     , 10      ,10  };
  G4bool CPipeFilledWater_[] = { true  , true   ,true   , true   , true  ,true   ,false  ,true   , false  , true         ,true          , true   ,true     , true    , true    , true    , true};
  G4double CPipeX0_[]        = {0.     ,0.      ,0.     , 0.     ,   0.  ,0.     , 0.    , 0.    ,  0.    , 0.           ,0.            , 0      , 0       , 0       , 0       , 0       , 0};
  G4double CPipeY0_[]        = {1.05e-2,-1.05e-2,1.05e-2,-1.05e-2,3.5e-2 ,-3.5e-2, 0.    , 0.    ,  0.    , 1.05e-2      ,-1.05e-2      , 5.95e-2, -5.95e-2, 5.95e-2 , -5.95e-2,5.95e-2  ,-5.95e-2};
  G4double CPipeZ0_[]        = {-0.275 , -0.275 ,-0.30  ,-0.30   , -.3001,-.3001 ,.969   , .9695 , .9735  , 0.955        ,0.955         , -0.151 , -0.151  ,-0.287   ,-0.287   , -.30    , -.30};
  G4double CPipeDXDZ_[]      = {    0  ,  0     , 0     ,  0     ,-99999 ,-99999 , 0     ,  0    , 0      , 0            ,0             , 0.     , 0.      , 0.      , 0.      , 0.      , 0.}; 
  G4double CPipeDYDZ_[]      = {    0  ,  0     , 0     ,  0     , 0     ,0      , 0     ,  0    , 0      , 0            ,0             , 0.     , 0.      , 0.      , 0.      , 0.      , 0.};
  G4double CPipeLength_[]    = {1.230  , 1.230  , .025  , .025   , 0     ,0      , 5e-4  ,  4e-3 , 3e-3   , 14.e-3      ,14.e-3         , 3e-2   , 3e-2    , 13.6e-2 ,13.6e-2  , 13e-3   ,13e-3};//PipeAdapter off by 4.5mm???
  G4double CPipeRadiusOut_[] = {3e-3   , 3e-3   ,3e-3   ,3e-3    , 3e-3  ,3e-3   ,13.52e-3,14.1e-3,14.1e-3, 3.4e-3       , 3.4e-3       , 5e-3   , 5e-3    , 5.2e-3   , 5.2e-3 , 3e-3    ,3e-3};
  G4double CPipeRadiusIn_[]  = {   0.  ,    0.  ,   0.  ,   0.   , 0.    ,   0.  ,7.52e-3,7.3e-3,7.3e-3   , 0.           , 0.           , 0      , 0       , 0      , 0        , 0.      ,0.}; 
  //CPipeRadiusIn is not the inner pipe radius - use radiusOut and wallthickness for that
  G4double CPipeWallThick_[] = {4e-4   ,4e-4    ,4e-4   ,4e-4    ,4e-4   ,4e-4   , 0     ,1.6e-3 , 0      , 0.8e-3       , 0.8e-3       , 1.e-3  , 1.e-3   , 1.1e-3  , 1.1e-3  , 4e-4    , 4e-4}; 
  G4double CPipeCurvRad_[]   = {0      , 0      , 0     , 0      ,2.45e-2,2.45e-2, 0     ,   0   , 0      , 0            , 0            , 0      , 0       , 0      , 0        , 0       , 0}; // straight tube if this is 0;
  G4double CPipeOpenAng_[]   = {0      , 0      , 0     , 0      ,90.    ,90.    , 0     ,   0   , 0      , 0            , 0            , 0      , 0       , 0      , 0        , 0       , 0   };
  G4double CPipeCloseAng_[]  = {0      , 0      , 0     , 0      ,270.   ,270.   , 0     ,   0   , 0      , 0            , 0            , 0      , 0       , 0      , 0        , 0       , 0  };
  G4String CPipeVolName_[]   = {"Pipe1","Pipe2" ,"Pipe3","Pipe4" ,"Pipe5","Pipe6","Pipe7","Pipe8","Pipe9" ,"PipeAdapter1","PipeAdapter2","PipeC1","PipeC2" ,"PipeB1","PipeB2"  ,"Pipe1tp","Pipe2btm" };
  //For bellows (PipeB1 & 2) I don't have real dimensions
  for (G4int ii=0;ii<NCPipeN;ii++){
    CPGeantMat.push_back(CPGeantMat_[ii]);
    CPipeFilledWater.push_back(CPipeFilledWater_[ii]);
    CPipeX0.push_back(CPipeX0_[ii]*m);
    CPipeY0.push_back(CPipeY0_[ii]*m);   
    CPipeZ0.push_back(CPipeZ0_[ii]*m);
    CPipeDXDZ.push_back(CPipeDXDZ_[ii]);
    CPipeDYDZ.push_back(CPipeDYDZ_[ii]);
    CPipeLength.push_back(CPipeLength_[ii]*m);
    CPipeRadiusOut.push_back(CPipeRadiusOut_[ii]*m);
    CPipeRadiusIn.push_back(CPipeRadiusIn_[ii]*m);
    CPipeWallThick.push_back(CPipeWallThick_[ii]*m);
    CPipeCurvRad.push_back(CPipeCurvRad_[ii]*m);
    CPipeOpenAng.push_back(CPipeOpenAng_[ii]*deg);
    CPipeCloseAng.push_back(CPipeCloseAng_[ii]*deg);
    CPipeVolName.push_back(CPipeVolName_[ii]);
  }
  
  //Container
  //=======================================================================
  // Z0 with respect to the first target fin
  NContainerN=11;
  
  G4double CTubeZ0_[]     ={0.9785  , -0.0025 ,  0.012  , -0.0492 ,-0.0572, -.36481,-.364556,-.35855,-0.3522,-.3332  ,-0.0802  };
  G4double CTubeLength_[] ={5e-4    , 14.5e-3 , 0.9785  , 17e-3   ,8e-3   , .254e-3,6.006e-3,6.35e-3 , 19e-3, 253e-3 , 23e-3 };
  G4double CTubeRin_[]    ={0.0     , 16e-3   ,14.6e-3  , 24e-3   , 74e-3 , 0.     ,1.27e-2,22.22e-3 ,17.5e-3,77e-3  , 77e-3};
  G4double CTubeRout_[]   ={14.15e-3,16.4e-3  , 15e-3   , 106e-3  ,106e-3 , 22e-3  ,34.67e-3,34.67e-3,107e-3 , 83e-3 , 120e-3  };
  G4int CTubeGeantMat_[]  ={5       ,  9      , 9       , 10      , 10    ,  5    ,  10    , 10     , 10    , 10    , 10};
  G4String CTubeVolName_[]={"BeDW"  ,"AlTube1","AlTube2", "CLid2" ,"CLid1","BeUp1", "BeUp2","BeUp3" ,"BFront","Body","BEnd" };

 for (G4int ii=0;ii<NContainerN;ii++){
    CTubeZ0.push_back(CTubeZ0_[ii]*m);
    CTubeLength.push_back(CTubeLength_[ii]*m);
    CTubeRin.push_back(CTubeRin_[ii]*m);
    CTubeRout.push_back(CTubeRout_[ii]*m);
    CTubeGeantMat.push_back(CTubeGeantMat_[ii]);
    CTubeVolName.push_back(CTubeVolName_[ii]);
 }
 //Rings holding target and cooling pipes
 NTgtRingN = 5;
 G4double TgtRingZ0_[]         = {-0.03   ,  0.212   , 0.455    , 0.697   , 0.939 };
 G4double TgtRingLength_[]     = {8e-3    ,  8e-3    , 8e-3     , 8e-3    , 8e-3 };
 G4double TgtRingRin_[]        = {9.5e-3  , 9.5e-3   , 9.5e-3   , 9.5e-3  , 9.5e-3};
 G4double TgtRingRout_[]          = {14.55e-3, 14.55e-3 , 14.55e-3 , 14.55e-3, 14.55e-3};
 G4int TgtRingGeantMaterial_[] = {9       , 9        , 9        , 9       , 9 };
 G4String TgtRingVolName_[]    = {"Ring1" , "Ring2"  , "Ring3"  , "Ring4" , "Ring5" };

for (G4int ii=0;ii<NTgtRingN;ii++){
    TgtRingZ0.push_back(TgtRingZ0_[ii]*m);
    TgtRingLength.push_back(TgtRingLength_[ii]*m);
    TgtRingRin.push_back(TgtRingRin_[ii]*m);
    TgtRingRout.push_back(TgtRingRout_[ii]*m);
    TgtRingGeantMaterial.push_back(TgtRingGeantMaterial_[ii]);
    TgtRingVolName.push_back(TgtRingVolName_[ii]);
 }


  //Tunnel         1
  //=======================================================================
  TunnelZ0       = 45.28*m;
  TunnelRadius   = 3.3*m+.5*m; //added .5m because hadron absorber does not fit entirely inside 3.3
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

  BlockNblock = 13;
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
  G4int BlockGeantMaterial_[] = { 10,  10   ,  10    ,  10    , 10    , 10    , 10    ,  10    , 10      , 17      , 9     ,  10    , 10    ,  18 };
  G4String BlockName_[]  = {"BLK1" , "BLK2" ,"BLK3"  ,"BLK4"  ,"BLK5" ,"BLK6" ,"BLK7" , "BLK8" , "BLK9"  , "BLK10" ,"BLK11","BLK12" ,"BLK13","BLK14" };
  G4String BlockMotherVolume_[]={"TGAR","TGAR","TGAR","TGAR"  ,"TGAR" ,"TGAR" ,"TUNE" , "TUNE" , "TUNE"  , "TUNE"  ,"TUNE" ,"TUNE"  ,"TUNE" ,"TGAR" };

  for (G4int ii=0;ii<BlockNblock;ii++){
    BlockX0.push_back(BlockX0_[ii]*m);
    BlockY0.push_back(BlockY0_[ii]*m);
    BlockZ0.push_back(BlockZ0_[ii]*m);
    BlockDxdz.push_back(BlockDxdz_[ii]);
    BlockDydz.push_back(BlockDydz_[ii]);
    BlockLength.push_back(BlockLength_[ii]*m);
    BlockHdx.push_back(BlockHdx_[ii]*m);
    BlockHdy.push_back(BlockHdy_[ii]*m);
    BlockGeantMaterial.push_back(BlockGeantMaterial_[ii]);
    BlockName.push_back(BlockName_[ii]);
    BlockMotherVolume.push_back(BlockMotherVolume_[ii]);
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
    PhornZ1.push_back(PhornZ1_[ii]*m);
    PhornZ2.push_back(PhornZ2_[ii]*m);  
    
    PhornNpoint.push_back(PhornNpoint_[ii]);
    PhornAin.push_back(PhornAin_[ii]*cm); 
    PhornBin.push_back(PhornBin_[ii]/cm);
    PhornCin.push_back(PhornCin_[ii]*cm);   
    PhornAout.push_back(PhornAout_[ii]*cm);   
    PhornBout.push_back(PhornBout_[ii]/cm);
    PhornCout.push_back(PhornCout_[ii]*cm);
    
    PhornROCin.push_back(PhornROCin_[ii]*m);  
    PhornROCout.push_back(PhornROCout_[ii]*m);
    
    PhornThickFront.push_back(PhornThickFront_[ii]*m);
    PhornThickEnd.push_back(PhornThickEnd_[ii]*m);
    
    PhornX0.push_back(PhornX0_[ii]*m);    
    PhornY0.push_back(PhornY0_[ii]*m);  
    PhornZ0.push_back(PhornZ0_[ii]*m);  
    PhornDXDZ.push_back(PhornDXDZ_[ii]);   
    PhornDYDZ.push_back(PhornDYDZ_[ii]); 
    PhornCurrent.push_back(PhornCurrent_[ii]);
    PhornGEANTmat.push_back(PhornGEANTmat_[ii]);
  }

 //Near & Far Detector location
  xdet_near.push_back(0);
  ydet_near.push_back(0);
  zdet_near.push_back(1040.*m);
  xdet_far.push_back(0);
  ydet_far.push_back(0);
  zdet_far.push_back(735000.*m);

  for(G4int ii=1;ii<10;ii++){
  xdet_near.push_back(0);
  ydet_near.push_back(0);
  zdet_near.push_back(0);
  xdet_far.push_back(0);
  ydet_far.push_back(0);
  zdet_far.push_back(0);
  }
  
}
NumiDataInput::~NumiDataInput()
{
}
