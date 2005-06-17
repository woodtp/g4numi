#include "NumiDataInput.hh"
#include "G4ThreeVector.hh"
#include "NumiHornSpiderSupport.hh"

NumiDataInput* NumiDataInput::fNumiDataInput = 0;

NumiDataInput* NumiDataInput::GetNumiDataInput()
{return fNumiDataInput;}

NumiDataInput::NumiDataInput()
{
  if (fNumiDataInput)
    { G4Exception("NumiDataInput constructed twice.");}
  fNumiDataInput = this;


  NImpWeightOn = true; CreateNuNtuple=true;CreateHadmmNtuple=false;
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
  TargetZ0           = -0.35*m-10.*cm;
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
  HPBaffleZ0         = -3.14*m;
  HPBaffleDXDZ       =  0.0;
  HPBaffleDYDZ       =  0.0;
  HPBaffleLength     =  1.20*m;
  HPBaffleRin        =  5.5*mm;
  HPBaffleRout       =  3.*cm;
 
  //Cooling pipes
  NCPipeN = 19;
  //=======================================================================
  // CPipeDXDZ=-99999 and CPipeDYDZ=0 puts curved part in z-y plane
  //=======================================================================
  G4int CPGeantMat_[]        = {10     ,   10   ,10     , 10     , 10    , 10    ,  10   , 10    , 10     , 10           ,10            , 31     , 31      ,10        ,10        ,  10     ,  10     , 10      ,10  };
  G4bool CPipeFilledWater_[] = { true  , true   ,true   , true   , true  ,true   ,false  ,true   , false  , true         ,true          , true   ,true     ,true      ,true      , true    , true    , true    , true};
  G4double CPipeX0_[]        = {0.     ,0.      ,0.     , 0.     ,   0.  ,0.     , 0.    , 0.    ,  0.    , 0.           ,0.            , 0      , 0       ,0         ,0         , 0       , 0       , 0       , 0};
  G4double CPipeY0_[]        = {1.05e-2,-1.05e-2,1.05e-2,-1.05e-2,3.5e-2 ,-3.5e-2, 0.    , 0.    ,  0.    , 1.05e-2      ,-1.05e-2      , 5.95e-2, -5.95e-2,5.95e-2   ,-5.95e-2  , 5.95e-2 , -5.95e-2,5.95e-2  ,-5.95e-2};
  G4double CPipeZ0_[]        = {-0.275 , -0.275 ,-0.30  ,-0.30   , -.3001,-.3001 ,.969   , .9695 , .9735  , 0.955        ,0.955         , -0.215 , -0.215  ,-.071     ,-.071     ,-0.287   ,-0.287   , -.30    , -.30};
  G4double CPipeDXDZ_[]      = {    0  ,  0     , 0     ,  0     ,-99999 ,-99999 , 0     ,  0    , 0      , 0            ,0             , 0.     , 0.      ,0.        ,0         , 0.      , 0.      , 0.      , 0.}; 
  G4double CPipeDYDZ_[]      = {    0  ,  0     , 0     ,  0     , 0     ,0      , 0     ,  0    , 0      , 0            ,0             , 0.     , 0.      ,0.        ,0         , 0.      , 0.      , 0.      , 0.};
  G4double CPipeLength_[]    = {1.230  , 1.230  , .025  , .025   , 0     ,0      , 5e-4  ,  4e-3 , 3e-3   , 14.e-3      ,14.e-3         , 14.4e-2, 14.4e-2 ,2e-2      ,2e-2      , 7.2e-2      ,7.2e-2  , 13e-3   ,13e-3};
  G4double CPipeRadiusOut_[] = {3e-3   , 3e-3   ,3e-3   ,3e-3    , 3e-3  ,3e-3   ,13.52e-3,14.1e-3,14.1e-3, 3.4e-3       , 3.4e-3       , 5e-3   , 5e-3    ,4e-3      ,4e-3      , 7.5e-3      , 7.5e-3 , 3e-3    ,3e-3};
  G4double CPipeRadiusIn_[]  = {   0.  ,    0.  ,   0.  ,   0.   , 0.    ,   0.  ,7.52e-3,7.3e-3,7.3e-3   , 0.           , 0.           , 0      , 0       ,0         ,0         , 0           , 0        , 0.      ,0.}; 
  //CPipeRadiusIn is not the inner pipe radius - use radiusOut and wallthickness for that
  G4double CPipeWallThick_[] = {4e-4   ,4e-4    ,4e-4   ,4e-4    ,4e-4   ,4e-4   , 0     ,1.6e-3 , 0      , 0.8e-3       , 0.8e-3       , 1.e-3  , 1.e-3   ,1e-3      ,1e-3      , 1.1e-3      , 1.1e-3  , 4e-4    , 4e-4}; 
  G4double CPipeCurvRad_[]   = {0      , 0      , 0     , 0      ,2.45e-2,2.45e-2, 0     ,   0   , 0      , 0            , 0            , 0      , 0       ,0         ,0         , 0           , 0        , 0       , 0}; // straight tube if this is 0;
  G4double CPipeOpenAng_[]   = {0      , 0      , 0     , 0      ,90.    ,90.    , 0     ,   0   , 0      , 0            , 0            , 0      , 0       ,0         ,0         , 0           , 0        , 0       , 0   };
  G4double CPipeCloseAng_[]  = {0      , 0      , 0     , 0      ,270.   ,270.   , 0     ,   0   , 0      , 0            , 0            , 0      , 0       ,0         ,0         , 0           , 0        , 0       , 0  };
  G4String CPipeVolName_[]   = {"Pipe1","Pipe2" ,"Pipe3","Pipe4" ,"Pipe5","Pipe6","Pipe7","Pipe8","Pipe9" ,"PipeAdapter1","PipeAdapter2","PipeC1","PipeC2" ,"PipeEndT","PipeEndB","PipeBellowT","PipeBellowB"  ,"Pipe1tp","Pipe2btm" };
  
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
  NContainerN=17;
  
  G4double CTubeZ0_[]     ={-.36481-.057, -.364556-.057, -.41555, -.35855 , -0.3522 ,-.3332 ,-0.0802, -0.123  ,-0.109  ,-0.0572 ,-0.0492 , -0.04  ,-0.042 ,-.022 , -0.013   , -0.0025  ,  0.011    , 0.9785  };
  G4double CTubeLength_[] ={.254e-3, 6.006e-3, 5.7e-2 , 6.35e-3 , 19e-3   , 253e-3, 23e-3 , 11e-2   , 66e-3  ,8e-3   , 17e-3  , 6e-3   ,0.02   ,15e-3  , 24e-3    , 14.5e-3  , 0.9785    , 5e-4 };
  G4double CTubeRin_[]    ={  0.   , 1.27e-2 , 1.7e-2 , 22.22e-3, 17.5e-3 , 77e-3 , 77e-3 , 14.6e-3 , 16.e-3 , 74e-3 , 24e-3  , 18.5e-3,16e-3  ,16e-3  , 14.6e-3  , 16e-3    , 14.6e-3   , 0.0  };
  G4double CTubeRout_[]   ={22e-3  , 34.67e-3, 1.9e-2 , 34.67e-3, 107e-3  , 83e-3 , 120e-3, 16.0e-3 , 21.4e-3,106e-3 , 106e-3 , 24e-3  ,18.5e-3,18.5e-3, 16e-3    , 16.4e-3  , 15e-3     , 14.15e-3  };
  G4int CTubeGeantMat_[]  ={   5   ,  10     ,   10   , 10      , 10      , 10    , 10    ,  31     , 10     , 10    , 10     ,  10    , 10    ,10     , 10       , 9        , 9         , 5 };
  G4String CTubeVolName_[]={"BeUp1", "BeUp2" , "Added", "BeUp3" , "BFront", "Body", "BEnd","CerTube", "Conn1","CLid1", "CLid2", "Conn2","Conn3","Tube1a" ,"Tube1b", "AlTube1", "AlTube2" , "BeDW"};

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
 G4double TgtRingRout_[]       = {14.55e-3, 14.55e-3 , 14.55e-3 , 14.55e-3, 14.55e-3};
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

  THBlockNblock = 6;
  //======================================================================= 
  // Target Hall shielding (1-6) ;
  G4double THBlockX0_[]    = { 1.3970, -1.3970,  0.    ,  0.    , 0     , 0     };
  G4double THBlockY0_[]    = { 0.8707,  0.8707, -1.4153,  1.3665, 1.3792+0.0325, 1.3665};
  G4double THBlockZ0_[]    = {-4.0   , -4.0   , -4.0   , -4.0   , 9.0   , 15.0  };
  G4double THBlockDxdz_[]  = { 0.0   ,  0.0   ,  0.0   ,  0.0   , 0     , 0     };
  G4double THBlockDydz_[]  = { 0.0   ,  0.0   ,  0.0   ,  0.0   , 0     , 0     };
  G4double THBlockLength_[]= { 49.28 ,  49.28 ,  49.28 ,  13.0  , 6.0   , 30.279};
  G4double THBlockHdx_[]   = { 0.8128,  0.8128,  2.2098,  0.5842, 0.5842, 0.5842};
  G4double THBlockHdy_[]   = { 1.8232,  1.8232,  0.4628,  1.0109, 0.9982-0.065, 1.0109};
  G4int THBlockGeantMaterial_[] = { 10,  10   ,  10    ,  10    , 10    , 10   };
  G4String THBlockName_[]  = {"BLK1" , "BLK2" ,"BLK3"  ,"BLK4"  ,"BLK5" ,"BLK6"};

  for (G4int ii=0;ii<THBlockNblock;ii++){
    THBlockX0.push_back(THBlockX0_[ii]*m);
    THBlockY0.push_back(THBlockY0_[ii]*m);
    THBlockZ0.push_back(THBlockZ0_[ii]*m);
    THBlockDxdz.push_back(THBlockDxdz_[ii]);
    THBlockDydz.push_back(THBlockDydz_[ii]);
    THBlockLength.push_back(THBlockLength_[ii]*m);
    THBlockHdx.push_back(THBlockHdx_[ii]*m);
    THBlockHdy.push_back(THBlockHdy_[ii]*m);
    THBlockGeantMaterial.push_back(THBlockGeantMaterial_[ii]);
    THBlockName.push_back(THBlockName_[ii]);
  }

  HABlockNblock = 7;
  //======================================================================= 
  // Hadron Absorber 
 
  G4double HABlockX0_[]    = { 1.651 , -1.3208, 0.0     , 0.3302  , 0.0   ,  0.0   , 0.0     };
  G4double HABlockY0_[]    = { 0.0   ,  0.0   , 0.0     , 0.0     , 0.0   , -1.651 , 1.651  };
  G4double HABlockZ0_[]    = { 723.48,  723.48, 726.4772, 728.7632, 723.48,  723.48, 723.48  };
  G4double HABlockDxdz_[]  = { 0.0   ,  0.0   , 0.0     , 0.0     , 0.0   ,  0.0   , 0.0     };
  G4double HABlockDydz_[]  = { 0.0   ,  0.0   , 0.0     , 0.0     , 0.0   ,  0.0   , 0.0     };
  G4double HABlockLength_[]= { 5.2832,  5.2832, 2.2860  , 0.9144  , 2.4384,  2.9972, 2.9972  };
  G4double HABlockHdx_[]   = { 0.9906,  0.6604, 0.6604  , 2.3114  , 0.6477,  0.6604, 0.6604  };
  G4double HABlockHdy_[]   = { 2.6416,  2.6416, 2.6416  , 2.6416  , 0.6477,  0.9906, 0.9906};
  G4int HABlockGeantMaterial_[] = {10,  10    , 10      , 17      , 9     ,  10    , 10    };
  G4String HABlockName_[]  = {"HABLK1" , "HABLK2" ,"HABLK3"   ,"HABLK4"   ,"HABLK5" ,"HABLK6"  ,"HABLK7"};
 
  for (G4int ii=0;ii<HABlockNblock;ii++){
    HABlockX0.push_back(HABlockX0_[ii]*m);
    HABlockY0.push_back(HABlockY0_[ii]*m);
    HABlockZ0.push_back(HABlockZ0_[ii]*m);
    HABlockDxdz.push_back(HABlockDxdz_[ii]);
    HABlockDydz.push_back(HABlockDydz_[ii]);
    HABlockLength.push_back(HABlockLength_[ii]*m);
    HABlockHdx.push_back(HABlockHdx_[ii]*m);
    HABlockHdy.push_back(HABlockHdy_[ii]*m);
    HABlockGeantMaterial.push_back(HABlockGeantMaterial_[ii]);
    HABlockName.push_back(HABlockName_[ii]);
  }
  
  
  HornCurrent=200.; //kA
  static const G4double in=2.54*cm;
  
  NPHorn2EndN=3;
  G4double PHorn2EndZ0_[]     ={135.861        ,137.611     ,139.486};
  G4double PHorn2EndLength_[] ={1.75           ,2.188       ,.625};
  G4double PHorn2EndRin_[]    ={12.719         ,12.532      ,11.};
  G4double PHorn2EndRout_[]   ={14.405         ,14.469      ,12.532};
  G4int PHorn2EndGeantMat_[]  ={31             ,  9         ,  9   };
  G4String PHorn2EndVolName_[]={"PHorn2InsRing","PHorn2CPB1","PHorn2CPB2"};

  for (G4int ii=0;ii<NPHorn2EndN;ii++){
    PHorn2EndZ0.push_back(PHorn2EndZ0_[ii]*in);
    PHorn2EndLength.push_back(PHorn2EndLength_[ii]*in);
    PHorn2EndRin.push_back(PHorn2EndRin_[ii]*in);
    PHorn2EndRout.push_back(PHorn2EndRout_[ii]*in);
    PHorn2EndGeantMat.push_back(PHorn2EndGeantMat_[ii]);
    PHorn2EndVolName.push_back(PHorn2EndVolName_[ii]);
  }

  NPHorn1EndN=3;
  G4double PHorn1EndZ0_[]     ={126.092        ,127.842     ,129.718};
  G4double PHorn1EndLength_[] ={1.75           ,2.188       ,.624};
  G4double PHorn1EndRin_[]    ={6.00           ,5.815       ,4.25};
  G4double PHorn1EndRout_[]   ={7.687          ,7.75        ,5.815};
  G4int PHorn1EndGeantMat_[]  ={31             ,  9         ,  9   };
  G4String PHorn1EndVolName_[]={"PHorn1InsRing","PHorn1CPB1","PHorn1CPB2"};

  for (G4int ii=0;ii<NPHorn1EndN;ii++){
    PHorn1EndZ0.push_back(PHorn1EndZ0_[ii]*in);
    PHorn1EndLength.push_back(PHorn1EndLength_[ii]*in);
    PHorn1EndRin.push_back(PHorn1EndRin_[ii]*in);
    PHorn1EndRout.push_back(PHorn1EndRout_[ii]*in);
    PHorn1EndGeantMat.push_back(PHorn1EndGeantMat_[ii]);
    PHorn1EndVolName.push_back(PHorn1EndVolName_[ii]);
  }

  //Spider supports for Horn1
  NHorn1SpiderSupportPlanesN=3;
  NHorn1SpidersPerPlaneN=3;
  Horn1SpiderSupportZ0.push_back(19.261*in);
  Horn1SpiderSupportZ0.push_back(42.862*in);
  Horn1SpiderSupportZ0.push_back(82.444*in);

  NumiHornSpiderSupport dummy=NumiHornSpiderSupport();
  dummy.stripW=.5*in;
  dummy.stripH=3.685*in;
  dummy.stripL=.076*in;
  dummy.topW=1.437*in;
  dummy.topH=.108*in;
  dummy.topL=.4*in;
  dummy.bottomW=1.437*in;
  dummy.bottomH=.542*in;
  dummy.bottomL=.4*in;
  dummy.bottomThickMid=.3*in;  //Thickness of bottom part in the middle
  dummy.bottomR=1.188*in;
  dummy.ceramicRodR=.491*in;
  
  Horn1SS.push_back(dummy);
  
  dummy.stripW=.5*in;
  dummy.stripH=3.123*in;
  dummy.stripL=.063*in;
  dummy.topW=1.437*in;
  dummy.topH=.108*in;
  dummy.topL=.4*in;
  dummy.bottomW=1.437*in;
  dummy.bottomH=.454*in;
  dummy.bottomL=.4*in;
  dummy.bottomThickMid=.3*in;  //Thickness of bottom part in the middle
  dummy.bottomR=1.75*in;
  dummy.ceramicRodR=.491*in;
  
  Horn1SS.push_back(dummy);
 
  dummy.stripW=.5*in;
  dummy.stripH=1.56*in;
  dummy.stripL=.031*in;
  dummy.topW=1.437*in;
  dummy.topH=.108*in;
  dummy.topL=.4*in;
  dummy.bottomW=1.437*in;
  dummy.bottomH=.379*in;
  dummy.bottomL=.4*in;
  dummy.bottomThickMid=.3*in;  //Thickness of bottom part in the middle
  dummy.bottomR=3.313*in;
  dummy.ceramicRodR=.491*in;
  
  Horn1SS.push_back(dummy);

  //Spider supports for Horn2
  NHorn2SpiderSupportPlanesN=1;
  NHorn2SpidersPerPlaneN=3;
  Horn2SpiderSupportZ0.push_back(59.809*in);

  dummy.stripW=.5*in;
  dummy.stripH=7.033*in;
  dummy.stripL=.125*in;
  dummy.topW=1.437*in;
  dummy.topH=.108*in;
  dummy.topL=.4*in;
  dummy.bottomW=1.437*in;
  dummy.bottomH=.345*in;
  dummy.bottomL=.4*in;
  dummy.bottomThickMid=.3*in;  //Thickness of bottom part in the middle
  dummy.bottomR=5.812*in;
  dummy.ceramicRodR=.491*in;
 
  Horn2SS.push_back(dummy);
  //---------------------------------------------------------
  
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
