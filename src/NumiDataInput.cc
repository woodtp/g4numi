
//----------------------------------------------------------------------
//
//
// $Id: NumiDataInput.cc,v 1.32.2.25 2018/03/01 03:37:05 kordosky Exp $
//----------------------------------------------------------------------

//C++
#include <string>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdio.h>

#include "NumiDataInput.hh"
#include "G4ThreeVector.hh"
#include "NumiHornSpiderSupport.hh"
#include "G4UserLimits.hh"
//#include "globals.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"

#include "CLHEP/Units/PhysicalConstants.h"

static const G4double in = 2.54*CLHEP::cm;
static const G4double fTargetZ0_ref   = -0.35*CLHEP::m; //LE000 position
static const G4double fHPBaffleZ0_ref = -2.53*CLHEP::m; //LE000 position
static const G4double fBeamZ0_ref     = -4.0*CLHEP::m; //this is arbitrary right now
//Note: The baffle LE000 position was -3.04m. This doesn't make sense. if the baffle
//is at -3.04 then the distance between the DS edge of the
//baffle and the US edge of the first Vertical target fin is 1.19m. The Documentation,
//http://www-numi.fnal.gov/numwork/tdh/TDH_V2_4.2.2-baffle.pdf, says it should be 68cm

NumiDataInput* NumiDataInput::fNumiDataInput = 0;

//---------------------------------------------------------------------------------
NumiDataInput* NumiDataInput::GetNumiDataInput()
{
  //G4cout << "Requesting NumiDataInput " << G4endl;
    if (!fNumiDataInput) {
        fNumiDataInput = new NumiDataInput();
    }
    return fNumiDataInput;
}

//---------------------------------------------------------------------------------
NumiDataInput::NumiDataInput()
   :debugOn(false),
    fPrintGeometry(false),
    fOkToRun(false),
    fUseNuBeam(false),
    fUseCorrHornCurrent(true),
    fUseDetailedProtonBeam(false),
    fUseWaterInTgt(false),
    fUseHornMisalign(false),
    fUseTgtDensity(false),

    fBeamConfig(""),
    fPlaylist(""),
    fTargetConfig(""),
    fIHornConfig(""),
    fHornConfig(""),
    fSimulation(""),
    fSubSimulation(""),
    fPhysicsListName("<unset>"),

    fRunPeriod(-999),
    fDebugLevel(0),

    fProton_outR(-999.0),
    fProton_inR(-999.0),
    fProtonDiv(-999.0),
    fProtonSpread(-999.0),
    fProton_cosx(-999.0),
    fProton_cosy(-999.0),

    fLengthOfWaterInTgt(4.0*CLHEP::cm),

    fPrintInfo(1)

{

   if(fPrintInfo > 0 || debugOn) G4cout << "NumiDataInput Constructor Called" << G4endl;

   if (fNumiDataInput)
   { G4Exception("NumiDataInput","NumiDataInput",FatalException,"NumiDataInput constructed twice"); }

//  fNumiDataInput = this;


   if(fPrintInfo > 0 || debugOn)
   {
      G4cout << G4endl;
      G4cout << G4endl;
      G4cout << "********************************************************************" << G4endl;
      G4cout << "********************************************************************" << G4endl;
      G4cout << "NumiDataInput - Initializing..." << G4endl;
      G4cout << "********************************************************************" << G4endl;

   }

   /*
  //
  //this is flugg stuff
  //
  bool useFile = false;
  int runPeriodFile = -999;
  double targetZFile = -999.;
  double hornCurrentFile = -999.;

  std::ifstream datafile("../runConfig.inp");
  if (datafile.is_open()) {
    useFile = true;
    datafile  >> runPeriodFile >> targetZFile >> hornCurrentFile;
    datafile.close();
    G4cout << G4endl << G4endl << G4endl << G4endl
    << G4endl << G4endl << G4endl << G4endl
    << G4endl << G4endl << G4endl << G4endl
    << G4endl << G4endl << G4endl << G4endl
        << "Using parameters from runConfig.inp: " << G4endl
           << "  Run Period = " << runPeriodFile << G4endl
           << "  Target Z = " << targetZFile << G4endl
           << "  Horn Current = " << hornCurrentFile << G4endl;
    if (hornCurrentFile == 185) hornCurrentFile = 182.1;
    else if (hornCurrentFile == 170) hornCurrentFile = 167.3;
    else if (hornCurrentFile == 200) hornCurrentFile = 196.9;
    else if (hornCurrentFile == 0) hornCurrentFile = 0;
    else {
        G4cout << "Unrecognized horn current " << hornCurrentFile << ", bailing." << G4endl;
        assert(false);
    }
  }
  //
  */



  NImpWeightOn = true;
  createNuNtuple=false;  createHadmmNtuple=false;
  createASCII=false;     createBXDRAW = false;
  useFlukaInput = false; useMarsInput=false;
  createTarNtuple=false;

  //
  //can set useMuonInupt and useMuonBeam from macro
  //
  useMuonInput = false;
  NInputParts = 0;
  NInputPart  = 0;
  reWeightDeltas = false;
  nSplitDeltas = 0;
  useMuonBeam = false;
  hadmmNtupleDir  = "./";
  hadmmNtupleName = "";
  fNEvents = -1;
  solidMuMons = false;
  absorberConfig = "None";
  fHorn1IsAlternate = false;
  fHorn1IsRefined = false;
  simAbsBkg = false;
  createAbsBkgNtuple = false;
  absbkgNtupleName = "";
  absbkgNtupleDir = "./";
  simDRays = false;
  useZPosCut = false;
  muonBeamMomentum = 10.0*CLHEP::GeV;
  muonBeamShape = "";
  muonBeamZPos = -1.0*CLHEP::mm;
  muonBeamGaussXsig = 0.0*CLHEP::mm;
  muonBeamGaussYsig = 0.0*CLHEP::mm;

  useTestBeam = false;
  useDecayPipeSelect = false;
  KillTracking = true; // false for ahimmel
  testTheta = M_PI/6.;

   StepLimit = 0.0;
   TimeLimit = 60.0*CLHEP::s;

  extNtupleFileName=""; //fluka or mars or muon ntuple with particles coming of the target
  //Set the energy threshold for 'killing' particles
   KillTrackingThreshold = 0.05*CLHEP::GeV; //for defaut neutrino MC g4numi
   //KillTrackingThreshold = 0.001*CLHEP::GeV; //for muon beam MC


   //base name for output files:
   nuNtupleName    = "nuNtuple";
   tarNtupleName    = "tarNtuple";
   //
   //can set hadmmNtupleDir and hadmmNtupleName from macro
   //
   //hadmmNtupleDir  = "./";
   //hadmmNtupleName = "";
   asciiName       = "asciiOut";
   bxdrawName      = "bxdrawOut";
   RunNumber       = "0000";
   geometry        = "";

   //
   // Denotes the change in sigma to the
   // rock density and muon alcove wall location
   //
   materialSigma   = 0;


   //===================================
   //------ Use Macro---------------
   // useMacro = true if you actually want to change beam parameters
   // from the macro

   useMacro = false;

  //====================================
  //-- Changes for "Air Horns"---------- >> G4numi validation
  //-- Changes for "Truncated Horns"-----
  airhrn =false; // airhrn must be changed before compilation
  vacuumworld=false;
  jCompare = false; // make horns have the same B field;
  g3Chase = true;

if(!vacuumworld && !airhrn){
  hrnmat = 9;   // Al
  hrnmatcr =31; // CT852
  hallmat=15;   // Air
}
 else {
   if(airhrn){
     hrnmat=15;
     hrnmatcr=15;
     hallmat=15;
    }
   if(vacuumworld){
     hrnmat=16;
     hrnmatcr=16;
     hallmat=16;
   }
 }

//======================================

//======================================
//--Ray Tracing parameter-------------
//--------------------------------------
 raytracing =false;
 //==============================
 //------------------------------

 //initialize zpoints
 NZpoint=48;
 Zpoint.push_back(0*CLHEP::m); //1
 Zpoint.push_back(0.5*CLHEP::m); //2
 Zpoint.push_back(0.7*CLHEP::m); //3
 Zpoint.push_back(0.75*CLHEP::m); //4
 Zpoint.push_back(0.8*CLHEP::m); //5
 Zpoint.push_back(0.85*CLHEP::m); //6
 Zpoint.push_back(0.9*CLHEP::m); //7
 Zpoint.push_back(1*CLHEP::m); //8
 Zpoint.push_back(1.25*CLHEP::m); //9
 Zpoint.push_back(1.5*CLHEP::m); //10
 Zpoint.push_back(1.75*CLHEP::m); //11
 Zpoint.push_back(2*CLHEP::m);//12
 Zpoint.push_back(2.25*CLHEP::m);//13
 Zpoint.push_back(2.5*CLHEP::m);//14
 Zpoint.push_back(2.75*CLHEP::m); //15
 Zpoint.push_back(3*CLHEP::m);//16
 Zpoint.push_back(3.25*CLHEP::m);//17
 Zpoint.push_back(3.5*CLHEP::m);//18
 Zpoint.push_back(4*CLHEP::m);//19
 Zpoint.push_back(6*CLHEP::m);//20
 Zpoint.push_back(9*CLHEP::m);//21
 Zpoint.push_back(9.5*CLHEP::m); //22
 Zpoint.push_back(9.25*CLHEP::m); //23
 Zpoint.push_back(10*CLHEP::m);//24
 Zpoint.push_back(10.25*CLHEP::m);//25
 Zpoint.push_back(10.5*CLHEP::m);//26
 Zpoint.push_back(10.75*CLHEP::m);//27
 Zpoint.push_back(11.80*CLHEP::m);//28
 Zpoint.push_back(11.9*CLHEP::m);//29
 Zpoint.push_back(11*CLHEP::m);//30
 Zpoint.push_back(11.1*CLHEP::m);//31
 Zpoint.push_back(11.2*CLHEP::m);//32
 Zpoint.push_back(11.25*CLHEP::m);//33
 Zpoint.push_back(11.5*CLHEP::m);//34
 Zpoint.push_back(11.75*CLHEP::m);//35
 Zpoint.push_back(12*CLHEP::m);//36
 Zpoint.push_back(12.25*CLHEP::m);//37
 Zpoint.push_back(12.5*CLHEP::m);//38
 Zpoint.push_back(12.75*CLHEP::m);//39
 Zpoint.push_back(13*CLHEP::m);//40
 Zpoint.push_back(13.25*CLHEP::m);//41
 Zpoint.push_back(13.5*CLHEP::m);//42
 Zpoint.push_back(13.75*CLHEP::m);//43
 Zpoint.push_back(14*CLHEP::m);//43
 Zpoint.push_back(14.5*CLHEP::m);//44
 Zpoint.push_back(15.*CLHEP::m);//45
 Zpoint.push_back(16.*CLHEP::m);//46
 Zpoint.push_back(20*CLHEP::m);//47
 Zpoint.push_back(45.7*CLHEP::m);//48
 createZpNtuple=false;
 zpNtupleName="zpNtuple";

/*
  This is flugg stuff
//==============================================================
  //Target Configuration (LE010 = -10.0, LE100 = -100.0, etc.)
  // runPeriod corresponds to Runs I, II, III, IV
  //=======================================================================
  double TargetConfigZ = -10.0*CLHEP::cm;
  int runPeriod = 0;

  if (useFile) {
      runPeriod = runPeriodFile;
      TargetConfigZ = -1*targetZFile*CLHEP::cm;
  }
*/

  G4float beam_x_dir = 0;
  G4float beam_y_dir = 0;
  G4float beam_z_dir = 1;//cos(.01*pi/180);
  //actual dm is 5/13e-4 radians
  G4float beam_x_pos = 0;
  G4float beam_y_pos = 0;
  G4float beam_z_pos = fBeamZ0_ref;

/*
  G4float beam_z_pos = -4.0*CLHEP::m;
  This is wrong somehow. If running G4NuMI with a proton beam
  need the beam to start before the baffle to get the effect of protons
  interacting with the baffle. - Laura

//G4float beam_z_pos = -4.0*CLHEP::m+ TargetConfigZ;
  // the reason for the beam_z_pos change was to move the beam to start
  // immediately before the target so that the beam spot interaction point
  // would remain constant, but the angle would change.
  */

  protonMomentum = 120.*CLHEP::GeV;
  //  beamSigmaY     = 1.1*CLHEP::mm;//1.25*CLHEP::mm;
  //  beamSigmaX     = 1.1*CLHEP::mm;//1.1*CLHEP::mm;
  // MAK 10/14/16 - hack in a wider ME beam
  beamSigmaY     = 1.4*CLHEP::mm;//1.25*CLHEP::mm;
  beamSigmaX     = 1.4*CLHEP::mm;//1.1*CLHEP::mm;

  beamDirection  = G4ThreeVector(beam_x_dir,beam_y_dir,beam_z_dir);
  beamPosition  = G4ThreeVector(beam_x_pos,beam_y_pos,beam_z_pos);

  protonKineticEnergy = sqrt(pow((.938*CLHEP::GeV),2)+pow(protonMomentum,2))-0.938*CLHEP::GeV;

    //Rock        1
  //=======================================================================
  RockRadius  = 10.0*CLHEP::m;
  RockHalfLen = 1200.0*CLHEP::m;
  RockDensity = 2.41*CLHEP::g/CLHEP::cm3; // not
  RockRadLen  = 0.0;        // used




  constructTarget = true;
  //TargetArea          1
  //=======================================================================
  TargetAreaZ0       = -6.7*CLHEP::m;  //was -4.0*m (08/09/05);
  TargetAreaLength   = 52.398*CLHEP::m;//was 49.28*m (08/09/05);

  // TargetAreaHeight and TargetAreaWidth were 6.0 meters in Zarko's older version,
  // but I had to extend them by 1 m to fit the concrete chase in TGAR
  // (working from dimensions and placement given in the Numi Technical Design Handbook)
  // I (Zarko) had to add 1.5m more to fit all blocks
  TargetAreaHeight   = 8.5*CLHEP::m;
  TargetAreaWidth    = 8.5*CLHEP::m;
  TargetAreaGEANTmat = 15;

  // Target   1
  //=======================================================================
  TargetX0           = 0.0;
  TargetY0           = 0.0;
//  TargetZ0           = -0.35*m + TargetConfigZ;
  TargetZ0           = fTargetZ0_ref; // this is LE000 position
  TargetDxdz         = 0.0; // doesn't
  TargetDydz         = 0.0; // work properly yet
  TargetSLength      = 20.*CLHEP::mm;
  TargetSWidth       = 6.40E-03*CLHEP::m;
  TargetSHeight      = 18.0E-03*CLHEP::m;
  TargetCPGRadius    = 3.2*CLHEP::mm; // Cooling pipe groove
  TargetCPGPosition  = 10.7*CLHEP::mm;
  TargetEndRounded   = true;
  TargetSegmentNo    = 47;
  TargetSegmentPitch = 0.3*CLHEP::mm;
  TargetA            = 12.01*CLHEP::g/CLHEP::mole;
  TargetZ            = 6.;
  TargetDensity      = 1.78*CLHEP::g/CLHEP::cm3; //1.815*CLHEP::g/CLHEP::cm3;//1.754*CLHEP::g/CLHEP::cm3;
  TargetRL           = 25.692;
  TargetGEANTmat     = 18;

  //Budal Monitor
  //this is with respect to the target -Laura
  //
  BudalX0 = 0.0;
  BudalY0 = 0.0;
  BudalZ0 = -16.72*CLHEP::cm;
  BudalDxdz = 0.0;
  BudalDydz = 0.0;

  /*

  this is flugg stuff. Note if the target moves 1.1cm so does the baffle!

  if (runPeriod == 1) {
    TargetY0 = -1.1*mm;
    BudalY0  = 2.26*mm;
  }
  else {
    // Align Target and Budal
    TargetY0 = 0.0;
    BudalY0 = 0.0;
  }

  if (runPeriod == 2 || runPeriod == 3) {
    // Change LE010 to LE009
    if (TargetConfigZ == -10.0*CLHEP::cm) {
        TargetZ0 += 1.1*CLHEP::cm;
    }
  }
  */

  //HPBaffle           1 // Only the length and position can be changed currently
  //=======================================================================
  HPBaffleGEANTMat   =  18;
  HPBaffleX0         =  0.00;
  HPBaffleY0         =  0.00;
//HPBaffleZ0         = -3.04*m + TargetConfigZ; //for flugg
  //This doesn't make sense. if the baffle is at -3.04 then the distance between the DS edge of the
  //baffle and the US edge of the first Vertical target fin is 1.19m. The Documentation,
  //http://www-numi.fnal.gov/numwork/tdh/TDH_V2_4.2.2-baffle.pdf, says it should be 68cm
  HPBaffleZ0         = fHPBaffleZ0_ref; //This is LE000 position
  HPBaffleDXDZ       =  0.0;
  HPBaffleDYDZ       =  0.0;
//  HPBaffleLength     =  1.20*m;
  HPBaffleLength     =  1.50*CLHEP::m; //Why was this 1.2m?!
  HPBaffleRin        =  5.5*CLHEP::mm;
  HPBaffleRout       =  3.*CLHEP::cm;

  //Cooling pipes
  NCPipeN = 19;
  //=======================================================================
  // CPipeDXDZ=-99999 and CPipeDYDZ=0 puts curved part in z-y plane
  //=======================================================================
  G4int CPGeantMat_[]        = {10     ,   10   ,10     , 10     , 10    , 10    ,  10   , 10    , 10     , 10           ,10            , 31     , 31      ,10        ,10        ,  10     ,  10     , 10      ,10  };
  G4bool CPipeFilledWater_[] = { true  , true   ,true   , true   ,  true , true  ,false  ,true   , false  , true         ,true          , true   ,true     ,true      ,true      , true    , true    , true    , true};
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
    if(airhrn){
      CPGeantMat.push_back(15);
    }
    else{
    CPGeantMat.push_back(CPGeantMat_[ii]);
    }
    CPipeFilledWater.push_back(CPipeFilledWater_[ii]);
    CPipeX0.push_back(CPipeX0_[ii]*CLHEP::m);
    CPipeY0.push_back(CPipeY0_[ii]*CLHEP::m);
    CPipeZ0.push_back(CPipeZ0_[ii]*CLHEP::m);
    CPipeDXDZ.push_back(CPipeDXDZ_[ii]);
    CPipeDYDZ.push_back(CPipeDYDZ_[ii]);
    CPipeLength.push_back(CPipeLength_[ii]*CLHEP::m);
    CPipeRadiusOut.push_back(CPipeRadiusOut_[ii]*CLHEP::m);
    CPipeRadiusIn.push_back(CPipeRadiusIn_[ii]*CLHEP::m);
    CPipeWallThick.push_back(CPipeWallThick_[ii]*CLHEP::m);
    CPipeCurvRad.push_back(CPipeCurvRad_[ii]*CLHEP::m);
    CPipeOpenAng.push_back(CPipeOpenAng_[ii]*CLHEP::deg);
    CPipeCloseAng.push_back(CPipeCloseAng_[ii]*CLHEP::deg);
    CPipeVolName.push_back(CPipeVolName_[ii]);
  }

  //Container
  //=======================================================================
  // Z0 with respect to the first target fin
  NContainerN=21;

  G4double CTubeZ0_[]     ={-.42181, -.421556, -.41555, -.35855 , -0.3522 ,-.3332 ,-0.0802, -0.123  ,-0.109  ,-0.0572 ,-0.0492 , -0.04  ,-0.042 ,-.022 , -0.013   , -0.0035  ,0.011     , 0.011       ,  0.011      , 0.9785 , 0.9896  };
  G4double CTubeLength_[] ={.254e-3, 6.006e-3, 5.7e-2 , 6.35e-3 , 19e-3   , 253e-3, 23e-3 , 11e-2   , 66e-3  ,8e-3   , 17e-3  , 6e-3   ,0.02   ,15e-3  , 24e-3    , 14.5e-3  , 0.9785   ,1e-6         , 0.9786      ,5e-4    , 1e-6};
  G4double CTubeRin_[]    ={  0.   , 1.27e-2 , 1.7e-2 , 22.22e-3, 17.5e-3 , 77e-3 , 77e-3 , 14.6e-3 , 16.e-3 , 74e-3 , 24e-3  , 18.5e-3,16e-3  ,16e-3  , 14.6e-3  , 16e-3    , 14.6e-3  ,15.101e-3    , 15.1e-3     , 0.0    , 0.};
  G4double CTubeRout_[]   ={22e-3  , 34.67e-3, 1.9e-2 , 34.67e-3, 107e-3  , 83e-3 , 120e-3, 16.0e-3 , 21.4e-3,106e-3 , 106e-3 , 24e-3  ,18.5e-3,18.5e-3, 16e-3    , 16.4e-3  , 15e-3    , .35         , 15.101e-3   ,14.15e-3, 15.1e-3  };
  G4int CTubeGeantMat_[]  ={   5   ,  10     ,   10   , 10      , 10      , 10    , 10    ,  31     , 10     , 10    , 10     ,  10    , 10    ,10     , 10       , 9        , 9        ,15           ,  15         ,  5     ,  15};
  G4String CTubeVolName_[]={"BeUp1", "BeUp2" , "Added", "BeUp3" , "BFront", "Body", "BEnd","CerTube", "Conn1","CLid1", "CLid2", "Conn2","Conn3","Tube1a" ,"Tube1b", "AlTube1", "AlTube2","TGTExitCyl1","TGTExitCyl2", "BeDW" ,"TGTExitTop"};

 TargetVirtualCanisterHeight = 63.0*CLHEP::mm + 0.075*CLHEP::mm; // Sett Mova Target Heigth...
 TargetVirtualCanisterWidth = 13.75*CLHEP::mm; // a bit bigger than LE as well,

 for (G4int ii=0;ii<NContainerN;ii++){
    CTubeZ0.push_back(CTubeZ0_[ii]*CLHEP::m);
    CTubeLength.push_back(CTubeLength_[ii]*CLHEP::m);
    CTubeRin.push_back(CTubeRin_[ii]*CLHEP::m);
    CTubeRout.push_back(CTubeRout_[ii]*CLHEP::m);
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
    TgtRingZ0.push_back(TgtRingZ0_[ii]*CLHEP::m);
    TgtRingLength.push_back(TgtRingLength_[ii]*CLHEP::m);
    TgtRingRin.push_back(TgtRingRin_[ii]*CLHEP::m);
    TgtRingRout.push_back(TgtRingRout_[ii]*CLHEP::m);
    TgtRingGeantMaterial.push_back(TgtRingGeantMaterial_[ii]);
    TgtRingVolName.push_back(TgtRingVolName_[ii]);
 }

  //Tunnel         1
  //=======================================================================
  TunnelZ0       = 45.6985*CLHEP::m; //was 45.28*m (08/09/05);
  TunnelRadius   = 3.3*CLHEP::m+.5*CLHEP::m; //added .5m because hadron absorber does not fit entirely inside 3.3
  TunnelLength   = 693.4415*CLHEP::m; //was 693.86*m (08/09/05);
  TunnelA        = 0.0;
  TunnelZ        = 0.0;
  TunnelGEANTmat = 15;
  BeamAngle      = 0.05835; // .05835 in radians.

  //ShieldNshield  5
  //=======================================================================
  ShieldX0       = 0.0;
  ShieldY0       = 0.0;
  ShieldZ0       = 45.699*CLHEP::m; //was 45.28*m (08/09/05);
  ShieldDxdz     = 0.0; // not
  ShieldDydz     = 0.0; // used
  ShieldLength   = 676.681*CLHEP::m; //was 677.1*m (08/09/05);
  ShieldRout     = 2.23*CLHEP::m;
  ShieldRin      = 1.0097*CLHEP::m;
  ShieldGEANTmat = 17;

  //DecayPipe          1
  //=======================================================================
  DecayPipeZ0        = 45.699*CLHEP::m; //was 45.28*CLHEP::m (08/09/05);
  //DecayPipeRadius    = 0.9716*CLHEP::m; // was 0.9906 but doesnt correlate w gnumi
  DecayPipeRadius    = 0.9906*CLHEP::m;
  DecayPipeLength    = 676.681*CLHEP::m; //was 677.1*CLHEP::m (08/09/05);
  DecayPipeFWinThick = 1.60E-3*CLHEP::m;
  DecayPipeEWinThick = 4.76E-3*CLHEP::m;
  DecayPipeWallThick = 1.905E-2*CLHEP::m;
  DecayPipeA         = 55.85;
  DecayPipeZ         = 26.0;
  DecayPipeGEANTmat  = 10;
  DecayPipeFWinmat   = 9;
  DecayPipeEWinmat   = 10;
  HeInDecayPipe      = true;
  applyDecayPipeMagneticField = false;
  // New Target Hall by Zach Barnett

  //==========================================================================
  // TargetHallChase       1

  NTHConcreteSectionsN = 6;

  // units in meters, dimensions based on Numi note.
  G4double THConcreteX0_[]= {3.36073, 2.60670, 0.0, -2.60670, -3.36073, 0.0};
  G4double THConcreteY0_[]= {2.801, 0.0, -2.801, 0, 2.801, 3.45665};
  // these Z0 values are set with respect to the ROCK volume, but NumiTargetHall.cc adjusts them to the TGAR volume
  G4double THConcreteZ0_[]= {19.499, 19.499, 19.499, 19.499, 19.499, 19.499};
  G4double THConcreteLength_[]= {52.397, 52.397, 52.397, 52.397, 52.397, 52.397};
  G4double THConcreteHdx_[]={.754025,.46675,3.07345,.46675,.754025, 2.921};
  G4double THConcreteHdy_[]={.42705,2.37395,.42705,2.37395,.42705, .2286};
  G4double THConcreteDxdz_[]={0,0,0,0,0, 0.0};
  G4double THConcreteDydz_[]={0,0,0,0,0, 0.0};
  G4int THConcreteGeantMaterial_[] = {17,17,17,17,17, 17};
  G4String THConcreteName_[]={"Section1","Section2","Section3","Section4","Section5","Section6(lid)"};

  for (G4int ii=0;ii<NTHConcreteSectionsN;ii++){
	  THConcreteX0.push_back(THConcreteX0_[ii]*CLHEP::m);
	  THConcreteY0.push_back(THConcreteY0_[ii]*CLHEP::m);
	  THConcreteZ0.push_back(THConcreteZ0_[ii]*CLHEP::m);
	  THConcreteDxdz.push_back(THConcreteDxdz_[ii]*CLHEP::m);
	  THConcreteDydz.push_back(THConcreteDydz_[ii]*CLHEP::m);
	  THConcreteLength.push_back(THConcreteLength_[ii]*CLHEP::m);
	  THConcreteHdx.push_back(THConcreteHdx_[ii]*CLHEP::m);
	  THConcreteHdy.push_back(THConcreteHdy_[ii]*CLHEP::m);
	  THConcreteGeantMaterial.push_back(THConcreteGeantMaterial_[ii]);
	  THConcreteName.push_back(THConcreteName_[ii]);
  }

  //==========================================================================

  //=======================================================================
  // Target Hall shielding (18 blocks, numbered 0-17)

  if(g3Chase){
    THBlockNblock = 24;
  }
  else{
    THBlockNblock=18;
  }
  // reminder: all length dimensions for the target shielding blocks are in meters.
  //These are the Duratek blocks mentioned in the Numi Technical Design Handbook.

  //Blocks 15-17 are not actually in the Numi Technical Design Handbook--they are used
  //to provide a covering over the top

  G4double THBlockX0_[]    = { 0.0,      //Block 0
							   1.3554,   //Block 1
							   -1.3554,  //Block 2
							   1.7158,    //Block 3
							   .6904,     //Block 4
							   -.6904,    //Block 5
							   -1.7158,   //Block 6
							   1.0254,      //Block 7
							   -1.0254,     //Block 8
							   1.7158,    //Block 9
							   -1.7158,   //Block 10
							   1.0254,      //Block 11
							   -1.0254,     //Block 12
							   1.7158,    //Block 13
							   -1.7158,    //Block 14
							   0.6904,        //Block15 (top block1--one of the slightly enlarged blocks)
							   -0.6904,        //Block16 (top block2--one of the slightly enlarged blocks)
							   0.0,        //Block17 (top block3)
			       0.0,        // Block 18 (bottom gnumi block)
			       0.6373,  // Block 19 ( side 1 gnumi block)
			       -0.6373, // Block 20 (Side 2 gnumi block)
			       0.0, // Block 21 (top gnumi block bf Horn2)
			       0.0,//
			       0.0//
  };

  G4double THBlockY0_[]    = { -2.03895,      //Block 0
							   -2.03895,      //Block 1
							   -2.03895,      //Block 2
							   -1.03895,      //Block 3
							   -1.36895,      //Block 4
							   -1.36895,      //Block 5
							   -1.03895,      //Block 6
							   -0.36895,     //Block 7
							   -0.36895,     //Block 8
							   0.29105,     //Block 9
							   0.29105,     //Block 10
							   0.96105,      //Block 11
							   0.96105,      //Block 12
							   1.62105,     //Block 13
							   1.62105,      //Block 14
							   1.98395,        //Block15 (top block1--one of the slightly enlarged blocks)
							   1.98395,        //Block16 (top block2--one of the slightly enlarged blocks)
							   1.29355,        //Block17 (top block3)
			       -0.993215, //Block 18 (bottom gnumi block)
                               -0.0377, //Block 19 (side 1 gnumi block)
                               -0.0377, //Block 20 (side 2 gnumi block)
                               0.657075,//Block 21 (top gnumi block)
			       0.657075, // Block 22 (top gnumi block after h2)
                               0.672075 //Block 23 above Horn2
  };



  // these Z0 values are set with respect to the ROCK volume, but NumiTargetHall.cc adjusts them to the TGAR volume

  G4double THBlockZ0_[]= {19.499};
  G4double THBlockDxdz_[] = {0.0};
  G4double THBlockDydz_[] = {0.0};

  /*
  // these Z0 values are set with respect to the ROCK volume, but NumiTargetHall.cc adjusts them to the TGAR volume
  G4double THBlockZ0_[]    = { 20.64,        //Block 0
  20.64,        //Block 1
							   20.64,        //Block 2
							   20.64,        //Block 3
							   20.64,        //Block 4
							   20.64,        //Block 5
							   20.64,        //Block 6
							   20.64,        //Block 7
							   20.64,        //Block 8
							   20.64,        //Block 9
							   20.64,        //Block 10
							   20.64,        //Block 11
							   20.64,        //Block 12
							   20.64,        //Block 13
							   20.64,         //Block 14
							   20.64,        //Block15 (top block1--one of the slightly enlarged blocks)
							   20.64,        //Block16 (top block2--one of the slightly enlarged blocks)
							   20.64,        //Block17 (top block3)
							   19.499,//Block 18 bottom
							   19.499,//Block 19 side
							   19.499,//Block 20 side
							   19.499,//Block 21 top
							   19.499,//Block 22 top
							   19.499,//Block 23 top

  };


  G4double THBlockDxdz_[]  = { 0.0,        //Block 0
							   0.0,        //Block 1
							   0.0,        //Block 2
							   0.0,        //Block 3
							   0.0,        //Block 4
							   0.0,        //Block 5
							   0.0,        //Block 6
							   0.0,        //Block 7
							   0.0,        //Block 8
							   0.0,        //Block 9
							   0.0,        //Block 10
							   0.0,        //Block 11
							   0.0,        //Block 12
							   0.0,        //Block 13
							   0.0,         //Block 14
							   0.0,        //Block15 (top block1--one of the slightly enlarged blocks)
							   0.0,        //Block16 (top block2--one of the slightly enlarged blocks)
							   0.0,        //Block17 (top block3)
                             };

  G4double THBlockDydz_[]  = { 0.0,        //Block 0
			       0.0,        //Block 1
			       0.0,        //Block 2
			       0.0,        //Block 3
			       0.0,        //Block 4
			       0.0,        //Block 5
			       0.0,        //Block 6
			       0.0,        //Block 7
			       0.0,        //Block 8
			       0.0,        //Block 9
			       0.0,        //Block 10
			       0.0,        //Block 11
			       0.0,        //Block 12
			       0.0,        //Block 13
			       0.0,         //Block 14
			       0.0,        //Block15 (top block1--one of the slightly enlarged blocks)
			       0.0,        //Block16 (top block2--one of the slightly enlarged blocks)
			       0.0,        //Block17 (top block3)
                             };
  */

  // all normal blocks have the same Length, Hdx, and Hdy values.  These values for the slightly larger blocks are set
  // in NumiTargetHall.cc, line 103
  G4double THBlockLength_[]= {52.397};

  G4double THBlockHdx_[]   = {1.33/2.0};

  G4double THBlockHdy_[]   = {0.67/2.0};

  G4int THBlockGeantMaterial_[] = {10};

  G4String THBlockName_[]  = {"BLK0",
			      "BLK1",
			      "BLK2",
			      "BLK3",
			      "BLK4",
			      "BLK5",
			      "BLK6",
			      "BLK7",
			      "BLK8",
			      "BLK9",
			      "BLK10",
			      "BLK11",
			      "BLK12",
			      "BLK13",
			      "BLK14",
			      "BLK15",
			      "BLK16",
			      "BLK17",
			      "BLK18",
			      "BLK19",
			      "BLK20",
			      "BLK21",
			      "BLK22",
			      "BLK23"

                              };

  //This next block puts all the block coordinates into vectors and assigns the unit "meters" to the values.

  for (G4int ii=0;ii<THBlockNblock;ii++){
    THBlockX0.push_back(THBlockX0_[ii]*CLHEP::m);
    THBlockY0.push_back(THBlockY0_[ii]*CLHEP::m);
    if(ii==21){
      THBlockZ0.push_back(1.15*CLHEP::m);
    }
    else if(ii==22) THBlockZ0.push_back(30.34*CLHEP::m);
    else if(ii==23) THBlockZ0.push_back(12*CLHEP::m);
    else    THBlockZ0.push_back(THBlockZ0_[0]*CLHEP::m);
    THBlockDxdz.push_back(THBlockDxdz_[0]*CLHEP::m);
    THBlockDydz.push_back(THBlockDydz_[0]*CLHEP::m);
    THBlockLength.push_back(THBlockLength_[0]*CLHEP::m);
    THBlockHdx.push_back(THBlockHdx_[0]*CLHEP::m);
    THBlockHdy.push_back(THBlockHdy_[0]*CLHEP::m);
    THBlockGeantMaterial.push_back(THBlockGeantMaterial_[0]);
    THBlockName.push_back(THBlockName_[ii]);
  }

  //=======================================================================
  // Hadron Box Dimensions
  // 55.75' is the surveyed distance between the
  // downstream wall of Muon Alcove 1 and the
  // nominal center of the Decay Pipe.
  HadrBox_width = 324*.0254*CLHEP::m;
  HadrBox_height = 6.6294*CLHEP::m;
  HadrBox_length = 55.75*12*in-(4*12*in+2.24*in+9*12*in+9*in);

  /*
    //
    //this stuff is for FLUGG
    //
    //it is horrible!
    //

  HornCurrent=182100.*CLHEP::ampere;
  if (useFile) {
    HornCurrent = hornCurrentFile*CLHEP::ampere;
    if (HornCurrent < 250*CLHEP::ampere) HornCurrent *= 1000.; // Convert to kA
  }

  if (runPeriod == 4) {
      HornCurrent *= -1;
  }


  G4cout << "Running with: " << G4endl
       << "  Run Period = " << runPeriod << G4endl
       << "  Target Z = " << TargetConfigZ/CLHEP::cm << " cm" << G4endl
       << "  Horn Current = " << HornCurrent/CLHEP::ampere/1000. << " kA" << G4endl;
  */


  //
  //Horns
  //
  Horn1X0 = 0.0*CLHEP::cm;
  Horn1Y0 = 0.0*CLHEP::cm;
  Horn1Z0 = 3.0*CLHEP::cm;

  Horn2X0 = 0.0*CLHEP::cm;
  Horn2Y0 = 0.0*CLHEP::cm;
  Horn2Z0 = 10.0*CLHEP::m;

  Horn1Phi   = 0.0*CLHEP::rad;
  Horn1Theta = 0.0*CLHEP::rad;
  Horn1Psi   = 0.0*CLHEP::rad;

  Horn2Phi   = 0.0*CLHEP::rad;
  Horn2Theta = 0.0*CLHEP::rad;
  Horn2Psi   = 0.0*CLHEP::rad;

//  fHornWaterLayerThick = 0.0*CLHEP::mm; // July 14 2014, P.L.  Default is zero, should be 0.5 mm
  fHornWaterLayerThick = 1.0*CLHEP::mm; // Oct 18 2014, Based on studies done
//  mid-October, after fixing bug in the Magnetic field pointer assignments,
//  default is now what it should be
  fHorn1ExtraLayerAlum = 0.0*CLHEP::mm;
  fDumpBFieldPlease = false;
//  fDumpBFieldPlease = true;// To get enhanced field.. Dirty back door...

  NPHorn2EndN=3;

  G4double PHorn2EndZ0_[]     ={135.861        ,137.611     ,139.486};
  if (jCompare) {
	PHorn2EndZ0_[0] = 118.11;
	PHorn2EndZ0_[1] = 119.86;
	PHorn2EndZ0_[2] = 122.048;
  }
  G4double PHorn2EndLength_[] ={1.75           ,2.188       ,.625};
  G4double PHorn2EndRin_[]    ={12.719         ,12.532      ,11.};
  G4double PHorn2EndRout_[]   ={14.405         ,14.469      ,12.532};
  G4int PHorn2EndGeantMat_[]  ={31             ,  9         ,  9 };
  G4String PHorn2EndVolName_[]={"PHorn2InsRing","PHorn2CPB1","PHorn2CPB2"};

  for (G4int ii=0;ii<NPHorn2EndN;ii++){
    PHorn2EndZ0.push_back(PHorn2EndZ0_[ii]*in);
    PHorn2EndLength.push_back(PHorn2EndLength_[ii]*in);
    PHorn2EndRin.push_back(PHorn2EndRin_[ii]*in);
    PHorn2EndRout.push_back(PHorn2EndRout_[ii]*in);
    if(airhrn){
      PHorn2EndGeantMat.push_back(15);
    }
    else{
    PHorn2EndGeantMat.push_back(PHorn2EndGeantMat_[ii]);
    }
    PHorn2EndVolName.push_back(PHorn2EndVolName_[ii]);
  }

  NPHorn1EndN=3;
  G4double PHorn1EndZ0_[]     ={126.092        ,127.842     ,129.718};
  if (jCompare) {
	PHorn1EndZ0_[0] = 118.11;
	PHorn1EndZ0_[1] = 119.86;
	PHorn1EndZ0_[2] = 122.048;
  }
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
    if(airhrn){
      PHorn1EndGeantMat.push_back(15);
    }
    else{
      PHorn1EndGeantMat.push_back(PHorn1EndGeantMat_[ii]);
    }
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

  //Notes about locations (Leo Aliaga, Oct 2014)
  //Near & Far Detector locations:

  //Miniboone, Microboone and Sciboone (Zarko Pavlovic private e-mail 2014-09-12)

  //Minos Near and Far: positions from R. Hatcher using GENIE GNuMIFlux R-2_8_4 transformation.

  //Nova Near and Far: positions from R. Hatcher (using GENIE GNuMIFlux R-2_8_4 transformation) and updated by
  //Alex Radovic and Ryan Nichol.

  //Minerva values are in the center of the tracker: (0,0,716.95cm) in the Minerva coordinate system.
  //We provide these numbers to Robert Hatcher and he calculated the numbers below using GENIE GNuMIFlux R-2_8_4 transformation.
  //See also M. Kordosky Minerva docdb 10414.

  nNear=7;
  nFar=2;
  //We are checking now Minerva positions.
  G4double xdetNear[]    = {-0.5628   ,     0         ,   11.7174545 ,  -0.29,    26.04      ,    53.    ,    197.6};
  G4double ydetNear[]    = {-0.5329317,     0         ,   -3.3151325 ,  92.21,    78.64      ,    76.    ,     53.4};
  G4double zdetNear[]    = {1032.319  ,  1036.488     ,  992.9347347 , 841.76,   744.87      ,   679.    ,    339.4};
  G4String detNameNear[] = {"MINERvA" ,"MINOS NearDet","NOvA NearDet", "NDOS", "MiniBooNE" , "MicroBooNE", "SciBooNE"};

  G4double xdetFar[]     = {0             ,   11037.29645};
  G4double ydetFar[]     = {0             ,   -4162.557681};
  G4double zdetFar[]     = {735340        ,  810422.3274};
  G4String detNameFar[]  = {"MINOS FarDet", "Nova FarDet"};

  for(G4int ii=0;ii<nNear;ii++){
    xdet_near.push_back(xdetNear[ii]*CLHEP::m);
    ydet_near.push_back(ydetNear[ii]*CLHEP::m);
    zdet_near.push_back(zdetNear[ii]*CLHEP::m);
    det_near_name.push_back(detNameNear[ii]);
  }
  for(G4int ii=0;ii<nFar;ii++){
    xdet_far.push_back(xdetFar[ii]*CLHEP::m);
    ydet_far.push_back(ydetFar[ii]*CLHEP::m);
    zdet_far.push_back(zdetFar[ii]*CLHEP::m);
    det_far_name.push_back(detNameFar[ii]);
  }

  //Variables needed fr custumized dk2nu storing (Leo Aliaga, Feb17, 2015)

  // vdbl:
  //For now, just storing values for absorption studies
  nGenAbs = 3;
  nVolAbs = 4;
  G4String name_gen[] = {"parent","granparent","greatgranparent"};
  G4String name_vol[] = {"PHorn1IC","PHorn2IC","DPIP","DVOL"};
  for(G4int ii=0;ii<nGenAbs;ii++){
    GenAbsName.push_back(name_gen[ii]);
  }
  for(G4int ii=0;ii<nVolAbs;ii++){
    VolAbsName.push_back(name_vol[ii]);
  }
  char namevar[50];
  for(int ii=0; ii<nVolAbs; ii++){
    for(int jj=0; jj<nGenAbs; jj++){
      sprintf(namevar,"Material_%s_%s",VolAbsName[ii].c_str(),GenAbsName[jj].c_str());
      VolVdblName.push_back(G4String(namevar));
    }
  }
  nVdblTot = int(VolVdblName.size());

  // vint:
  //For now, storing two values
  nVintTot = 2;
  G4String name_vint[] = {"Index_Tar_In_Ancestry","playlistID"};
  for(G4int ii=0;ii<nVintTot;ii++){
    VolVintName.push_back(name_vint[ii]);
  }



     // =================================================
  // Nova Medium Energy Target
  // =================================================

  //  48 Target segments
  // + 1 Budal VFHS (Vertical Fin for Horizontal Scan)
  // + 1 Budal HFVS (Horizontal Fin for Vertical Scan)

  // Martens 3/26/10
  // Change the target width to 7.4 mm
  pSurfChk = true;
  // TargetDxdz           = 0.0; // doesn't
  // TargetDydz           = 0.0; // work properly yet
  TargetSegLength      = 24.0*CLHEP::mm;
  TargetSegWidth       = 7.4*CLHEP::mm;
  TargetSegHeight      = 63.0*CLHEP::mm;
  TargetSegPitch       = 0.5*CLHEP::mm;
  TargetGraphiteHeight = 150.0*CLHEP::mm;
  // TargetEndRounded     = true;
  // TargetSegmentNo      = 48;

  BudalVFHSLength      = TargetSegLength;
  BudalVFHSWidth       = TargetSegWidth;
  BudalVFHSHeight      = TargetSegHeight;
  BudalVFHSPitch       = 4.5*CLHEP::mm;
  BudalVFHSEndRounded  = TargetEndRounded;

  BudalHFVSLength      = TargetSegLength;
  BudalHFVSWidth       = TargetSegWidth;
  BudalHFVSHeight      = TargetSegHeight;
  BudalHFVSPitch       = 5.0*CLHEP::mm;
  BudalHFVSEndRounded  = TargetEndRounded;

  // The material "Target" is defined in NumiMaterials.cc
  // using these properties
  // TargetA              = 12.01*CLHEP::g/mole;
  // TargetZ              = 6.0;
  // TargetDensity        = 1.78*CLHEP::g/CLHEP::cm3;

  //  Z=4.,A=9.01*CLHEP::g/mole, density=1.848*CLHEP::g/CLHEP::cm3
  //Be Target
  //TargetA              = 9.01*CLHEP::g/mole;
  //TargetZ              = 4.0;
  //TargetDensity        = 1.848*CLHEP::g/CLHEP::cm3;

  // Default location of the Target wrt the MCZERO location
  // Downstream end of target is fixed at -20 cm with respect to MCZERO
  // MCZERO is at location (0,0,0) in the ROCK volume.
  // Work upstream from there to find the upstream end of the Target
  TotalTargetLength =
    TargetSegmentNo*TargetSegLength
    + (TargetSegmentNo-1)*TargetSegPitch
    + BudalVFHSLength + BudalVFHSPitch
    + BudalHFVSLength + BudalHFVSPitch;
  // TargetX0           = 0.0;
  // TargetY0           = 0.0;
  // TargetZ0           = (-20*CLHEP::cm - TotalTargetLength) + TargetConfigZ;
  // Allows TARGETZ to shift position of the target

  G4cout << "TargetZ0 = " << TargetZ0/CLHEP::cm << " cm. (Should be equal to -143.3 cm for the NOvA design.)" << G4endl;;


  // The distance between the downstream end of the target and the
  // downstream canister flange is 8 cm.
  // The centerline of the target canister is 4.15 cm below the beamline
  TargetEndtoDnFlange = 8.0*CLHEP::cm;
  TargetCanisterCenterOffset = -4.15*CLHEP::cm;

  // The X0, YO, and Z0 for the flanges and cansiter are with respect to the
  // the canister center for X0 and Y0 and the start of the target material for Z0
  // X0 and Y0 are the locations of the center of the volume, and Z0 is the location of the upstream end.
  TargetDnFlangeLength = 2.5*CLHEP::cm;
  TargetDnFlangeOutRad = 15.0*CLHEP::cm;
  TargetDnFlangeX0     = 0.0;
  TargetDnFlangeY0     = TargetCanisterCenterOffset;
  TargetDnFlangeZ0     = TotalTargetLength + TargetEndtoDnFlange;

  // Part of downstream flange is cutout to remove material
  // Part of downstream flange is cutout for the Be window
  TargetDnFlangeCutoutLength = 1.75*CLHEP::cm;
  TargetDnBeWindowRadius     = 60*CLHEP::mm;
  TargetDnBeWindowLength     = 1*CLHEP::mm;


  // The body of the target canister is modeled with three layers
  // to approximate the actual geometry given in nova-doc 3681-v3 :
  //   The outer shell is 3mm thick aluminum
  TargetOutsideCasingOutRad = 15.0*CLHEP::cm;
  TargetOutsideCasingInRad  = 14.7*CLHEP::cm;
  TargetCasingWaterOutRad   = TargetOutsideCasingInRad;
  TargetCasingWaterInRad    = 13.37*CLHEP::cm;
  TargetInsideCasingOutRad  = TargetCasingWaterInRad;
  TargetInsideCasingInRad   = 12.3*CLHEP::cm;

  TargetCasingLength = 145*CLHEP::cm;
  TargetCasingX0     = 0.0;
  TargetCasingY0     = TargetCanisterCenterOffset;
  TargetCasingZ0     = TargetDnFlangeZ0 - TargetCasingLength;


  // Upstream flange
  TargetUpFlangeLength = 2.5 * CLHEP::cm;
  TargetUpFlangeOutRad = 15.0 * CLHEP::cm;
  TargetUpFlangeX0     = 0.0;
  TargetUpFlangeY0     = TargetCanisterCenterOffset;
  TargetUpFlangeZ0     = TargetCasingZ0 - TargetUpFlangeLength;


  // Upstream Be window flange
  TargetUpBeFlangeLength = 0.5*(CLHEP::cm * 2.54) ;
  TargetUpBeFlangeOutRad = 2.73*(CLHEP::cm * 2.54) / 2.0;
  TargetUpBeFlangeX0 = 0.0;
  TargetUpBeFlangeY0 = 0.0;
  TargetUpBeFlangeZ0 =  TargetUpFlangeZ0 - TargetUpBeFlangeLength;

  TargetUpBeFlangeCutoutLength = 0.25*(CLHEP::cm * 2.54);
  TargetUpBeFlangeCutoutRadius = 1.75/2.0*(CLHEP::cm * 2.54);

  TargetUpBeWindowRadius = 0.5*(CLHEP::cm * 2.54);
  TargetUpBeWindowLength = 0.01*(CLHEP::cm * 2.54);

  // The target graphite material is clamped between the pressing plate and cooling plate
  // Pressing plate
  PressingPlateLength = TotalTargetLength;
  PressingPlateHeight = TargetGraphiteHeight - TargetSegHeight;
  PressingPlateWidth  = 20.0*CLHEP::mm;
  PressingPlateX0     = -PressingPlateWidth/2.0 - TargetSegWidth/2.0;
  PressingPlateY0     = -PressingPlateHeight/2.0 -TargetSegHeight + TargetSegWidth/2.0;
  PressingPlateZ0     = 0.0;

  PressingPlateCutoutWidth  = 5.0*CLHEP::mm;
  PressingPlateCutoutHeight = 60.0*CLHEP::mm;

  // Cooling plate
  CoolingPlateLength = TotalTargetLength;
  CoolingPlateHeight = PressingPlateHeight;
  CoolingPlateWidth  = 20.0*CLHEP::mm;
  CoolingPlateX0     = CoolingPlateWidth/2.0 + TargetSegWidth/2.0;
  CoolingPlateY0     = -CoolingPlateHeight/2.0 - TargetSegHeight + TargetSegWidth/2.0;
  CoolingPlateZ0     = 0.0;

  CoolingPlateCutoutWidth  = PressingPlateCutoutWidth;
  CoolingPlateCutoutHeight = PressingPlateCutoutHeight;

  // The water cooling runs through the cooling plate
  CoolingWaterPipeOutRad = 5.0*CLHEP::mm;
  CoolingWaterPipeX0     = 0.0*CLHEP::mm;
  CoolingWaterPipeY0     = CoolingPlateCutoutHeight/2.0 + (CoolingPlateHeight - CoolingPlateCutoutHeight)/4.0;



  NumiDataInput::Print();

  usePosLocalCoordInMagField = true;
  useRotLocalCoordInMagField = true;

  if(fPrintInfo > 0 || debugOn)
  {
     G4cout << "********************************************************************" << G4endl;
     G4cout << "...NumiDataInput Initialization Complete." << G4endl;
     G4cout << "********************************************************************" << G4endl;
     G4cout << "********************************************************************" << G4endl;
     G4cout << G4endl;
     G4cout << G4endl;

  }
}
//---------------------------------------------------------------------------------
NumiDataInput::~NumiDataInput()
{
}


//---------------------------------------------------------------------------------
void NumiDataInput::Print()
{
   G4cout << "NumiDataInput::Print() - Printing Parameters..." << G4endl;
   G4cout << " SIMULATION                                  = " << fSimulation << G4endl;
   if(fUseNuBeam || fSimulation.empty())
   {
      if(fUseWaterInTgt)
      {
         G4cout << " SUBSIMULATION                               = " << fSubSimulation << G4endl;
         G4cout << "    Simulating " << fLengthOfWaterInTgt/CLHEP::cm << " cm of water filling the "
                << "target starting from the downstream tip." << G4endl;

      }
      G4cout << " Run Period                                  = " << fRunPeriod << G4endl
             << " Beam Configuration                          = " << fBeamConfig << G4endl
	     << " Playlist                                    = " << fPlaylist << G4endl
             << " Horn Configuration                          = " << fHornConfig << G4endl
             << " Target Configuration                        = " << fTargetConfig << G4endl
             << " Horn Current Configuration                  = " << fIHornConfig << G4endl
             << " Target               (X0, Y0, Z0) m         = (" << TargetX0/CLHEP::m << ", " << TargetY0/CLHEP::m << ", " << TargetZ0/CLHEP::m << ") m" << G4endl
             << " Baffle               (X0, Y0, Z0) m         = (" << HPBaffleX0/CLHEP::m << ", " << HPBaffleY0/CLHEP::m << ", " << HPBaffleZ0/CLHEP::m << ") m" << G4endl
             << " Horizontal Fin w.r.t. Target (X0, Y0, Z0) m = (" << BudalX0/CLHEP::m << ", " << BudalY0/CLHEP::m << ", " << BudalZ0/CLHEP::m << ") m" << G4endl
             << " Horn Current                                = " << HornCurrent/CLHEP::ampere << " A" << G4endl
             << " Horn 1               (X0, Y0, Z0) m         = (" << Horn1X0/CLHEP::m << ", " << Horn1Y0/CLHEP::m << ", " << Horn1Z0/CLHEP::m << ") m" << G4endl
	     << " Horn 1               (Phi, Theta, Psi) rad  = (" << Horn1Phi/CLHEP::rad << ", " << Horn1Theta/CLHEP::rad << ", " << Horn1Psi/CLHEP::rad << ") rad" << G4endl
             << " Horn 2               (X0, Y0, Z0) m         = (" << Horn2X0/CLHEP::m << ", " << Horn2Y0/CLHEP::m << ", " << Horn2Z0/CLHEP::m << ") m" << G4endl
	     << " Horn 2               (Phi, Theta, Psi) rad  = (" << Horn2Phi/CLHEP::rad << ", " << Horn2Theta/CLHEP::rad << ", " << Horn2Psi/CLHEP::rad << ") rad" << G4endl
             << " Target Density                              = " << TargetDensity/CLHEP::g*CLHEP::cm3 << " g/CLHEP::cm^3" << G4endl
             << " Proton Beam Position (X0, Y0, Z0) m         = (" << beamPosition[0]/CLHEP::m << ", " << beamPosition[1]/CLHEP::m << ", " << beamPosition[2]/CLHEP::m << ") m" << G4endl
             << " Proton Beam Momentum                        = " << protonMomentum/CLHEP::GeV << " GeV/c" << G4endl
             << " Proton Beam X-Sigma                         = " << beamSigmaX/CLHEP::mm << " mm" << G4endl
             << " Proton Beam Y-Sigma                         = " << beamSigmaY/CLHEP::mm << " mm" << G4endl
	     << " Decay Pipe Helium                           = " << std::boolalpha << HeInDecayPipe<<G4endl
	     << " Decay Pipe Magnetic Field                   = " << std::boolalpha << applyDecayPipeMagneticField << G4endl;

      if(fUseDetailedProtonBeam)
      {
         G4cout << " Using Detailed Proton Beam with parameters..." << G4endl
                << "    Proton Beam Mean Momentum = " << protonMomentum/CLHEP::GeV << " GeV/c" << G4endl
                << "    Proton out_R              = " << fProton_outR << " probably mm " << G4endl
                << "    Proton in_R               = " << fProton_inR << " probably mm " << G4endl
                << "    Proton Divergence         = " << fProtonDiv << " probably rad " << G4endl
                << "    Proton Spread             = " << fProtonSpread  << " probably GeV/c " << G4endl
                << "    Proton cosx               = " << fProton_cosx  << G4endl
                << "    Proton cosy               = " << fProton_cosy  << G4endl;
      }



   }
   if(useMuonBeam || fSimulation.empty())
   {
   }
   if(simAbsBkg || fSimulation.empty())
   {
   }
}

//---------------------------------------------------------------------------------
bool NumiDataInput::SetBeamConfig(G4String config)
{

   if(fRunPeriod == -999)
   {
      G4cout << "NumiDataInput::SetBeamConfig() - PROBLEM: Run Period not set. MUST SET RUN PERIOD BEFORE SETTING BEAM CONFIGURATION. " << G4endl;
      return false;
   }

   std::string::size_type beg_loc_tgt = std::string::npos;
   std::string::size_type size_horn;

   vstring_t tgtvec;
   tgtvec.push_back("le");
   tgtvec.push_back("me");

   for(vstring_t::const_iterator git = tgtvec.begin(); git != tgtvec.end(); ++git)
   {
      beg_loc_tgt = config.find(*git, 0);
      if(beg_loc_tgt != std::string::npos){ size_horn = (*git).size(); break;}
   }

   if(beg_loc_tgt == std::string::npos)
   {
      G4cout << "NumiDataInput::SetBeamConfig() - PROBLEM. Can't get beam configuration. "
             << "Beam configuration must be in the form \"LE#z#i\" or \"le#z#i\", where # is any number. "
             << "Examples are le010z185i, LE025.3z-200i, LE250z185.6i....etc."
             << "Note that I don't know how to simulate \"ME\" beam configurations right now." << G4endl;
      return false;
   }

   //
   //set the horn config
   //

   if(!NumiDataInput::SetHornConfig(config.substr(0,2))) return false;

   //
   //get the target config
   //

   std::string::size_type end_loc_tgt = config.find("z", 0);

   if(end_loc_tgt == std::string::npos)
   {
      G4cout << "NumiDataInput::SetBeamConfig() - PROBLEM. Can't find \"z\". Can't get beam configuration. "
             << "Beam configuration must be in the form \"LE#z#i\" or \"le#z#i\", where # is any number. "
             << "Examples are le010z185i, LE025.3z-200i, LE250z185.6i....etc." << G4endl;
      return false;
   }


   std::string::size_type size_tgt = end_loc_tgt-(beg_loc_tgt+size_horn) +1;
   std::string tgtzstr = config.substr(beg_loc_tgt+size_horn, size_tgt);

   //
   //set the target config and location
   //

   if(!NumiDataInput::SetTargetConfig(tgtzstr)) return false;


   //
   //get horn current
   //

   std::string::size_type end_loc_ihorn = config.find("i", 0);

   if(end_loc_ihorn == std::string::npos)
   {
      G4cout << "NumiDataInput::SetBeamConfig() - PROBLEM. Can't find \"i\". Can't get beam configuration. "
             << "Beam configuration must be in the form \"LE#z#i\" or \"le#z#i\", where # is any number. "
             << "Examples are le010z185i, LE025.3z-200i, LE250z185.6i....etc." << G4endl;
      return false;
   }

   std::string ihornstr = config.substr(end_loc_tgt+1, end_loc_ihorn-end_loc_tgt);


   if(!NumiDataInput::SetHornCurrentConfig(ihornstr)) return false;
//These lines are to change the horizontal transverse position of the position of the horns and the target density --Bruce Howard, Jr.
   if(fUseHornMisalign){
      //
      //get horn 1 X0
      //
      std::string::size_type end_loc_h1 = config.find("hot",0);

      if(end_loc_h1 == std::string::npos)
      {
         G4cout << "\n\n\n\n\n HORN 1 POSITION NOT PROPERLY ENTERED!! MUST ENTER AS ###hot. PROGRAM ABORTED." << G4endl;
         std::exit (EXIT_FAILURE);
      }

      std::string ih1str = config.substr(end_loc_h1-3,3); //Grabbing the characters before end_loc_h1 (ie USE  ###hot)

      if(!NumiDataInput::SetHornOnePos(ih1str)) return false;

      //
      //get horn 2 X0
      //
      std::string::size_type end_loc_h2 = config.find("htt",0);

      if(end_loc_h2 == std::string::npos)
      {
         G4cout << "\n\n\n\n\n HORN 2 POSITION NOT PROPERLY ENTERED!! MUST ENTER AS ###htt. PROGRAM ABORTED." << G4endl;
         std::exit (EXIT_FAILURE);
      }

      std::string ih2str = config.substr(end_loc_h2-3,3);

      if(!NumiDataInput::SetHornTwoPos(ih2str)) return false;
   }
   else
   {
     Horn1X0=0.0;
     Horn2X0=0.0;
   }
   if(fUseTgtDensity){
      std::string::size_type end_loc_tgd = config.find("tgd",0);

      if(end_loc_tgd==std::string::npos)
      {
         G4cout << "\n\n\n\n\n TARGET DENSITY NOT PROPERLY ENTERED!! MUST ENTER AS #######tgd, in ug/CLHEP::cm3. PROGRAM ABORTED." <<G4endl;
         std::exit (EXIT_FAILURE);
      }

      std::string targetdensity = config.substr(end_loc_tgd-7,7);
      if(!NumiDataInput::SetTargetDensity(targetdensity)) return false;
   }
   if(!fUseTgtDensity)
   {
      TargetDensity=1.78*CLHEP::g/CLHEP::cm3;
   }
//--------------------------------------------------------------------------------------------------------------


   if(!NumiDataInput::ConfigureRunPeriod(config)) return false;


   fBeamConfig = config;



   return true;

}

//---------------------------------------------------------------------------------
G4bool NumiDataInput::SetPlaylist(G4String pl)
{
   fPlaylist = pl;

   return true;
}

//---------------------------------------------------------------------------------
G4bool NumiDataInput::SetTargetConfig(G4String config)
{
   const std::string::size_type loc_z = config.find("z", 0);

   const std::string::size_type size_config = config.size();

   if(loc_z != size_config-1)
   {
      G4cout << "NumiDataInput::SetTargetConfig() - PROBLEM. Target Configuration, " << config << " is not "
             << " in the form \"#z\", where # is a number. Can't get target configuration. " << G4endl;
      return false;
   }

   const std::string tgtzstr = config.substr(0,size_config-1);

   const std::string number = "0123456789.";

   for(unsigned int i = 0; i < tgtzstr.size(); ++i)
   {
      bool isnumber = false;
      for(unsigned int j = 0; j < number.size(); ++j)
      {
         if(tgtzstr[i] == number[j])
         {
            isnumber = true;
            break;
         }
      }
      if(isnumber == false)
      {
         G4cout << "NumiDataInput::SetTargetConfig() - PROBLEM. Target Z positon, " << tgtzstr << " is not a number. "
                << " Can't get target configuration. "
                << "Target configuration must be in the form \"#z\", where # is any number. "
                << "Examples are 010z, 100z, 025.3z, 250.5z....etc." << G4endl;
         return false;
      }
   }


   G4double tgtz;
   std::istringstream tgtzstrm(tgtzstr);
   tgtzstrm >> tgtz;

   NumiDataInput::SetTargetZ0(fTargetZ0_ref   - tgtz*CLHEP::cm);
   NumiDataInput::SetBaffleZ0(fHPBaffleZ0_ref - tgtz*CLHEP::cm);
   NumiDataInput::SetBeamZ0  (fBeamZ0_ref     - tgtz*CLHEP::cm);

   fTargetConfig = config;


   return true;
}

//---------------------------------------------------------------------------------
G4bool NumiDataInput::SetHornCurrentConfig(G4String config)
{
   const std::string::size_type loc_i = config.find("i", 0);

   const std::string::size_type size_config = config.size();

   if(loc_i != size_config-1)
   {
      G4cout << "NumiDataInput::SetHornCurrentConfig() - PROBLEM. Horn Current Configuration, " << config << " is not "
             << " in the form \"#i\", where # is a number. Can't get horn current configuration. " << G4endl;
      return false;
   }

   std::string ihornstr = config.substr(0,size_config-1);

   const std::string number = "-0123456789.";

   for(unsigned int i = 0; i < ihornstr.size(); ++i)
   {
      bool isnumber = false;
      for(unsigned int j = 0; j < number.size(); ++j)
      {
         if(ihornstr[i] == number[j])
         {
            isnumber = true;
            break;
         }
      }
      if(isnumber == false)
      {
         G4cout << "NumiDataInput::SetBeamConfig() - PROBLEM. Horn Current, " << ihornstr << " is not a number. "
                << "Can't get horn current configuration. "
                << "Horn current configuration must be in the form \"#i\", where # is any number. "
                << "Examples are 185i, -059i, -200i, 20.6i....etc." << G4endl;
         return false;
      }
   }

   if(fUseCorrHornCurrent && fHornConfig == "le")
   {
      //
      //calibration corrections from MINOS-doc-1303
      //
           if (ihornstr == "200")   ihornstr = "196.8";
      else if (ihornstr == "-200")  ihornstr = "-196.8";
      else if (ihornstr == "185")   ihornstr = "182.1";
      else if (ihornstr == "-185")  ihornstr = "-182.1";
      else if (ihornstr == "170")   ihornstr = "167.3";
      else if (ihornstr == "-170")  ihornstr = "-167.3";
      else if (ihornstr == "0"   ||
               ihornstr == "000" ||
               ihornstr == "0.0" ||
               ihornstr == "000.0")  ihornstr = "0";
//--------------------------------The following ihornstr added by Bruce Howard, Jr.---------------------------------
//-----------These provide changes up and down from positive and negative 185kA, stepped by 1% of 185kA-------------
      else if (ihornstr == "172")   ihornstr = "172.995";
      else if (ihornstr == "174")   ihornstr = "174.816";
      else if (ihornstr == "176")   ihornstr = "176.637";
      else if (ihornstr == "178")   ihornstr = "178.458";
      else if (ihornstr == "180")   ihornstr = "180.279";
      else if (ihornstr == "183")   ihornstr = "183.921";
      else if (ihornstr == "186")   ihornstr = "185.742";
      else if (ihornstr == "187")   ihornstr = "187.563";
      else if (ihornstr == "189")   ihornstr = "189.384";
      else if (ihornstr == "191")   ihornstr = "191.205";
      else if (ihornstr == "-172")   ihornstr = "-172.995";
      else if (ihornstr == "-174")   ihornstr = "-174.816";
      else if (ihornstr == "-176")   ihornstr = "-176.637";
      else if (ihornstr == "-178")   ihornstr = "-178.458";
      else if (ihornstr == "-180")   ihornstr = "-180.279";
      else if (ihornstr == "-183")   ihornstr = "-183.921";
      else if (ihornstr == "-186")   ihornstr = "-185.742";
      else if (ihornstr == "-187")   ihornstr = "-187.563";
      else if (ihornstr == "-189")   ihornstr = "-189.384";
      else if (ihornstr == "-191")   ihornstr = "-191.205";
//-------------------------------------------------------------------------------------------------------------------
      else
      {
         G4cout << "NumiDataInput::SetHornConfig() - PROBLEM. Requesting to use the Corrected Horn Current but I do not know "
                << "the correction for a horn current of " << ihornstr << ". I only know the corrections for horn currents of "
                << "+/- 200kA, +/- 185kA, and +/- 170kA (and 0kA). Can't get horn current configuration. " << G4endl;
         return false;
      }
   }

   G4double ihorn;
   std::istringstream ihornstrm(ihornstr);
   ihornstrm >> ihorn;
   //
   //note the horn current must be stored in amps
   //

   NumiDataInput::SetHornCurrent(ihorn*1000.*CLHEP::ampere);
   //
   // Do something naughty: Multiply by 1000 for BFieldMu geantino analysis.
   //
   if (fDumpBFieldPlease) {
       NumiDataInput::SetHornCurrent(ihorn*1.e6*CLHEP::ampere);
       std::cerr << " Horn Current has been multiplied by 1 e3 for Muon Geantino use, Hor Current is now  "
                 << this->HornCurrent << std::endl;
//       exit(2);
   }
   fIHornConfig = ihornstr + "i";

   return true;
}

//---------------------------------------------------------------------------------
G4bool NumiDataInput::SetHornOnePos(G4String config)
{
   G4double hornpos;
   std::istringstream posstream(config);
   posstream >> hornpos;

   NumiDataInput::SetHorn1X0((hornpos/10.)*CLHEP::cm); //Input as mm but we want cm

   return true;
}
//------------------------------------------------------------------------------------
G4bool NumiDataInput::SetHornTwoPos(G4String config)
{
   G4double hornpos;
   std::istringstream posstream(config);
   posstream >> hornpos;

   std::cerr << " Setting Horn2 position at " << hornpos << std::endl;

   NumiDataInput::SetHorn2X0((hornpos/10.)*CLHEP::cm); //Input as mm, but we want cm

   return true;
}
//------------------------------------------------------------------------------------
G4bool NumiDataInput::SetTargetDensity(G4String config)
{
   G4double targetdensity;
   std::istringstream posstream(config);
   posstream >> targetdensity;

   NumiDataInput::SetTargetDensityVal((targetdensity/1000000.)*CLHEP::g/CLHEP::cm3);

   return true;
}
//------------------------------------------------------------------------------------
G4bool NumiDataInput::SetHornConfig(G4String config)
{
   if(config != "le" && config != "me")
   {
      G4cout << "NumiDataInput::SetHornConfig() - PROBLEM. Can't set horn configuration. "
             << "Valid horn configs are \"LE\" or \"le\"."
             << "Note that I don't know how to simulate \"ME\" beam configurations right now." << G4endl;

      return false;
   }

   fHornConfig = config;

   return true;
}

//---------------------------------------------------------------------------------
G4bool NumiDataInput::ConfigureRunPeriod(G4String &beamconfig)
{
   if(fRunPeriod == 0)
   {
      return true;
   }
   else if(fRunPeriod == 1)
   {
      if(beamconfig != "LE010z000i" && beamconfig != "le010z000i" &&
         beamconfig != "LE010z170i" && beamconfig != "le010z170i" &&
         beamconfig != "LE010z185i" && beamconfig != "le010z185i" &&
         beamconfig != "LE010z200i" && beamconfig != "le010z200i" &&
         beamconfig != "LE100z200i" && beamconfig != "le100z200i" &&
         beamconfig != "LE250z200i" && beamconfig != "le250z200i")
      {
         G4cout << "NumiDataInput::ConfigureRunPeriod() - PROBLEM: No beam configuration " << beamconfig
                << " during run period " << fRunPeriod << G4endl;
         return false;
      }

      TargetY0 = -1.1*CLHEP::mm;
      BudalY0  = 2.26*CLHEP::mm;

      return true;
   }
   else if (fRunPeriod == 2)
   {
      if(beamconfig != "LE010z000i" && beamconfig != "le010z000i" &&
         beamconfig != "LE010z185i" && beamconfig != "le010z185i" &&
         beamconfig != "LE150z200i" && beamconfig != "le150z200i" &&
         beamconfig != "LE250z200i" && beamconfig != "le250z200i")
      {
         G4cout << "NumiDataInput::ConfigureRunPeriod() - PROBLEM: No beam configuration " << beamconfig
                << " during run period " << fRunPeriod << G4endl;
         return false;
      }

      // Change LE010 to LE008.9
      if (beamconfig.find("LE010") != std::string::npos ||
          beamconfig.find("le010") != std::string::npos )
      {
         if(!NumiDataInput::SetTargetConfig("008.9z")) return false;
      }

      return true;
   }
   else if(fRunPeriod == 3)
   {
      if(beamconfig != "LE010z000i" && beamconfig != "le010z000i" &&
         beamconfig != "LE010z185i" && beamconfig != "le010z185i" )
      {
         G4cout << "NumiDataInput::ConfigureRunPeriod() - PROBLEM: No beam configuration " << beamconfig
                << " during run period " << fRunPeriod << G4endl;
         return false;
      }

      // Change LE010 to LE008.9
      if (beamconfig.find("LE010") != std::string::npos ||
          beamconfig.find("le010") != std::string::npos )
      {
         if(!NumiDataInput::SetTargetConfig("008.9z")) return false;
      }

      return true;
   }
   else if(fRunPeriod == 4)
   {
      if(beamconfig != "LE010z-185i" && beamconfig != "le010z-185i" )
      {
         G4cout << "NumiDataInput::ConfigureRunPeriod() - PROBLEM: No beam configuration " << beamconfig
                << " during run period " << fRunPeriod << G4endl;
         return false;
      }

      return true;
   }
   else
   {

      G4cout << "NumiDataInput::ConfigureRunPeriod() - PROBLEM: Don't know about Run Period " << fRunPeriod << G4endl;
      return false;
   }



}

//---------------------------------------------------------------------------------
void NumiDataInput::SetDetailedProtonBeam(G4bool val)
{
   if(val)
   {
      //
      //These values are in g4numi_flugg/scripts/g4numi_fluka.sh
      //
      //protonMomentum = 120.0*CLHEP::GeV;
      fProton_outR  = -0.26350;  //this is probably mm
      fProton_inR   = -0.26139; //this is probably mm
      fProtonDiv    = -0.02355; //this is probably rads
      fProtonSpread = -0.000001; //this is probably GeV?

      //
      //these values are in NumiPrimaryMessenger.cc
      //
      fProton_cosx = 0.0;
      fProton_cosy = 0.0;
   }

   fUseDetailedProtonBeam = val;
}

//---------------------------------------------------------------------------------
void NumiDataInput::SetLengthOfWaterInTgt(G4double val)
{
   if(val < 3.0*CLHEP::cm)
   {
      G4cout << " NumiDataInput::SetLengthOfWaterInTgt() - PROBLEM: Can't have less than "
             << "3 cm of water filling the end of the target. Setting water length to 3 cm" << G4endl;
      fLengthOfWaterInTgt = 3.0*CLHEP::cm;
   }
   else
      fLengthOfWaterInTgt = val;
}



//---------------------------------------------------------------------------------
void NumiDataInput::ApplyStepLimits(G4LogicalVolume *vol)
{
  if (StepLimit == 0.0) return;
  vol->SetUserLimits(new G4UserLimits(StepLimit));
}

//---------------------------------------------------------------------------------
void NumiDataInput::ApplyTimeLimits(G4LogicalVolume *vol)
{
  if (TimeLimit == 0.0) return;
  vol->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,TimeLimit));
}

//---------------------------------------------------------------------------------
void NumiDataInput::SetAbsorberMaterial(G4Material* mat, G4int mon)
{
   if(mon == 0)
      Mon1AbsorberMaterial = mat;
   else if(mon == 1)
      Mon2AbsorberMaterial = mat;
   else if(mon == 2)
      Mon3AbsorberMaterial = mat;
   else
   {
      G4cout << "******Problem : NumiDataInput::SetAbsorberMaterial - Invalid Monitor number " << mon
             << " Setting absorber material to " << DefaultMaterial->GetName() << G4endl;
      Mon1AbsorberMaterial = DefaultMaterial;
      Mon2AbsorberMaterial = DefaultMaterial;
      Mon3AbsorberMaterial = DefaultMaterial;
   }

}

//---------------------------------------------------------------------------------
void NumiDataInput::SetAbsorberThickness(G4double val, G4int mon)
{
   if(mon == 0)
      Mon1AbsorberThickness = val;
   else if(mon == 1)
      Mon2AbsorberThickness = val;
   else if(mon == 2)
      Mon3AbsorberThickness = val;
   else
   {
      G4cout << "******Problem : NumiDataInput::SetAbsorberThickness - Invalid Monitor number " << mon
             << " Setting absorber thickness to 0.0 mm" << G4endl;
      Mon1AbsorberThickness = 0.0;
      Mon2AbsorberThickness = 0.0;
      Mon3AbsorberThickness = 0.0;
   }

}

//----------------------------------------------------------------------------------------
void NumiDataInput::SetAbsorberMonDist(G4double val, G4int mon)
{
   if(mon == 0)
      Mon1AbsorberDist = val;
   else if(mon == 1)
      Mon2AbsorberDist = val;
   else if(mon == 2)
      Mon3AbsorberDist = val;
   else
   {
      G4cout << "******Problem : NumiDataInput::SetAbsorberMonDist - Invalid Monitor number " << mon
             << " Setting absorber distance to 1 m" << G4endl;
      Mon1AbsorberDist = 1.0*CLHEP::m;
      Mon2AbsorberDist = 1.0*CLHEP::m;
      Mon3AbsorberDist = 1.0*CLHEP::m;
   }

}


//----------------------------------------------------------------------
G4Material* NumiDataInput::GetAbsorberMaterial(G4int mon)
{
   if(mon == 0) return Mon1AbsorberMaterial;
   else if(mon == 1) return Mon2AbsorberMaterial;
   else if(mon == 2) return Mon3AbsorberMaterial;
   else
   {
      G4cout << "******Problem : NumiDataInput::GetAbsorberMaterial - Invalid Monitor number " << mon
             << " returing " << DefaultMaterial->GetName() << G4endl;
      return DefaultMaterial;
   }

}

//----------------------------------------------------------------------
G4double NumiDataInput::GetAbsorberThickness(G4int mon)
{
   if(mon == 0) return Mon1AbsorberThickness;
   else if(mon == 1) return Mon2AbsorberThickness;
   else if(mon == 2) return Mon3AbsorberThickness;
   else
   {
      G4cout << "******Problem : NumiDataInput::GetAbsorberThickness - Invalid Monitor number " << mon
             << " returning 1e-6 mm" << G4endl;
      return 0.0;
   }

}

//----------------------------------------------------------------------------------------
G4double NumiDataInput::GetAbsorberMonDist(G4int mon)
{
   if(mon == 0) return Mon1AbsorberDist;
   else if(mon == 1) return Mon2AbsorberDist;
   else if(mon == 2) return Mon3AbsorberDist;
   else
   {
      G4cout << "******Problem : NumiDataInput::GetAbsorberMonDist - Invalid Monitor number " << mon
             << " returning  1 m" << G4endl;
      return 1.0*CLHEP::m;
   }

}
//--------------------------------------------------------------------------------
void NumiDataInput::SetjCompare(G4bool _jc) {
    jCompare = _jc;
    G4double PHorn1EndZ0_[]     ={126.092        ,127.842     ,129.718};
    if (jCompare) {
        PHorn1EndZ0_[0] = 118.11;
        PHorn1EndZ0_[1] = 119.86;
        PHorn1EndZ0_[2] = 122.048;
    }
    PHorn1EndZ0.clear();
    for (G4int ii=0;ii<NPHorn1EndN;ii++){
        PHorn1EndZ0.push_back(PHorn1EndZ0_[ii]*in);
    }


    G4double PHorn2EndZ0_[]     ={135.861        ,137.611     ,139.486};
    if (jCompare) {
        PHorn2EndZ0_[0] = 118.11;
        PHorn2EndZ0_[1] = 119.86;
        PHorn2EndZ0_[2] = 122.048;
    }
    PHorn2EndZ0.clear();
    for (G4int ii=0;ii<NPHorn2EndN;ii++){
        PHorn2EndZ0.push_back(PHorn2EndZ0_[ii]*in);
    }
}

//---------------------------------------------------------------------------------
void NumiDataInput::Setg3Chase(G4bool _gc)
{
    g3Chase = _gc;

    if(g3Chase){
        THBlockNblock = 24;
    }
    else{
        THBlockNblock=18;
    }

}

//---------------------------------------------------------------------------------
void NumiDataInput::SetBeamSigmaX(G4double val)
{
  beamSigmaX = val;
}

//---------------------------------------------------------------------------------
void NumiDataInput::SetBeamSigmaY(G4double val)
{
  beamSigmaY = val;
}


