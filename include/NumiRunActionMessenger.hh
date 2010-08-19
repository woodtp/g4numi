#ifndef NumiRunActionMessenger_h
#define NumiRunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ios.hh"

class NumiRunAction;
class NumiDataInput;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithoutParameter;

class NumiRunActionMessenger: public G4UImessenger
{
public:
  NumiRunActionMessenger(NumiRunAction* );
  ~NumiRunActionMessenger();
  
  void SetNewValue(G4UIcommand* ,G4String );
 
private:
   NumiRunAction*              runAction;
   NumiDataInput*              NumiData;
   G4UIdirectory*              RndmDir;
   G4UIcmdWithAString*         readRndmCmd;  
   G4UIcmdWithoutParameter*    showRndmCmd;
   G4UIcmdWithAnInteger*       setRndmSeedCmd;
   G4UIdirectory*              NumiRunDir;
   G4UIcmdWithABool*           debugOn;
   G4UIcmdWithAnInteger*       setRunID;
   G4UIcmdWithADouble*         setMaterialSigma;
   G4UIcmdWithABool*           useNImpWeight;
   G4UIcmdWithABool*           useTestBeam;
   G4UIcmdWithABool*           useFlukaInput;
   G4UIcmdWithABool*           useMarsInput;
   G4UIcmdWithABool*           useMuonBeam;
   G4UIcmdWithABool*           simAbsBkg;
   G4UIcmdWithABool*           simDRays;
   G4UIcmdWithABool*           useZPosCut;
   G4UIcmdWithAString*         muonBeamShape;
   G4UIcmdWithADoubleAndUnit*  muonBeamMomentum;
   G4UIcmdWithADoubleAndUnit*  muonBeamZPos;
   G4UIcmdWithADoubleAndUnit*  muonBeamGaussXsig;
   G4UIcmdWithADoubleAndUnit*  muonBeamGaussYsig;
   G4UIcmdWithABool*           useMuonInput;
   G4UIcmdWithAnInteger*       NInputParts;
   G4UIcmdWithAnInteger*       NInputPart;
   G4UIcmdWithABool*           reWeightDeltas;
   G4UIcmdWithAnInteger*       nSplitDeltas;
   G4UIcmdWithAString*         extNtupleFileName;
   G4UIcmdWithAString*         setHadmmNtupleDir;
   G4UIcmdWithAString*         setAbsBkgNtupleDir;
   G4UIdirectory*              NumiOutputDir;
   G4UIcmdWithAString*         setNuNtupleFile; 
   G4UIcmdWithAString*         setHadmmNtupleFile;
   G4UIcmdWithAString*         setAbsBkgNtupleFile;
   G4UIcmdWithAnInteger*       setNEvents;
   G4UIcmdWithAString*         setASCIIFile; 
   G4UIcmdWithAString*         setBXDRAWFile; 
   G4UIcmdWithABool*           outputNuNtuple;
   G4UIcmdWithABool*           outputHadmmNtuple;
   G4UIcmdWithABool*           outputAbsBkgNtuple;
   G4UIcmdWithABool*           outputASCIIFile;
   G4UIcmdWithABool*           outputBXDRAWFile;
   G4UIcmdWithADouble*         setTestTheta;
   G4UIcmdWithABool*           useDecayPipeSelect;
   G4UIcmdWithADoubleAndUnit*  setStepLimit;
   
   
   G4UIcmdWithABool*           UseMacro;
   G4UIcmdWithABool*           outputZpNtuple;
   G4UIcmdWithAString*         setZpNtupleFile;
   G4UIcmdWithABool*           KillTracking;
   G4UIcmdWithADoubleAndUnit*  KillTrackingThreshold;
   
};

#endif
