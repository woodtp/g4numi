// --------------------------------------------------------------
//NumiMagneticField by Yuki 7/12/04
//modified by Yuki 8/2/04

#include "NumiMagneticField.hh"
#include "G4ios.hh"
#include "G4TransportationManager.hh"
#include "NumiDataInput.hh" 
#include "G4RunManager.hh"
#include "G4VPhysicalVolume.hh"

//magnetic field between conductors ====================================================
NumiMagneticField::NumiMagneticField()
{
  NumiData=NumiDataInput::GetNumiDataInput();
  if(NumiData) current = NumiData->HornCurrent;
  else current=205.; //kA
}

NumiMagneticField::~NumiMagneticField(){;}

void NumiMagneticField::GetFieldValue(const double Point[3],double *Bfield) const
{
  G4double radius = sqrt(Point[0]*Point[0]+Point[1]*Point[1]);    
  G4double B = current / (5.*radius/cm)/10*tesla;  //B(kG)=i(kA)/[5*r(cm)], 1T=10kG

  Bfield[0] = -B*Point[1]/radius;
  Bfield[1] = B*Point[0]/radius;
  Bfield[2] = 0.;
}

//magnetic field in inner conductor ====================================================
NumiMagneticFieldIC::NumiMagneticFieldIC()
{
  NumiData=NumiDataInput::GetNumiDataInput();
  if(NumiData) current = NumiData->HornCurrent;
  else current=205.; //kA 
}

NumiMagneticFieldIC::~NumiMagneticFieldIC(){;}

void NumiMagneticFieldIC::GetFieldValue(const double Point[3],double *Bfield) const
{
  G4Navigator* numinavigator=new G4Navigator(); //geometry navigator
  G4Navigator* theNavigator=G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  numinavigator->SetWorldVolume(theNavigator->GetWorldVolume());
  G4ThreeVector Position=G4ThreeVector(Point[0],Point[1],Point[2]); 
  G4VPhysicalVolume* myVolume = numinavigator->LocateGlobalPointAndSetup(Position);
  G4TouchableHistoryHandle aTouchable = numinavigator->CreateTouchableHistoryHandle();
  G4ThreeVector localPosition = aTouchable->GetHistory()->GetTopTransform().TransformPoint(Position);

  delete numinavigator;

  G4VSolid * solid=myVolume->GetLogicalVolume()->GetSolid();

  G4double radius = sqrt(Point[0]*Point[0]+Point[1]*Point[1]);
  G4double d_out=0.;
  G4double d_in=0.;
  G4double Bfield_mag = 0.;  

  if (myVolume->GetName().contains("PI")){
  d_out=solid->DistanceToOut(localPosition,G4ThreeVector(Point[0]/radius,Point[1]/radius,0)); //distance to outer boundary
  d_in=solid->DistanceToOut(localPosition,G4ThreeVector(-Point[0]/radius,-Point[1]/radius,0));//distance to inner boundary
  if (d_out<1.*cm&&d_in<1.*cm) 
    {
    Bfield_mag = current / (5.*radius/cm)/10*tesla; //B(kG)=i(kA)/[5*r(cm)], 1T=10kG
    Bfield_mag=Bfield_mag*(1-(radius*radius-(radius-d_in)*(radius-d_in))/((radius+d_out)*(radius+d_out)-(radius-d_in)*(radius-d_in)));// linear distribution of current
    }
  }
 
  Bfield[0] = -Bfield_mag*Point[1]/radius;
  Bfield[1] = Bfield_mag*Point[0]/radius;
  Bfield[2] = 0.;
}


//magnetic field in outter conductor====================================================
NumiMagneticFieldOC::NumiMagneticFieldOC()
{
  NumiData=NumiDataInput::GetNumiDataInput();
  if(NumiData) current = NumiData->HornCurrent;
  else current=205.; //kA

 
}

NumiMagneticFieldOC::~NumiMagneticFieldOC(){;}

void NumiMagneticFieldOC::GetFieldValue(const double Point[3],double *Bfield) const
{
  G4Navigator* numinavigator=new G4Navigator(); //geometry navigator
  G4Navigator* theNavigator=G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  numinavigator->SetWorldVolume(theNavigator->GetWorldVolume());
  G4ThreeVector Position=G4ThreeVector(Point[0],Point[1],Point[2]); 
  G4VPhysicalVolume* myVolume = numinavigator->LocateGlobalPointAndSetup(Position);
  G4TouchableHistoryHandle aTouchable = numinavigator->CreateTouchableHistoryHandle();
  G4ThreeVector localPosition = aTouchable->GetHistory()->GetTopTransform().TransformPoint(Position);

  delete numinavigator;

  G4double rin=15.33*cm;
  G4double rout=16.2*cm;
  G4double radius = sqrt(Point[0]*Point[0]+Point[1]*Point[1]); //cm
  if ((myVolume->GetName().contains("PI06"))||
      (myVolume->GetName().contains("PI07"))||
      (myVolume->GetName().contains("PI08")))
    {
      rin=37.*cm;
      rout=37.86*cm; 
    }
  G4double B = current / (5.*radius/cm)/10*tesla;  //B(kG)=i(kA)/[5*r(cm)], 1T=10kG
  B=B*(1-(radius*radius-rin*rin)/(rout*rout-rin*rin)); // linear distribution of current

  Bfield[0] = -B*Point[1]/radius;
  Bfield[1] = B*Point[0]/radius;
  Bfield[2] = 0.;
}
