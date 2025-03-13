//
// NumiTrajectory.cc
//

#include "NumiTrajectory.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4Polyline.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"
#include "G4VProcess.hh"
#include "NumiTrackInformation.hh"
#include "NumiDataInput.hh"


#include "G4ProcessType.hh"
#include "G4HadronicProcess.hh"
// for "EM" nuclear interactions, as "muonNuclear" Process occurs.
#include "G4MuonNuclearProcess.hh"
#include "G4ElectronNuclearProcess.hh"
#include "G4PositronNuclearProcess.hh"

G4Allocator<NumiTrajectory> myTrajectoryAllocator;

NumiTrajectory::NumiTrajectory()
{
   fParticleDefinition = 0;
   fParticleName       = "";
   fPDGCharge          = 0;
   fPDGEncoding        = 0;
   fTrackID            = 0;
   fParentID           = 0;
   fPositionRecord     = 0;
   fMomentum           = G4ThreeVector(0.,0.,0.);
   fMomentumRecord     = 0;
   fVertexPosition     = G4ThreeVector(0.,0.,0.);
   fParticleMass       = 0.;
   fDecayCode          = 0;
   fTgen               = 0;
   fNImpWt             = 1.;
   fPreStepVolume      = 0;
   fStepLength         = 0;
   fMaterialNumber1rst = 0;
   fVolName1rst        = std::string("");
   fMaterialNumberLast = 0;
   fTimeStart          = 0.;
   fPolarization       = G4ThreeVector(0.,0.,0.);
   fParentMomentumAtThisProduction = G4ThreeVector(0.,0.,0.); // for ppfx....A. Bashyal
   fPDGNucleus         = 0;
   fProcessName        = std::string("");
   fMaterialName1rst    = std::string("");

   fND=NumiDataInput::GetNumiDataInput();

}

NumiTrajectory::NumiTrajectory(const G4Track* aTrack)
{
   fParticleDefinition = aTrack->GetDefinition();
   fParticleName       = fParticleDefinition->GetParticleName();
   fPDGCharge          = fParticleDefinition->GetPDGCharge();
   fPDGEncoding        = fParticleDefinition->GetPDGEncoding();
   fTrackID            = aTrack->GetTrackID();
   fParentID           = aTrack->GetParentID();
   fPositionRecord     = new NumiTrajectoryPointContainer();
   fPositionRecord->push_back(new G4TrajectoryPoint(aTrack->GetPosition()));
   fMomentumRecord     = new NumiTrajectoryMomentumContainer();
   fMomentumRecord->push_back(aTrack->GetMomentum());
   fPreStepVolume      = new NumiTrajectoryVolumeName();
   fPreStepVolume->push_back(aTrack->GetVolume()->GetName());
   fMaterialName       = new NumiTrajectoryMaterialName();
   fMaterialName->push_back(aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName());
   fStepLength         = new DVec();
   fStepLength->push_back(aTrack->GetStepLength());
   fMomentum = aTrack->GetMomentum();
   fVertexPosition = aTrack->GetPosition();
   fParticleMass = aTrack->GetDefinition()->GetPDGMass();
   if (!aTrack->GetLogicalVolumeAtVertex()) {
//     std::cerr << " NumiTrajectory::NumiTrajectory No logical volume at vertex... ???    " << std::endl;
     fMaterialNumber1rst = 0;
     fMaterialName1rst = std::string("Unknown");
     fVolName1rst = std::string("Unknown");
    } else {
     fMaterialNumber1rst = aTrack->GetLogicalVolumeAtVertex()->GetMaterial()->GetIndex();
     fMaterialName1rst = aTrack->GetLogicalVolumeAtVertex()->GetMaterial()->GetName();
     fVolName1rst = aTrack->GetLogicalVolumeAtVertex()->GetName();
   }
   const G4VProcess* pCreatorProcess = aTrack->GetCreatorProcess();
   if (pCreatorProcess) {
       fProcessName = pCreatorProcess->GetProcessName();
   } else {
       fProcessName = "Primary";
   }
   fTimeStart = aTrack->GetGlobalTime();
   fPolarization = aTrack->GetPolarization();

   NumiTrackInformation* info = (NumiTrackInformation*)(aTrack->GetUserInformation());
   if (info!=0)
   {
     fDecayCode = info->GetDecayCode();
     fTgen = info->GetTgen();
     fNImpWt = info->GetNImpWt();
     fParentMomentumAtThisProduction = info->GetParentMomentumAtThisProduction();
   }
   else
  {
     fDecayCode = 0;
     fTgen      = 0;
     fNImpWt    = 1.;
     fParentMomentumAtThisProduction = G4ThreeVector(0.,0.,0.);
   }

   // Now the tricky bit: Dk2Nu package wants the process name and the nucleus
   // if (aTrack->GetCreatorProcess() != 0)     fProcessName = aTrack->GetCreatorProcess()->GetProcessName(); // Original 2014 code
   if (aTrack->GetCreatorProcess() != 0){ // attempt to get based on: https://geant4-forum.web.cern.ch/t/access-information-about-interaction-nucleus-of-material-in-step/5743
     fProcessName = aTrack->GetCreatorProcess()->GetProcessName();

     // if it derives from G4HadronicProcess then we can ask for the nucleus
     auto hadronicProcess =
       dynamic_cast<const G4HadronicProcess*>( aTrack->GetCreatorProcess() );

     if ( hadronicProcess ) {
       /* getting nucleus at this point is not necessarily getting the one associated with this track
	  this gets the nucleus associate with the last time this hadronicProcess was used.
	  So instead will store this info in NumiTrackInfo within NumiStackingAction and then retrieve here 
      
       G4Nucleus nucleus = *(hadronicProcess->GetTargetNucleus());
       fPDGNucleus = 1000000000 +
         nucleus.GetZ_asInt()*10000 + nucleus.GetA_asInt()*10;
       */
       NumiTrackInformation* trackInfo=(NumiTrackInformation*)(aTrack->GetUserInformation());  
       if (trackInfo) fPDGNucleus = trackInfo->GetPDGNucleus();
       
       // flag it as EM process with negative isotope code
       if ( fProcessName == "muonNuclear"     ||
            fProcessName == "electronNuclear" ||
            fProcessName == "positronNuclear"   ) {
         fPDGNucleus *= -1;
       }
     } else {
       // not derived from G4HadronicProcess

       if (        fProcessName == "Decay"         ) {
         fPDGNucleus =  1000000000;
       } else if ( fProcessName == "gamma"         ) {
         fPDGNucleus = -1000220220;
       } else if ( fProcessName == "PhotoNuclearXS") {
         fPDGNucleus = -1000220220;
       } else if ( fProcessName == "nCapture"      ) {
         fPDGNucleus =  1000210120;
       } else if ( fProcessName == "nKiller"       ) {
         fPDGNucleus = -1000210120;
       } else {
         std::cerr << " NumiTrajectory::NumiTrajectory process "
                   << fProcessName << " not a G4HadronicProcess ... "
                   << " and not special ProcessName, "
                   << "can't get nucleus  " << std::endl;
         fPDGNucleus = -1000000010;
       }
     }
   } else {
     // no CreatorProcess
     fProcessName = std::string("BeamParticle");
     fPDGNucleus = 0;
   }

   fND=NumiDataInput::GetNumiDataInput();
}

NumiTrajectory::NumiTrajectory(NumiTrajectory & right)
  : G4VTrajectory()
{
  fParticleName       = right.fParticleName;
  fParticleDefinition = right.fParticleDefinition;
  fPDGCharge          = right.fPDGCharge;
  fPDGEncoding        = right.fPDGEncoding;
  fTrackID            = right.fTrackID;
  fParentID           = right.fParentID;
  fMomentum           = right.fMomentum;
  fVertexPosition     = right.fVertexPosition;
  fParticleMass       = right.fParticleMass;
  fDecayCode          = right.fDecayCode;
  fTgen               = right.fTgen;
  fNImpWt             = right.fNImpWt;
  fPositionRecord     = new NumiTrajectoryPointContainer();
  fMomentumRecord     = new NumiTrajectoryMomentumContainer();
  fPreStepVolume      = new NumiTrajectoryVolumeName();
  fStepLength         = new DVec();
  fParentMomentumAtThisProduction = right.fParentMomentumAtThisProduction;
  fMaterialName       = new NumiTrajectoryMaterialName();

  for(size_t i=0;i<right.fPositionRecord->size();i++)
    {
      G4TrajectoryPoint* rightPoint = (G4TrajectoryPoint*)((*(right.fPositionRecord))[i]);
      fPositionRecord->push_back(new G4TrajectoryPoint(*rightPoint));
    }
  for(size_t i=0;i<right.fMomentumRecord->size();i++)
    {
      G4ThreeVector rightMomentum = (G4ThreeVector)((*(right.fMomentumRecord))[i]);
      fMomentumRecord->push_back(rightMomentum);
    }
  for(size_t i=0;i<right.fPreStepVolume->size();i++)
   {
     G4String rightPreStepVolume=(G4String)((*(right.fPreStepVolume))[i]);
     fPreStepVolume->push_back(rightPreStepVolume);
   }
  for(size_t i=0;i<right.fStepLength->size();i++)
   {
      G4double rightsteplength =(G4double)((*(right.fStepLength))[i]);
      fStepLength->push_back(rightsteplength);
   }
  for(size_t i=0;i<right.fMaterialName->size();i++)
    {
      G4String rightMaterialName = (G4String)((*(right.fMaterialName))[i]);
      fMaterialName->push_back(rightMaterialName);
    }


  fMaterialNumber1rst = right.fMaterialNumber1rst;
  fVolName1rst        = right.fVolName1rst;
  fMaterialNumberLast = right.fMaterialNumberLast;
  fTimeStart          = right.fTimeStart;
  fPolarization       = right.fPolarization;
  fPDGNucleus         = right.fPDGNucleus;
  fProcessName        = right.fProcessName;
  fMaterialName1rst   = right.fMaterialName1rst;
  fMaterialName       = right.fMaterialName;
  fMomentum           = right.fMomentum;
  fVertexPosition     = right.fVertexPosition;
  fParticleMass       = right.fParticleMass;
  fDecayCode          = right.fDecayCode;
  fTgen               = right.fTgen;
  fNImpWt             = right.fNImpWt;


  //
  fND=NumiDataInput::GetNumiDataInput();
}

NumiTrajectory::~NumiTrajectory()
{
  size_t i;
  for(i=0;i<fPositionRecord->size();i++){
    delete  (*fPositionRecord)[i];
  }
  fPositionRecord->clear();

  delete fPositionRecord;

  fMomentumRecord->clear();

  delete fMomentumRecord;
  fPreStepVolume->clear();
  delete fPreStepVolume;

  fStepLength->clear();
  delete fStepLength;

  fMaterialName->clear();
  delete fMaterialName;

}


void NumiTrajectory::ShowTrajectory() const
{
   G4cout << G4endl << "TrackID =" << fTrackID
        << " : ParentID=" << fParentID << G4endl;
   G4cout << "Particle name : " << fParticleName
        << "  Charge : " << fPDGCharge << G4endl;
   G4cout << "Original momentum : " <<
G4BestUnit(fMomentum,"Energy") << G4endl;
   G4cout << "Vertex : " << G4BestUnit(fVertexPosition,"Length") << G4endl;
   G4cout << "  Current trajectory has " << fPositionRecord->size()
        << " points." << G4endl;

   for( size_t i=0 ; i < fPositionRecord->size() ; i++){
       G4TrajectoryPoint* aTrajectoryPoint = (G4TrajectoryPoint*)((*fPositionRecord)[i]);
       G4cout << "Point[" << i << "]"
            << " Position= " << aTrajectoryPoint->GetPosition() << G4endl;
   }
}

void NumiTrajectory::ShowTrajectory(std::ostream& o) const
{
    G4VTrajectory::ShowTrajectory(o);
}

void NumiTrajectory::DrawTrajectory(G4int /*i_mode*/) const
{


   G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
   G4ThreeVector pos;

   G4Polyline pPolyline;
   for (size_t i = 0; i < fPositionRecord->size() ; i++) {
     G4TrajectoryPoint* aTrajectoryPoint = (G4TrajectoryPoint*)((*fPositionRecord)[i]);
     pos = aTrajectoryPoint->GetPosition();
     pPolyline.push_back( pos );
   }


   G4Colour color(1.0,1.0,1.0); //white

#ifdef SCHEME_LBNE
   if(fPDGCharge ==  0.)   color = G4Colour::White();
   if(fPDGCharge ==  1.0)  color = G4Colour(1.0, 0.5, 0.0); //orange
   if(fPDGCharge == -1.0)  color = G4Colour(0.5, 0.0, 0.5); //purple


   if (fParticleDefinition==G4Proton::ProtonDefinition()) color = G4Color::Green();
   if (fParticleDefinition==G4Neutron::NeutronDefinition()) color = G4Color::Grey();
   if (fParticleDefinition==G4PionMinus::PionMinusDefinition()) color = G4Color::Cyan();
   if (fParticleDefinition==G4PionPlus::PionPlusDefinition()) color = G4Color::Blue();
   if (fParticleDefinition==G4KaonMinus::KaonMinusDefinition()) color = G4Color::Yellow();
   if (fParticleDefinition==G4KaonPlus::KaonPlusDefinition()) color = G4Color::Yellow();
   if (fParticleDefinition==G4MuonMinus::MuonMinusDefinition()) color = G4Color::Red();
   if (fParticleDefinition==G4MuonPlus::MuonPlusDefinition()) color = G4Color::Red();

   if ((fParticleDefinition == G4NeutrinoE::NeutrinoEDefinition())  ||
       (fParticleDefinition == G4NeutrinoMu::NeutrinoMuDefinition()) ||
       (fParticleDefinition == G4NeutrinoTau::NeutrinoTauDefinition()) ||
       (fParticleDefinition == G4AntiNeutrinoE::AntiNeutrinoEDefinition()) ||
       (fParticleDefinition == G4AntiNeutrinoMu::AntiNeutrinoMuDefinition()) ||
       (fParticleDefinition == G4AntiNeutrinoTau::AntiNeutrinoTauDefinition()))
      color = G4Color::Magenta();

   G4VisAttributes attribs(color);
   pPolyline.SetVisAttributes(attribs);
   if(pVVisManager) pVVisManager->Draw(pPolyline);
#else
   G4Colour colour(0.2,0.2,0.2);
   if(fParticleDefinition==G4Gamma::GammaDefinition())
      colour = G4Colour(0.,0.,1.);
   else if(fParticleDefinition==G4Electron::ElectronDefinition()
         ||fParticleDefinition==G4Positron::PositronDefinition())
      colour = G4Colour(1.,1.,0.);
   else if(fParticleDefinition==G4MuonMinus::MuonMinusDefinition()
         ||fParticleDefinition==G4MuonPlus::MuonPlusDefinition())
      colour = G4Colour(0.,1.,0.);
   else if(fParticleDefinition->GetParticleType()=="meson")
   {
      if(fPDGCharge!=0.)
         colour = G4Colour(1.,0.,0.);
      else
         colour = G4Colour(0.5,0.,0.);
   }
   else if(fParticleDefinition->GetParticleType()=="baryon")
   {
      if(fPDGCharge!=0.)
         colour = G4Colour(0.,1.,1.);
      else
         colour = G4Colour(0.,0.5,0.5);
   }

   //G4VisAttributes attribs(colour);

   //draw only protons,pi+ and pi-
   G4VisAttributes attribs;
   if (fParticleDefinition==G4Proton::ProtonDefinition()) {
     colour=G4Colour(0.,0.,1.);
     attribs=G4VisAttributes(colour);
     pPolyline.SetVisAttributes(attribs);
     if(pVVisManager) pVVisManager->Draw(pPolyline);
   }
   if (fParticleDefinition==G4PionMinus::PionMinusDefinition()) {
     colour=G4Colour(1.,0.,0.);
     attribs=G4VisAttributes(colour);
     pPolyline.SetVisAttributes(attribs);
     if(pVVisManager) pVVisManager->Draw(pPolyline);
   }
   if (fParticleDefinition==G4PionPlus::PionPlusDefinition()) {
     colour=G4Colour(0.,1.,0.);
     attribs=G4VisAttributes(colour);
     pPolyline.SetVisAttributes(attribs);
     if(pVVisManager) pVVisManager->Draw(pPolyline);
   }

   if(fND->useMuonBeam)
   {
      if (fParticleDefinition==G4MuonMinus::MuonMinusDefinition()
          || fParticleDefinition==G4MuonPlus::MuonPlusDefinition())
      {
         //colour=G4Colour(0.,1.,0.5);
         colour=G4Colour(1.,1.,1.);
         attribs=G4VisAttributes(colour);
         pPolyline.SetVisAttributes(&attribs);
         if(pVVisManager) pVVisManager->Draw(pPolyline);
      }
      if (fParticleDefinition==G4Gamma::GammaDefinition())
      {
         colour=G4Colour(0.,1.,0.);
         attribs=G4VisAttributes(colour);
         pPolyline.SetVisAttributes(&attribs);
         if(pVVisManager) pVVisManager->Draw(pPolyline);
      }
      if (fParticleDefinition==G4Electron::ElectronDefinition()
          || fParticleDefinition==G4Positron::PositronDefinition())
      {
         colour=G4Colour(1.,0.,1.);
         attribs=G4VisAttributes(colour);
         pPolyline.SetVisAttributes(&attribs);
         if(pVVisManager) pVVisManager->Draw(pPolyline);
      }
   }

   //pPolyline.SetVisAttributes(attribs);
   //if(pVVisManager) pVVisManager->Draw(pPolyline);
#endif
}

//-------------------------------------------------------------------------
const std::map<G4String,G4AttDef>* NumiTrajectory::GetAttDefs() const
{
    G4bool isNew;
    std::map<G4String,G4AttDef>* store
       = G4AttDefStore::GetInstance("NumiTrajectory",isNew);
    if (isNew)
    {

       G4String ID("TrkID");
       (*store)[ID] = G4AttDef(ID,"Track ID","Physics","","G4int");

       G4String PID("ParTrkID");
       (*store)[PID] = G4AttDef(PID,"Parent ID","Physics","","G4int");

       G4String PN("Name");
       (*store)[PN] = G4AttDef(PN,"Particle Name","Physics","","G4String");

       G4String Ch("Q");
       (*store)[Ch] = G4AttDef(Ch,"Charge","Physics","e+","G4double");

       G4String PDG("PDG");
       (*store)[PDG] = G4AttDef(PDG,"PDG Encoding","Physics","","G4int");

       G4String IMom("IMom");
       (*store)[IMom] = G4AttDef(IMom, "Momentum at start of trajectory",
                                 "Physics","G4BestUnit","G4ThreeVector");

       G4String IMag("IMomMag");
       (*store)[IMag] =
         G4AttDef(IMag, "Magnitude of momentum at start of trajectory",
                  "Physics","G4BestUnit","G4double");

       G4String NTP("NPts");
       (*store)[NTP] = G4AttDef(NTP,"No. of points","Physics","","G4int");

    }
    return store;
}
//-------------------------------------------------------------------------
std::vector<G4AttValue>* NumiTrajectory::CreateAttValues() const
{
   std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

   values->push_back
      (G4AttValue("TrkID",G4UIcommand::ConvertToString(fTrackID),""));

   values->push_back
      (G4AttValue("ParTrkID",G4UIcommand::ConvertToString(fParentID),""));

   values->push_back(G4AttValue("Name",fParticleName,""));

   values->push_back
      (G4AttValue("Q",G4UIcommand::ConvertToString(fPDGCharge),""));

   values->push_back
      (G4AttValue("PDG",G4UIcommand::ConvertToString(fPDGEncoding),""));

   values->push_back
      (G4AttValue("IMom",G4BestUnit(GetInitialMomentum(),"Energy"),""));

   values->push_back
      (G4AttValue("IMomMag",G4BestUnit((GetInitialMomentum()).mag(),"Energy"),""));

   values->push_back
      (G4AttValue("NPts",G4UIcommand::ConvertToString(GetPointEntries()),""));

   return values;
}


//-------------------------------------------------------------------------

void NumiTrajectory::AppendStep(const G4Step* aStep)
{
  fPositionRecord
    ->push_back(new G4TrajectoryPoint(aStep->GetPostStepPoint()
                                      ->GetPosition() ));
   fMomentumRecord->push_back(aStep->GetPostStepPoint()->GetMomentum());

   G4Track* aTrack=aStep->GetTrack();
   NumiTrackInformation* info=(NumiTrackInformation*)(aTrack->GetUserInformation());
   if (info!=0) {
     fDecayCode=info->GetDecayCode();
     fParentMomentumAtThisProduction = info->GetParentMomentumAtThisProduction();
     fTgen=info->GetTgen();
   }
   else
   {
     fDecayCode=-1;
     fParentMomentumAtThisProduction = G4ThreeVector(0.,0.,0.);
   }


   G4StepPoint * steppoint=aStep->GetPreStepPoint();
   G4String PreVolumeName=steppoint->GetPhysicalVolume()->GetName();
   G4String MaterialName =
     steppoint->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetName();
   fPreStepVolume->push_back(PreVolumeName);
   fStepLength->push_back(aStep->GetStepLength());
   fMaterialName->push_back(MaterialName);
}

G4ParticleDefinition* NumiTrajectory::GetParticleDefinition()
{
   return (G4ParticleTable::GetParticleTable()->FindParticle(fParticleName));
}

void NumiTrajectory::MergeTrajectory(G4VTrajectory* secondTrajectory)
{
  if(!secondTrajectory) return;

  NumiTrajectory* seco = (NumiTrajectory*)secondTrajectory;
  G4int ent = seco->GetPointEntries();
  for(int i=1;i<ent;i++) // initial point of the second trajectory should not be merged
  {
    fPositionRecord->push_back((*(seco->fPositionRecord))[i]);
  }
  delete (*seco->fPositionRecord)[0];
  seco->fPositionRecord->clear();

}
