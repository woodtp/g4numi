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
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"
#include "NumiTrackInformation.hh"

G4Allocator<NumiTrajectory> myTrajectoryAllocator;

NumiTrajectory::NumiTrajectory()
{
   fpParticleDefinition = 0;
   ParticleName = "";
   PDGCharge = 0;
   PDGEncoding = 0;
   fTrackID = 0;
   fParentID = 0;
   positionRecord = 0;
   momentum = G4ThreeVector(0.,0.,0.);
   momentumRecord = 0;
   vertexPosition = G4ThreeVector(0.,0.,0.);
   ParticleMass = 0.;
   decaycode=0;
   tgen=0;
   PreStepVolume=0;

}

NumiTrajectory::NumiTrajectory(const G4Track* aTrack)
{
   fpParticleDefinition = aTrack->GetDefinition();
   ParticleName = fpParticleDefinition->GetParticleName();
   PDGCharge = fpParticleDefinition->GetPDGCharge();
   PDGEncoding = fpParticleDefinition->GetPDGEncoding();
   fTrackID = aTrack->GetTrackID();
   fParentID = aTrack->GetParentID();
   positionRecord = new NumiTrajectoryPointContainer();
   positionRecord->push_back(new G4TrajectoryPoint(aTrack->GetPosition()));
   momentumRecord = new NumiTrajectoryMomentumContainer();
   momentumRecord->push_back(aTrack->GetMomentum());
   
   PreStepVolume = new NumiTrajectoryVolumeName();  
   PreStepVolume->push_back(aTrack->GetVolume()->GetName());
   
     
   momentum = aTrack->GetMomentum();
   vertexPosition = aTrack->GetPosition();
   ParticleMass = aTrack->GetDefinition()->GetPDGMass();

   NumiTrackInformation* info=(NumiTrackInformation*)(aTrack->GetUserInformation());
   if (info!=0) {
     decaycode=info->GetDecayCode();
     tgen = info->Gettgen();
     nimpwt = info->GetNImpWt();}
   else { 
     decaycode=0;
     tgen=0; 
     nimpwt=1.;
   }
}

NumiTrajectory::NumiTrajectory(NumiTrajectory & right)
  : G4VTrajectory()
{
  ParticleName = right.ParticleName;
  fpParticleDefinition = right.fpParticleDefinition;
  PDGCharge = right.PDGCharge;
  PDGEncoding = right.PDGEncoding;
  fTrackID = right.fTrackID;
  fParentID = right.fParentID;
  positionRecord = new NumiTrajectoryPointContainer();
  for(size_t i=0;i<right.positionRecord->size();i++)
    {
      G4TrajectoryPoint* rightPoint = (G4TrajectoryPoint*)((*(right.positionRecord))[i]);
      positionRecord->push_back(new G4TrajectoryPoint(*rightPoint));
    }
  for(size_t i=0;i<right.momentumRecord->size();i++)
    {
      G4ThreeVector rightMomentum = (G4ThreeVector)((*(right.momentumRecord))[i]);
      momentumRecord->push_back(rightMomentum);
    }
  for(size_t i=0;i<right.PreStepVolume->size();i++)
   {
     G4String rightPreStepVolume=(G4String)((*(right.PreStepVolume))[i]);
     PreStepVolume->push_back(rightPreStepVolume);
   }
  momentum = right.momentum;
  vertexPosition = right.vertexPosition;
  ParticleMass = right.ParticleMass;
  decaycode = right.decaycode;
  tgen = right.tgen;
  nimpwt = right.nimpwt;
}

NumiTrajectory::~NumiTrajectory()
{
  size_t i;
  for(i=0;i<positionRecord->size();i++){
    delete  (*positionRecord)[i];
  }
  positionRecord->clear();

  delete positionRecord;

  momentumRecord->clear();
  
  delete momentumRecord;
  PreStepVolume->clear();
  delete PreStepVolume;
 
}


void NumiTrajectory::ShowTrajectory() const
{
   G4cout << G4endl << "TrackID =" << fTrackID 
        << " : ParentID=" << fParentID << G4endl;
   G4cout << "Particle name : " << ParticleName 
        << "  Charge : " << PDGCharge << G4endl;
   G4cout << "Original momentum : " <<
G4BestUnit(momentum,"Energy") << G4endl;
   G4cout << "Vertex : " << G4BestUnit(vertexPosition,"Length") << G4endl;
   G4cout << "  Current trajectory has " << positionRecord->size() 
        << " points." << G4endl;

   for( size_t i=0 ; i < positionRecord->size() ; i++){
       G4TrajectoryPoint* aTrajectoryPoint = (G4TrajectoryPoint*)((*positionRecord)[i]);
       G4cout << "Point[" << i << "]" 
            << " Position= " << aTrajectoryPoint->GetPosition() << G4endl;
   }
}

void NumiTrajectory::ShowTrajectory(std::ostream& o) const
{
    G4VTrajectory::ShowTrajectory(o);
}

void NumiTrajectory::DrawTrajectory(G4int i_mode) const
{

   G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
   G4ThreeVector pos;

   G4Polyline pPolyline;
   for (size_t i = 0; i < positionRecord->size() ; i++) {
     G4TrajectoryPoint* aTrajectoryPoint = (G4TrajectoryPoint*)((*positionRecord)[i]);
     pos = aTrajectoryPoint->GetPosition();
     pPolyline.push_back( pos );
   }

   G4Colour colour(0.2,0.2,0.2);
   if(fpParticleDefinition==G4Gamma::GammaDefinition())
      colour = G4Colour(0.,0.,1.);
   else if(fpParticleDefinition==G4Electron::ElectronDefinition()
         ||fpParticleDefinition==G4Positron::PositronDefinition())
      colour = G4Colour(1.,1.,0.);
   else if(fpParticleDefinition==G4MuonMinus::MuonMinusDefinition()
         ||fpParticleDefinition==G4MuonPlus::MuonPlusDefinition())
      colour = G4Colour(0.,1.,0.);
   else if(fpParticleDefinition->GetParticleType()=="meson")
   {
      if(PDGCharge!=0.)
         colour = G4Colour(1.,0.,0.);
      else
         colour = G4Colour(0.5,0.,0.);
   }
   else if(fpParticleDefinition->GetParticleType()=="baryon")
   {
      if(PDGCharge!=0.)
         colour = G4Colour(0.,1.,1.);
      else
         colour = G4Colour(0.,0.5,0.5);
   }

   //G4VisAttributes attribs(colour);

   //draw only protons,pi+ and pi-
   G4VisAttributes attribs;
   if (fpParticleDefinition==G4Proton::ProtonDefinition()) {
     colour=G4Colour(0.,0.,1.);
     attribs=G4VisAttributes(colour);
     pPolyline.SetVisAttributes(attribs);
     if(pVVisManager) pVVisManager->Draw(pPolyline);
   }
   if (fpParticleDefinition==G4PionMinus::PionMinusDefinition()) {
     colour=G4Colour(1.,0.,0.);
     attribs=G4VisAttributes(colour);
     pPolyline.SetVisAttributes(attribs);
     if(pVVisManager) pVVisManager->Draw(pPolyline);
   }
   if (fpParticleDefinition==G4PionPlus::PionPlusDefinition()) {
     colour=G4Colour(0.,1.,0.);
     attribs=G4VisAttributes(colour);
     pPolyline.SetVisAttributes(attribs);
     if(pVVisManager) pVVisManager->Draw(pPolyline);
   }

   //pPolyline.SetVisAttributes(attribs);
   //if(pVVisManager) pVVisManager->Draw(pPolyline);
}

void NumiTrajectory::AppendStep(const G4Step* aStep)
{
   positionRecord->push_back( new G4TrajectoryPoint(aStep->GetPostStepPoint()->
GetPosition() ));
   momentumRecord->push_back(aStep->GetPostStepPoint()->GetMomentum());

   G4Track* aTrack=aStep->GetTrack();
   NumiTrackInformation* info=(NumiTrackInformation*)(aTrack->GetUserInformation());
   if (info!=0) {
     decaycode=info->GetDecayCode();
     tgen=info->Gettgen();
   }
   else decaycode=-1;
   
   G4StepPoint * steppoint=aStep->GetPreStepPoint(); 
   G4String PreVolumeName=steppoint->GetPhysicalVolume()->GetName(); 
   PreStepVolume->push_back(PreVolumeName); 
   
}
  
G4ParticleDefinition* NumiTrajectory::GetParticleDefinition()
{
   return (G4ParticleTable::GetParticleTable()->FindParticle(ParticleName));
}

void NumiTrajectory::MergeTrajectory(G4VTrajectory* secondTrajectory)
{
  if(!secondTrajectory) return;

  NumiTrajectory* seco = (NumiTrajectory*)secondTrajectory;
  G4int ent = seco->GetPointEntries();
  for(int i=1;i<ent;i++) // initial point of the second trajectory should not be merged
  {
    positionRecord->push_back((*(seco->positionRecord))[i]);
  }
  delete (*seco->positionRecord)[0];
  seco->positionRecord->clear();

}

