//
// NumiTrajectory.hh
//

#ifndef NumiTrajectory_h
#define NumiTrajectory_h 1

#include "G4VTrajectory.hh"
#include "G4Allocator.hh"
#include <stdlib.h>
#include "G4ThreeVector.hh"
#include "G4ios.hh"     
#include <vector>
#include <iostream>
#include "globals.hh" 
#include "G4ParticleDefinition.hh" 
#include "G4TrajectoryPoint.hh"   
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4UserSteppingAction.hh"

class NumiTrajectory;
class G4Polyline;

typedef std::vector<G4VTrajectoryPoint*> NumiTrajectoryPointContainer;
typedef std::vector<G4ThreeVector> NumiTrajectoryMomentumContainer;
typedef std::vector<G4String> NumiTrajectoryVolumeName;

class NumiTrajectory : public G4VTrajectory
{
 public:
   NumiTrajectory();
   NumiTrajectory(const G4Track* aTrack);
   NumiTrajectory(NumiTrajectory &);
   virtual ~NumiTrajectory();

   inline void* operator new(size_t);
   inline void  operator delete(void*);
   inline int operator == (const NumiTrajectory& right) const
   {return (this==&right);} 

   inline G4int GetTrackID() const
   { return fTrackID; }
   inline G4int GetParentID() const
   { return fParentID; }
   inline G4String GetParticleName() const
   { return ParticleName; }
   inline G4double GetCharge() const
   { return PDGCharge; }
   inline G4double GetMass() const
   { return ParticleMass; }
   inline G4int GetPDGEncoding() const
   { return PDGEncoding; }
   inline const G4ThreeVector& GetVertexPosition() const
   { return vertexPosition; }
   virtual int GetPointEntries() const
   { return positionRecord->size(); }
   virtual G4VTrajectoryPoint* GetPoint(G4int i) const 
   { return (*positionRecord)[i]; }
   virtual G4ThreeVector GetMomentum(G4int i) const
   { return (*momentumRecord)[i]; }
   virtual G4String GetPreStepVolumeName(G4int i) const
   { return (*PreStepVolume)[i]; }
   inline G4ThreeVector GetInitialMomentum() const 
   { return momentum; }                            
   virtual G4int Gettgen() const
   { return tgen;}
   inline G4int GetDecayCode() const
   {return decaycode;}
   virtual G4double GetNImpWt() const
   {return nimpwt;}
  
   virtual void ShowTrajectory() const;
   virtual void ShowTrajectory(std::ostream& o) const;
   virtual void DrawTrajectory(G4int i_mode=0) const;
   virtual void AppendStep(const G4Step* aStep);
   virtual void MergeTrajectory(G4VTrajectory* secondTrajectory);

   G4ParticleDefinition* GetParticleDefinition();

 private:
   NumiTrajectoryMomentumContainer* momentumRecord;
   NumiTrajectoryPointContainer*    positionRecord;
   G4int                            decaycode;
   G4int                            eventno;
   G4int                            tgen;
   G4double                         nimpwt;
   G4int                            fTrackID;
   G4int                            fParentID;
   G4ParticleDefinition*            fpParticleDefinition;
   G4String                         ParticleName;
   G4double                         PDGCharge;
   G4int                            PDGEncoding;
   G4ThreeVector                    momentum;
   G4ThreeVector                    vertexPosition;
   G4ThreeVector                    InitialPolarization;
   G4ThreeVector                    FinalPolarization;
   G4double                         ParticleMass;
   NumiTrajectoryVolumeName*        PreStepVolume;
};

extern G4Allocator<NumiTrajectory> myTrajectoryAllocator;

inline void* NumiTrajectory::operator new(size_t)
{
  void* aTrajectory;
  aTrajectory = (void*)myTrajectoryAllocator.MallocSingle();
  return aTrajectory;
}

inline void NumiTrajectory::operator delete(void* aTrajectory)
{
  myTrajectoryAllocator.FreeSingle((NumiTrajectory*)aTrajectory);
}
#endif



