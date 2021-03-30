
#ifndef G4MinervaProtonIXS_h
#define G4MinervaProtonIXS_h

#include "globals.hh"
#include "G4Proton.hh"
#include "G4DynamicParticle.hh"
#include "G4VCrossSectionDataSet.hh"


#include "CLHEP/Units/PhysicalConstants.h"

class G4MinervaProtonIXS : public G4VCrossSectionDataSet
{
   public:
  G4MinervaProtonIXS();
   virtual
   G4bool IsApplicable(const G4DynamicParticle* aPart, const G4Element* aEle)
   {
     G4bool result = false;
     if(( aPart->GetDefinition()==G4Proton::Proton()) &&
        ( aPart->GetKineticEnergy()<1*CLHEP::TeV) ) result = true;
     if(aEle->GetZ()<3) result = false;
     return result;
   }

   G4bool IsZAApplicable(const G4DynamicParticle* aParticle,
                         G4double ZZ, G4double /*AA*/)
   {
     G4bool result = false;
     if (( aParticle->GetDefinition() == G4Proton::Proton()) &&
         ( aParticle->GetKineticEnergy() < 1*CLHEP::TeV) ) result = true;
     if (ZZ < 3) result = false;
     return result;
   }

   virtual
   G4double GetCrossSection(const G4DynamicParticle*, 
                            const G4Element*, G4double aTemperature);
   
   G4double GetIsoZACrossSection(const G4DynamicParticle* aParticle, 
                                 G4double ZZ, G4double AA, 
                                 G4double /*aTemperature*/)
   {
     return GetCrossSection(aParticle->GetKineticEnergy(), AA, ZZ);
   }
 

   G4double GetCrossSection(G4double anEnergy, G4double anA, G4double aZ);

   virtual
   void BuildPhysicsTable(const G4ParticleDefinition&)
   {}

   virtual
   void DumpPhysicsTable(const G4ParticleDefinition&) 
   {G4cout << "G4MinervaProtonIXS: uses formula"<<G4endl;}
 
  //Set scale value: 
  inline void SetScale(G4double value)
  {
    scaleVal = value;
  }

  inline G4double GetScale()
  {
    return scaleVal;
  }
 
  G4double scaleVal;

 
};

#endif
