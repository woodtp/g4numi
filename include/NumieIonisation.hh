//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: NumieIonisation.hh,v 1.1.2.1 2010/08/19 19:50:54 minervacvs Exp $
// GEANT4 tag $Name:  $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     NumieIonisation
//
// Author:        Laszlo Urban
//
// Creation date: 30.05.1997
//
// Modifications:
//
// corrected by L.Urban on 24/09/97
// corrected by L.Urban on 13/01/98
// bugs fixed by L.Urban on 02/02/99
// 10/02/00 modifications , new e.m. structure, L.Urban
// 10-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 14-08-01 new function ComputeRestrictedMeandEdx() + 'cleanup' (mma)
// 19-09-01 come back to previous process name "hIoni"
// 29-10-01 all static functions no more inlined
// 10-05-02 V.Ivanchenko update to new design
// 09-12-02 V.Ivanchenko remove warning
// 26-12-02 Secondary production moved to derived classes (VI)
// 24-01-03 Make models region aware (V.Ivanchenko)
// 05-02-03 Fix compilation warnings (V.Ivanchenko)
// 13-02-03 SubCutoff regime is assigned to a region (V.Ivanchenko)
// 23-05-03 Add fluctuation model as a member function (V.Ivanchenko)
// 03-06-03 Add SetIntegral method to choose fluctuation model (V.Ivanchenko)
// 03-06-03 Fix initialisation problem for STD ionisation (V.Ivanchenko)
// 08-08-03 STD substitute standard  (V.Ivanchenko)
// 12-11-03 G4EnergyLossSTD -> G4EnergyLossProcess (V.Ivanchenko)
// 21-01-04 Migrade to G4ParticleChangeForLoss (V.Ivanchenko)
// 17-08-04 Rename the process "Mu" -> "mu" (V.Ivanchenko)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
//
// Class Description:
//
// This class manages the ionisation process for muons.
// it inherites from G4VContinuousDiscreteProcess via G4VEnergyLossProcess.
//

// -------------------------------------------------------------------
//

#ifndef NumieIonisation_h
#define NumieIonisation_h 1

#include "G4VEnergyLossProcess.hh"
#include "G4VEnergyLoss.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "globals.hh"
#include "G4VEmModel.hh"

#include "G4VContinuousDiscreteProcess.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Track.hh"
#include "G4EmModelManager.hh"
#include "G4UnitsTable.hh"
#include "G4VParticleChange.hh"
#include "G4EmTableType.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"

#include "G4RegionStore.hh"


class G4Material;
class G4VEmFluctuationModel;
class G4ParticleChangeForLoss;

class NumieIonisation : public G4VEnergyLoss
{

public:

   NumieIonisation(const G4String& name = "numieIoni",
                    G4ProcessType aType = fElectromagnetic );

  virtual ~NumieIonisation();

  G4bool IsApplicable(const G4ParticleDefinition& p);

  G4double MinPrimaryEnergy(const G4ParticleDefinition* p,
			    const G4Material*, G4double cut);

  // Print out of the class parameters
  void PrintInfo();

   void PreparePhysicsTable(const G4ParticleDefinition& part);

   inline G4double GetContinuousStepLimit(const G4Track& track,
                                   G4double previousStepSize,
                                   G4double currentMinimumStep,
                                   G4double& currentSafety);
   
   inline G4double GetMeanFreePath(const G4Track& track,
                                   G4double previousStepSize,
                                   G4ForceCondition* condition);

   inline G4double AlongStepGetPhysicalInteractionLength(const G4Track&,
                                                         G4double  previousStepSize,
                                                         G4double  currentMinimumStep,
                                                         G4double& currentSafety,
                                                         G4GPILSelection* selection);

   inline G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                        G4double   previousStepSize,
                                                        G4ForceCondition* condition);

   G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&);
   
   G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

   void SetStepFunction(G4double v1, G4double v2);

   G4PhysicsTable* BuildDEDXTable();

   // Add EM model coupled with fluctuation model for the region
   inline void AddEmModel(G4int, G4VEmModel*, G4VEmFluctuationModel* fluc = 0,
                          const G4Region* region = 0);

   inline const G4String& GetProcessName() const;
   



protected:

  virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
                                           const G4ParticleDefinition*);

private:

   inline G4VEmModel* SelectModel(G4double kinEnergy);

   inline G4VEmModel* SelectModelForMaterial(G4double kinEnergy,
                                             size_t& idx) const;

   inline void InitialiseStep(const G4Track&);

   inline void DefineMaterial(const G4MaterialCutsCouple* couple);

   inline void InitialiseMassCharge(const G4Track&);

   inline G4double GetScaledRangeForScaledEnergy(G4double e);

   inline G4double GetDEDXForScaledEnergy(G4double e);

   inline G4double ScaledKinEnergyForLoss(G4double r);

private:

   // hide assignment operator
   NumieIonisation & operator=(const NumieIonisation &right);
   NumieIonisation(const NumieIonisation&);

   G4EmModelManager*                     modelManager;

   G4String    fName;
   
   G4double    mass;
   G4double    ratio;

   G4double lowestKinEnergy;
   G4double minKinEnergy;
   G4double maxKinEnergy;
   G4double maxKinEnergyCSDA;
   
   G4double massRatio;
   G4double reduceFactor;
   //G4double chargeSquare;
   G4double chargeSqRatio;
   
   G4double preStepLambda;
   G4double fRange;
   G4double preStepKinEnergy;
   G4double preStepScaledEnergy;
   
   G4double minSubRange;
   G4double dRoverRange;
   G4double finalRange;
   G4double lambdaFactor;
   G4double mfpKinEnergy;

   const G4ParticleDefinition* theElectron;
   
  G4bool isElectron;

   const G4ParticleDefinition* theParticle;
   const G4ParticleDefinition* theSecondaryParticle;
   const G4ParticleDefinition* theBaseParticle;
   G4VEmModel*                 EmModel;
   G4VEmFluctuationModel*      flucModel;
   
   G4ParticleChangeForLoss*    fParticleChange;
   
   G4bool                      isInitialised;
   G4bool                      isIonisation;

   G4int                       verboseLevel;
   
   G4PhysicsTable*             theDEDXTable;
   G4PhysicsTable*             theRangeTableForLoss;
   G4PhysicsTable*             theInverseRangeTable;

   const G4DataVector*         theCuts;

   const G4Material*           currentMaterial;
   const G4MaterialCutsCouple* currentCouple;
   size_t                      currentMaterialIndex;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool NumieIonisation::IsApplicable(const G4ParticleDefinition& p)
{
   return (&p == G4Electron::Electron() || &p == G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double NumieIonisation::MinPrimaryEnergy(const G4ParticleDefinition*,
						 const G4Material*,
						 G4double cut)
{
   G4double x = cut;
   if(isElectron) x += cut;
   return x;
}

inline void NumieIonisation::InitialiseStep(const G4Track& track)
{
  InitialiseMassCharge(track);
  preStepKinEnergy = track.GetKineticEnergy();
  preStepScaledEnergy = preStepKinEnergy*massRatio;
  DefineMaterial(track.GetMaterialCutsCouple());
//  if (theNumberOfInteractionLengthLeft < 0.0) mfpKinEnergy = DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void NumieIonisation::InitialiseMassCharge(const G4Track&)
{
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void NumieIonisation::DefineMaterial(const G4MaterialCutsCouple* couple)
{
  if(couple != currentCouple)
  {
    currentCouple   = couple;
    currentMaterial = couple->GetMaterial();
    currentMaterialIndex = couple->GetIndex();
    mfpKinEnergy = DBL_MAX;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* NumieIonisation::SelectModel(G4double kinEnergy)
{
  return modelManager->SelectModel(kinEnergy, currentMaterialIndex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* NumieIonisation::SelectModelForMaterial(G4double kinEnergy,
                                                            size_t& idx) const
{
  return modelManager->SelectModel(kinEnergy, idx);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
inline G4double NumieIonisation::GetContinuousStepLimit(const G4Track& ,
                                                  G4double,
                                                  G4double,
                                                  G4double& )
{

  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
inline G4double  NumieIonisation::GetMeanFreePath(const G4Track&,
                                                   G4double,
                                                   G4ForceCondition* condition)
{
   // return infinity so that it does nothing.
   *condition = NotForced;
   return DBL_MAX;
   
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double NumieIonisation::AlongStepGetPhysicalInteractionLength(const G4Track&,
                                                                        G4double,
                                                                        G4double,
                                                                        G4double&,
                                                                        G4GPILSelection*)
{

   return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double NumieIonisation::GetScaledRangeForScaledEnergy(G4double e)
{
  G4bool b;
  G4double x = ((*theRangeTableForLoss)[currentMaterialIndex])->GetValue(e, b);
  if(e < minKinEnergy) x *= std::sqrt(e/minKinEnergy);
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double NumieIonisation::GetDEDXForScaledEnergy(G4double e)
{
  G4bool b;
  G4double x = 
    ((*theDEDXTable)[currentMaterialIndex]->GetValue(e, b))*chargeSqRatio;
  if(e < minKinEnergy) x *= std::sqrt(e/minKinEnergy);
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double NumieIonisation::ScaledKinEnergyForLoss(G4double r)
{
  G4PhysicsVector* v = (*theInverseRangeTable)[currentMaterialIndex];
  G4double rmin = v->GetLowEdgeEnergy(0);
  G4double e = 0.0; 
  if(r >= rmin)
  {
    G4bool b;
    e = v->GetValue(r, b);
  }
  else if(r > 0.0)
  {
    G4double x = r/rmin;
    e = minKinEnergy*x*x;
  }
  return e;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double NumieIonisation::PostStepGetPhysicalInteractionLength(
                             const G4Track&,
                             G4double,
                             G4ForceCondition* condition)
{

   
  // condition is set to "Not Forced"
  *condition = NotForced;
  /*
  G4double x = DBL_MAX;
  if(previousStepSize <= DBL_MIN) theNumberOfInteractionLengthLeft = -1.0;
  InitialiseStep(track);

  if(preStepScaledEnergy < mfpKinEnergy) {
    if (integral) ComputeLambdaForScaledEnergy(preStepScaledEnergy);
    else  preStepLambda = GetLambdaForScaledEnergy(preStepScaledEnergy);
    if(preStepLambda <= DBL_MIN) mfpKinEnergy = 0.0;
  }

  if(preStepLambda > DBL_MIN) {

    if (theNumberOfInteractionLengthLeft < 0.0) {
      // beggining of tracking (or just after DoIt of this process)
      ResetNumberOfInteractionLengthLeft();
    } else if(previousStepSize > DBL_MIN) {
      // subtract NumberOfInteractionLengthLeft
      SubtractNumberOfInteractionLengthLeft(previousStepSize);
      if(theNumberOfInteractionLengthLeft<0.)
	theNumberOfInteractionLengthLeft=perMillion;
    }

    // get mean free path
    currentInteractionLength = 1.0/preStepLambda;
    x = theNumberOfInteractionLengthLeft * currentInteractionLength;
#ifdef G4VERBOSE
    if (verboseLevel>2){
      G4cout << "G4VEnergyLossProcess::PostStepGetPhysicalInteractionLength ";
      G4cout << "[ " << GetProcessName() << "]" << G4endl; 
      G4cout << " for " << particle->GetParticleName() 
             << " in Material  " <<  currentMaterial->GetName()
	     << " Ekin(MeV)= " << preStepKinEnergy/MeV 
	     <<G4endl;
      G4cout << "MeanFreePath = " << currentInteractionLength/cm << "[cm]" 
	     << "InteractionLength= " << x/cm <<"[cm] " <<G4endl;
    }
#endif
  }
  return x;
  */


  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void NumieIonisation::AddEmModel(G4int order, G4VEmModel* p, 
                                  G4VEmFluctuationModel* fluc,
                                  const G4Region* region)
{
   
  modelManager->AddEmModel(order, p, fluc, region); 
  if(p) p->SetParticleChange(fParticleChange, fluc);
 
  
}

const G4String& NumieIonisation::GetProcessName() const
{
   return fName;
}



#endif
