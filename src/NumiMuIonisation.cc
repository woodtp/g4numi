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
// $Id: NumiMuIonisation.cc,v 1.1.2.1 2010/08/19 19:50:54 minervacvs Exp $
// GEANT4 tag $Name:  $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     NumiMuIonisation
//
// Author:        Laszlo Urban
//
// Creation date: 30.09.1997
//
// Modifications:
//
// 08-04-98 remove 'tracking cut' of the ionizing particle (mma)
// 26-10-98 new stuff from R.Kokoulin + cleanup , L.Urban
// 10-02-00 modifications , new e.m. structure, L.Urban
// 23-03-01 R.Kokoulin's correction is commented out, L.Urban
// 29-05-01 V.Ivanchenko minor changes to provide ANSI -wall compilation
// 10-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 28-08-01 new function ComputeRestrictedMeandEdx() + 'cleanup' (mma)
// 17-09-01 migration of Materials to pure STL (mma)
// 26-09-01 completion of RetrievePhysicsTable (mma)
// 29-10-01 all static functions no more inlined (mma)
// 07-11-01 correction(Tmax+xsection computation) L.Urban
// 08-11-01 particleMass becomes a local variable (mma)
// 10-05-02 V.Ivanchenko update to new design
// 04-12-02 V.Ivanchenko the low energy limit for Kokoulin model to 10 GeV
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 26-12-02 Secondary production moved to derived classes (V.Ivanchenko)
// 13-02-03 SubCutoff regime is assigned to a region (V.Ivanchenko)
// 23-05-03 Define default integral + BohrFluctuations (V.Ivanchenko)
// 03-06-03 Add SetIntegral method to choose fluctuation model (V.Ivanchenko)
// 03-06-03 Fix initialisation problem for STD ionisation (V.Ivanchenko)
// 04-08-03 Set integral=false to be default (V.Ivanchenko)
// 08-08-03 STD substitute standard  (V.Ivanchenko)
// 12-11-03 G4EnergyLossSTD -> G4EnergyLossProcess (V.Ivanchenko)
// 10-02-04 Calculation of radiative corrections using R.Kokoulin model (V.Ivanchenko)
// 27-05-04 Set integral to be a default regime (V.Ivanchenko)
// 17-08-04 Utilise mu+ tables for mu- (V.Ivanchenko)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 12-08-05 SetStepLimits(0.2, 0.1*mm) (mma)
// 02-09-05 SetStepLimits(0.2, 1*mm) (V.Ivantchenko)
// 12-08-05 SetStepLimits(0.2, 0.1*mm) + integral off (V.Ivantchenko)
// 10-01-06 SetStepLimits -> SetStepFunction (V.Ivantchenko)
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "NumiMuIonisation.hh"
#include "G4LossTableManager.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VEmModel.hh"
#include "G4VEmFluctuationModel.hh"
#include "G4DataVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4VParticleChange.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4ProcessManager.hh"
#include "G4UnitsTable.hh"
#include "G4GenericIon.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4SafetyHelper.hh"
#include "G4TransportationManager.hh"

#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4BraggModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4MuBetheBlochModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4BohrFluctuations.hh"
#include "G4UnitsTable.hh"

#include "G4LossTableBuilder.hh"

#include "NumiAnalysis.hh"

#include "NumiMuIonisation.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

NumiMuIonisation::NumiMuIonisation(const G4String& name, G4ProcessType aType)
   : //G4VEnergyLossProcess(name),
   G4VEnergyLoss(name, aType),
   fName(name),
    theParticle(0),
    theBaseParticle(0),
    isInitialised(false),
    isIonisation(true),
    verboseLevel(1)
{
   SetStepFunction(0.2, 1*mm);
   //SetIntegral(true);

   verboseLevel = 1;

   
   
   lowestKinEnergy  = 1.*eV;
   minKinEnergy     = 0.1*keV;
   maxKinEnergy     = 100.0*TeV;

   modelManager = new G4EmModelManager();
  
   // default dRoverRange and finalRange
   SetStepFunction(0.2, 1.0*mm);
   
   //fParticleChange = new G4ParticleChangeForLoss();
   fParticleChange = new G4ParticleChangeForLoss();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

NumiMuIonisation::~NumiMuIonisation()
{
   if(theDEDXTable)         theDEDXTable->clearAndDestroy();
   if(theRangeTableForLoss) theRangeTableForLoss->clearAndDestroy();
   if(theInverseRangeTable) theInverseRangeTable->clearAndDestroy();
   delete modelManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void NumiMuIonisation::PreparePhysicsTable(const G4ParticleDefinition& part)
{
   NumiMuIonisation::InitialiseEnergyLossProcess(&part, theBaseParticle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void NumiMuIonisation::InitialiseEnergyLossProcess(const G4ParticleDefinition* part,
                                                 const G4ParticleDefinition* bpart)
{
   if(1 < verboseLevel) 
      G4cout << "NumiMuIonisation::InitialiseEnergyLossProcess ...." << G4endl;
   
  if(!isInitialised)
  {
     theParticle = part;
     //theBaseParticle = bpart;
     theBaseParticle = bpart;
     
     mass = theParticle->GetPDGMass();
     theSecondaryParticle = G4Electron::Electron();
     
     flucModel = new G4UniversalFluctuation();
     
     G4VEmModel* em = new G4BraggModel();
     em->SetLowEnergyLimit(0.1*keV);
     em->SetHighEnergyLimit(0.2*MeV);
     NumiMuIonisation::AddEmModel(1, em, flucModel);
     
     G4VEmModel* em1 = new G4BetheBlochModel();
     em1->SetLowEnergyLimit(0.2*MeV);
     em1->SetHighEnergyLimit(1.0*GeV);
     NumiMuIonisation::AddEmModel(2, em1, flucModel);

     G4VEmModel* em2 = new G4MuBetheBlochModel();
     em2->SetLowEnergyLimit(1.0*GeV);
     em2->SetHighEnergyLimit(100.0*TeV);
     NumiMuIonisation::AddEmModel(3, em2, flucModel);
     
    ratio = electron_mass_c2/mass;
    isInitialised = true;
  }

  G4double initialCharge = theParticle->GetPDGCharge();
  G4double initialMass   = theParticle->GetPDGMass();
  //chargeSquare = initialCharge*initialCharge/(eplus*eplus);
  chargeSqRatio = 1.0;
  massRatio = 1.0;
  reduceFactor = 1.0;
  
  if (theBaseParticle)
  {
     massRatio = (theBaseParticle->GetPDGMass())/initialMass;
     G4double q = initialCharge/theBaseParticle->GetPDGCharge();
     chargeSqRatio = q*q;
     if(chargeSqRatio > 0.0) reduceFactor = 1.0/(chargeSqRatio*massRatio);
  }


  theCuts = modelManager->Initialise(theParticle, G4Electron::Electron(), 0.1, 0);

   G4cout << "  Building Physics Tables ... " << G4endl;

   G4LossTableBuilder *tableBuilder = new G4LossTableBuilder();
   
   if(!theDEDXTable) theDEDXTable = G4PhysicsTableHelper::PreparePhysicsTable(theDEDXTable);
   theDEDXTable = NumiMuIonisation::BuildDEDXTable();
   
   if(!theRangeTableForLoss) theRangeTableForLoss  = G4PhysicsTableHelper::PreparePhysicsTable(theRangeTableForLoss);
   if(!theInverseRangeTable) theInverseRangeTable  = G4PhysicsTableHelper::PreparePhysicsTable(theInverseRangeTable);

   tableBuilder -> BuildRangeTable(theDEDXTable, theRangeTableForLoss, isIonisation);
   tableBuilder -> BuildInverseRangeTable(theRangeTableForLoss, theInverseRangeTable, isIonisation);


  if(1 < verboseLevel)
     G4cout << "NumiMuIonisation::InitialiseEnergyLossProcess - done" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4PhysicsTable* NumiMuIonisation::BuildDEDXTable()
{
  if(1 < verboseLevel) {
    G4cout << "NumiMuIonisation::BuildDEDXTable() for " << NumiMuIonisation::GetProcessName()
           << " and particle " << theParticle->GetParticleName()
           << G4endl;
  }
  G4PhysicsTable* table = 0;
  G4double emin = minKinEnergy;
  G4double emax = maxKinEnergy;
  G4int bin = 120;

  table = theDEDXTable;
  

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  if(1 < verboseLevel) {
    G4cout << numOfCouples << " materials"
           << " minKinEnergy= " << minKinEnergy
           << " maxKinEnergy= " << maxKinEnergy
           << " table= " << table
           << G4endl;
  }


  for(size_t i=0; i<numOfCouples; i++)
  {
     
    if(1 < verboseLevel) 
      G4cout << "NumiMuIonisation::BuildDEDXVector flag=  " 
	     << table->GetFlag(i) << G4endl;

    if (table->GetFlag(i))
    {

       // create physics vector and fill it
       const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
       G4PhysicsVector* aVector = new G4PhysicsLogVector(emin, emax, bin);
       modelManager->FillDEDXVector(aVector, couple, fRestricted);
       
      // Insert vector for this material into the table
      G4PhysicsTableHelper::SetPhysicsVector(table, i, aVector);
    }
  }

  if(1 < verboseLevel) {
    G4cout << "NumiMuIonisation::BuildDEDXTable(): table is built for "
           << theParticle->GetParticleName()
           << " and process " << NumiMuIonisation::GetProcessName()
           << G4endl;
    //    if(2 < verboseLevel) G4cout << (*table) << G4endl;
  }

  return table;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void NumiMuIonisation::PrintInfo()
{
  G4cout << "      Bether-Bloch model for E > 0.2 MeV, "
         << "parametrisation of Bragg peak below, "
         << G4endl;
  G4cout << "      radiative corrections for E > 1 GeV" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4VParticleChange* NumiMuIonisation::PostStepDoIt(const G4Track& track,
                                                  const G4Step&)
{
  fParticleChange ->InitializeForPostStep(track);

  return fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* NumiMuIonisation::AlongStepDoIt(const G4Track& track,
                                                   const G4Step& step)
{

   if(1 < verboseLevel) 
      G4cout << "NumiMuIonisation::AlongStepDoIt - called" << G4endl;

   G4double linLossLimit = 0.05;
   
  fParticleChange -> InitializeForAlongStep(track);

  NumiMuIonisation::InitialiseStep(track);

  G4String preStepName = step.GetPreStepPoint()->GetPhysicalVolume()->GetName();
  if(preStepName != "MuCell")
     return fParticleChange;

  G4cout << "****Pre step point = " << step.GetPreStepPoint()->GetPhysicalVolume()->GetName() << G4endl;

  G4double range = NumiMuIonisation::GetScaledRangeForScaledEnergy(preStepScaledEnergy)*reduceFactor;

  G4cout << "**********************range = " << range << G4endl;

  // Get the actual (true) Step length
  G4double length = step.GetStepLength();
  if(length <= DBL_MIN) return fParticleChange;

  G4cout << "**********************length = " << length << G4endl;
  
  // stopping
  if (length >= range)
  {
     G4cout << "NumiMuIonisation::AlongStepDoIt - PROBLEM: Stopping particle: range = " << fRange << G4endl;
     return fParticleChange;
  }
  

  NumiAnalysis* analysis = NumiAnalysis::getInstance();

  
  if(3 < verboseLevel)
  {
     const G4ParticleDefinition* d = track.GetDefinition();
     G4cout << "AlongStepDoIt for "
            << NumiMuIonisation::GetProcessName() << " and particle "
            << d->GetParticleName()
            << "  eScaled(MeV)= " << preStepScaledEnergy/MeV
            << "  range(mm)= " << fRange/mm
            << "  s(mm)= " << length/mm
            << "  q^2= " << chargeSqRatio
            << " md= " << d->GetPDGMass()
            << "  status= " << track.GetTrackStatus()
            << G4endl;
  }
  
      


 
  G4double eloss  = 0.0;
  
  
  // Short step
  eloss = NumiMuIonisation::GetDEDXForScaledEnergy(preStepScaledEnergy)*length;
  
  //  if( length <= linLossLimit * fRange ) {
  //  eloss = GetDEDXForScaledEnergy(preStepScaledEnergy)*length;
  
  // Long step
  //} else {
  if(eloss > preStepKinEnergy*linLossLimit)
  {
     
     G4cout << "  ****LONG STEP****" << G4endl;
     
     // G4double r = GetScaledRangeForScaledEnergy(preStepScaledEnergy);
     G4double r = range/reduceFactor;
     G4double x = r - length/reduceFactor;
     eloss = (NumiMuIonisation::ScaledKinEnergyForLoss(r) - NumiMuIonisation::ScaledKinEnergyForLoss(x))/massRatio;
     
     
     if(3 < verboseLevel)
     {
        G4cout << "Long STEP: rPre(mm)= " << r/mm
               << " rPost(mm)= " << x/mm
               << " ePre(MeV)= " << preStepScaledEnergy/MeV
               << " eloss(MeV)= " << eloss/MeV
               << " eloss0(MeV)= "
               << NumiMuIonisation::GetDEDXForScaledEnergy(preStepScaledEnergy)*length/MeV
               << G4endl;
     }
  }
  
  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();
  G4VEmModel* currentModel = NumiMuIonisation::SelectModel(preStepScaledEnergy);

  if(!currentModel)
  {
     G4cout << "NumiMuIonisation::AlongStepDoIt - PROBLEM Can't get current model" << G4endl;
  }
  
  if(3 < verboseLevel )
  {
     G4cout << "Before fluct: eloss(MeV)= " << eloss/MeV
            << " e-eloss= " << preStepKinEnergy-eloss
            << " step(mm)= " << length/mm
            << " range(mm)= " << range/mm
            << G4endl;
  }
  
  
  G4double cut  = (*theCuts)[currentMaterialIndex];
  
  
  
  // Corrections, which cannot be tabulated
  
  //
  //this function is empty in G4VEnergyLossProcess
  //
  //CorrectionsAlongStep(currentCouple, dynParticle, eloss, length);
  
  G4cout << "mean eloss  = " << eloss << G4endl;
    
  
  if((eloss + lowestKinEnergy) < preStepKinEnergy)
  {
     
     G4double tmax = std::min(currentModel->MaxSecondaryKinEnergy(dynParticle),cut);
     
     //
     // Sample fluctuations
     //
     G4cout << "fluctuations:  " ;
     for(int i = 1; i < 100; ++i)
     {
        G4double MeanELoss = eloss;
        G4double eflucloss = -999.0;
        
        G4VEmFluctuationModel* flucMod = currentModel->GetModelOfFluctuations();
        if(!flucMod)
        {
           G4cout << "NumiMuIonisation::AlongStepDoIt - PROBLEM Can't get current fluc model" << G4endl;
        }
        eflucloss = flucMod ->
           SampleFluctuations(currentMaterial,dynParticle,tmax,length,MeanELoss);
           
        if(3 < verboseLevel)
        {
           G4cout << "After fluct: eloss(MeV)= " << eflucloss/MeV
                  << " fluc= " << (eflucloss-MeanELoss)/MeV
              //<< " currentChargeSquare= " << chargeSquare
                  << " massRatio= " << massRatio
                  << " tmax= " << tmax
                  << G4endl;
        }
        
        
        G4cout << eflucloss << ", ";
     }
     G4cout << G4endl;
     //
     //
     //
  }
  else
     G4cout << "NumiMuIonisation::AlongStepDoIt - PROBLEM in trying to sample fluctuations" << G4endl;
  
  
  
  
  if(3 < verboseLevel)
  {
     G4cout << "Final value eloss(MeV)= " << eloss/MeV
            << " preStepKinEnergy= " << preStepKinEnergy
            << "  status= " << track.GetTrackStatus()
            << G4endl;
  }
  
  
  
  if(1 < verboseLevel) 
     G4cout << "NumiMuIonisation::AlongStepDoIt - done" << G4endl;
  
  
  
  return fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void NumiMuIonisation::SetStepFunction(G4double v1, G4double v2)
{
  dRoverRange = v1;
  finalRange = v2;
  if (dRoverRange > 0.999) dRoverRange = 1.0;
  currentCouple = 0;
  mfpKinEnergy  = DBL_MAX;
}

