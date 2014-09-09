
#include "G4MinervaProtonIXS.hh"
#include "globals.hh"

G4MinervaProtonIXS::G4MinervaProtonIXS(){
  scaleVal = 1.;
}

G4double G4MinervaProtonIXS::
GetCrossSection(const G4DynamicParticle* aPart, 
                const G4Element* anEle, G4double /*aTemperature*/)
{
  
  G4int nIso = anEle->GetNumberOfIsotopes();
  G4double KE = aPart->GetKineticEnergy(); 
  G4double cross_section = 0;
  
  if (nIso) {
    G4double psig;
    G4IsotopeVector* isoVector = anEle->GetIsotopeVector();
    G4double* abundVector = anEle->GetRelativeAbundanceVector();
    G4double ZZ;
    G4double AA;
 
    for (G4int i = 0; i < nIso; i++) {
      ZZ = G4double( (*isoVector)[i]->GetZ() );
      AA = G4double( (*isoVector)[i]->GetN() );
      psig = GetCrossSection(KE, AA, ZZ);
      cross_section += psig*abundVector[i];
    }
 
  } else {
    cross_section = GetCrossSection(KE, anEle->GetN(), anEle->GetZ());
  }
  return cross_section;
}

   
G4double G4MinervaProtonIXS::
GetCrossSection(G4double kineticEnergy, G4double atomicNumber, G4double nOfProtons)
{   

  if (kineticEnergy > 19.9*GeV ) 
  { // constant cross section above ~20GeV.
    return  GetCrossSection(19.8*GeV,atomicNumber,nOfProtons);
  } 
  G4double nOfNeutrons = atomicNumber-nOfProtons;
  kineticEnergy /=GeV;
  G4double a = atomicNumber;
  const G4double nuleonRadius=1.36E-15;
  const G4double pi=3.14159265;
  G4double fac=pi*nuleonRadius*nuleonRadius;
  G4double b0=2.247-0.915*(1-std::pow(a,-0.3333));
  G4double fac1=b0*(1-std::pow(a,-0.3333));
  G4double fac2=1.;
  if(nOfNeutrons>1.5) fac2=std::log((nOfNeutrons));
  G4double crossSection = 1E31*fac*fac2*(1+std::pow(a,0.3333)-fac1);

  // high energy correction

  crossSection = (1-0.15*std::exp(-kineticEnergy))*crossSection/(1.00-0.0007*a);
  // first try on low energies: rise

  G4double ff1= 0.70-0.002*a;  // slope of the drop at medium energies.
  G4double ff2= 1.00+1/a;  // start of the slope.
  G4double ff3= 0.8+18/a-0.002*a; // step height
  fac = 1.0;
  if (kineticEnergy > DBL_MIN) 
     fac= 1.0 - (1.0/(1+std::exp(-8*ff1*(std::log10(kineticEnergy)+1.37*ff2))));
  crossSection = crossSection*(1+ff3*fac);

  // low energy return to zero

  ff1=1.-1/a-0.001*a; // slope of the rise
  ff2=1.17-2.7/a-0.0014*a; // start of the rise
  fac = 0.0;
  if (kineticEnergy > DBL_MIN) {
    fac=-8.*ff1*(std::log10(kineticEnergy)+2.0*ff2);
    fac=1/(1+std::exp(fac));
  }
  crossSection = crossSection*fac;
  //  G4cout<<"scaleVal prt IXS "<<scaleVal<<" "<<scaleVal*crossSection<<G4endl;
  return scaleVal*crossSection*millibarn;
}

