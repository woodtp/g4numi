#include "NumiDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
#include "NumiDataInput.hh"

void NumiDetectorConstruction::DefineMaterials()
{  
  //------------------------------------------------------ materials

  G4double A;  // atomic mass
  G4double Z;  // atomic number
  G4String symbol,name;
  G4double density;

  // Air and Vacuum
  A = 1.01*g/mole;
  G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , Z= 1, A);

  A = 12.01*g/mole;
  G4Element* elC  = new G4Element(name="Carbon"  ,symbol="C" , Z = 6, A);

  A = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", Z=7.,A);

  A = 16.00*g/mole;
  G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , Z= 8, A);

  A = 22.99*g/mole; 
  G4Element* elNa  = new G4Element(name="Natrium"  ,symbol="Na" , Z=11 , A);

  A = 26.98*g/mole; 
  G4Element* elAl  = new G4Element(name="Aluminum"  ,symbol="Al" , Z=13, A);

  A = 28.09*g/mole;
  G4Element* elSi  = new G4Element(name="Silicon", symbol="Si", Z=14, A);

  A = 30.974*g/mole; 
  G4Element* elP  = new G4Element(name="Phosphorus"  ,symbol="P" , Z=15 , A);

  A = 32.065*g/mole; 
  G4Element* elS  = new G4Element(name="Sulfur"  ,symbol="S" , Z=16 , A);

  A = 39.948*g/mole;
  G4Element* elAr = new G4Element(name="Argon" , symbol="Ar", Z=18, A);

  A = 39.1*g/mole; 
  G4Element* elK  = new G4Element(name="Potassium"  ,symbol="K" , Z=19 , A);

  A = 47.867*g/mole; 
  G4Element* elTi  = new G4Element(name="Titanium"  ,symbol="Ti" , Z=22 , A);

  A = 51.9961*g/mole; 
  G4Element* elCr  = new G4Element(name="Chromium"  ,symbol="Cr" , Z=24 , A);

  A = 54.938*g/mole; 
  G4Element* elMn  = new G4Element(name="Manganese"  ,symbol="Mn" , Z=25 , A);

  A = 55.85*g/mole;
  G4Element* elFe = new G4Element(name="Iron"    ,symbol="Fe", Z=26, A);

  A = 58.6934*g/mole;
  G4Element* elNi = new G4Element(name="Nickel"    ,symbol="Ni", Z=28, A);

  A = 63.546*g/mole;
  G4Element* elCu = new G4Element(name="Copper"    ,symbol="Cu", Z=29, A);

  A = 69.72*g/mole; 
  G4Element* elCa  = new G4Element(name="Calcium"  ,symbol="Ca" , Z=31 , A);

  A = 200.59*g/mole; 
  G4Element* elHg  = new G4Element(name="Mercury"  ,symbol="Hg" , Z=80, A);

  density = 1.29*mg/cm3; 
  Air = new G4Material("Air", density, 2); //number of components =2
  Air->AddElement(elN, 70*perCent); //mass fraction =70%
  Air->AddElement(elO, 30*perCent); // 30%

  density = 7.75*g/cm3;  
  CT852 = new G4Material("CT852", density, 10); 
  CT852->AddElement(elC,  0.1*perCent); 
  CT852->AddElement(elSi, 0.8*perCent); 
  CT852->AddElement(elMn, 0.8*perCent); 
  CT852->AddElement(elCr, 13*perCent); 
  CT852->AddElement(elS,  0.025*perCent); 
  CT852->AddElement(elP,  0.03*perCent); 
  CT852->AddElement(elTi, 0.2*perCent); 
  CT852->AddElement(elCu, 0.3*perCent); 
  CT852->AddElement(elNi, 0.6*perCent); 
  CT852->AddElement(elFe, 84.145*perCent); 

  density = 1.*g/cm3; 
  G4int natoms;
  Water = new G4Material("Water", density, 2); //number of components =2
  Water->AddElement(elH, natoms=2); 
  Water->AddElement(elO, natoms=1); 

  density=2.376e-15*g/cm3;
  G4double temperature=300.*kelvin;
  G4double pressure=2.0e-7*bar;
  Vacuum = new G4Material("Vacuum", density, 1, kStateGas,temperature,pressure);
  Vacuum-> AddMaterial(Air, 1.);

  density=1.29/760.*mg/cm3;
  temperature=300.*kelvin;
  pressure=atmosphere/760.;
  DecayPipeVacuum = new G4Material("DecayPipeVacuum", density, 1, kStateGas,temperature,pressure);
  DecayPipeVacuum-> AddMaterial(Air, 1.);

  //other materials  
  Be = new G4Material("Berillium", Z=4.,A=9.01*g/mole, density=1.848*g/cm3);
  C =  new G4Material("Carbon", Z=6., A=12.01*g/mole, density= 1.83*g/cm3);
  Al = new G4Material("Aluminum", Z= 13., A= 26.98*g/mole, density= 2.7*g/cm3);
  Ar = new G4Material("Argon", Z= 18, A=39.948*g/mole,1.784*kg/m3,kStateGas,300*kelvin,atmosphere);
  Pb = new G4Material("Lead", Z= 82., A= 207.19*g/mole, density= 11.35*g/cm3);
  Fe = new G4Material("Iron", Z= 26., A=55.85*g/mole, density= 7.86999*g/cm3);
  Target =  new G4Material("Target", Z=NumiData->TargetZ, A=NumiData->TargetA, density= NumiData->TargetDensity);
  
  density = 2.03*g/cm3;
  G4double fractionmass;
  Concrete = new G4Material("Concrete", density, 10);
  Concrete->AddElement(elH , fractionmass= 0.01);
  Concrete->AddElement(elO , fractionmass= 0.529);
  Concrete->AddElement(elNa , fractionmass= 0.016);
  Concrete->AddElement(elHg , fractionmass= 0.002);
  Concrete->AddElement(elAl , fractionmass= 0.034);
  Concrete->AddElement(elSi , fractionmass= 0.337);
  Concrete->AddElement(elK , fractionmass= 0.013);
  Concrete->AddElement(elCa , fractionmass= 0.044);
  Concrete->AddElement(elFe , fractionmass= 0.014);
  Concrete->AddElement(elC , fractionmass= 0.001);
}

G4Material* NumiDetectorConstruction::GetMaterial(G4int matcode)
{

  if (matcode==15) return Air;
  if (matcode==6) return C; //?
  if (matcode==9 || matcode==20) return Al;
  if (matcode==5) return Be; 
  if (matcode==10) return Fe;
  if (matcode==31) return CT852;
  if (matcode==17) return Concrete;
  if (matcode==18) return Target;
  if (matcode==16) return Vacuum;
  if (matcode==25) return Water;
  G4cout << "Wrong material code " << matcode << G4endl;
  return Vacuum;
}
G4VisAttributes* NumiDetectorConstruction::GetMaterialVisAttrib(G4int matCode)
{ 
  G4VisAttributes* visAttrib=new G4VisAttributes(G4Color(1.,0.,0.));
  if (matCode==15) visAttrib=new G4VisAttributes(false); //Air
  if (matCode==6) visAttrib=new G4VisAttributes(G4Color(0.7,0.7,0.7)); //C
  if (matCode==9 || matCode==20) visAttrib=new G4VisAttributes(G4Color(0.2,0.8,1.));//Al
  if (matCode==5) visAttrib=new G4VisAttributes(G4Color(0.1,0.2,0.95));//Be 
  if (matCode==10) visAttrib=new G4VisAttributes(G4Color(0.5,0.5,0.5));//Fe
  if (matCode==31) visAttrib=new G4VisAttributes(G4Color(1.,1.,1.));//CT852
  if (matCode==17) visAttrib=new G4VisAttributes(G4Color(0.85,0.85,0.85));//Concrete
  if (matCode==18) visAttrib=new G4VisAttributes(G4Color(0.6,0.6,0.7));//Target
  if (matCode==16) visAttrib=new G4VisAttributes(false);//Vacuum
  if (matCode==25) visAttrib=new G4VisAttributes(G4Color(0.,0.,1.));//Water
  
  return visAttrib;
}
G4VisAttributes* NumiDetectorConstruction::GetMaterialVisAttrib(G4String matName){
  G4VisAttributes* visAttrib=new G4VisAttributes(G4Color(1.,0.,0.));
  if (matName=="Air") return GetMaterialVisAttrib(15);
  if (matName=="Carbon") return GetMaterialVisAttrib(6);
  if (matName=="Aluminum") return GetMaterialVisAttrib(9);
  if (matName=="Berillium") return GetMaterialVisAttrib(5);
  if (matName=="Iron") return GetMaterialVisAttrib(10);
  if (matName=="CT852") return GetMaterialVisAttrib(31);
  if (matName=="Concrete") return GetMaterialVisAttrib(17);
  if (matName=="Target") return GetMaterialVisAttrib(18);
  if (matName=="Vacuum") return GetMaterialVisAttrib(16);
  if (matName=="DecayPipeVacuum") return GetMaterialVisAttrib(16);
  if (matName=="Argon") return GetMaterialVisAttrib(16);
  if (matName=="Water") return GetMaterialVisAttrib(25);
  return visAttrib;
}

void NumiDetectorConstruction::DestroyMaterials()
{
  // Destroy all allocated elements and materials
  size_t i;
  G4MaterialTable* matTable = (G4MaterialTable*)G4Material::GetMaterialTable();
  for(i=0;i<matTable->size();i++)
  { delete (*(matTable))[i]; }
  matTable->clear();
  G4ElementTable* elemTable = (G4ElementTable*)G4Element::GetElementTable();
  for(i=0;i<elemTable->size();i++)
  { delete (*(elemTable))[i]; }
  elemTable->clear();
}
