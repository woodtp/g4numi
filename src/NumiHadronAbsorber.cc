// The new and improved Hadron Absorber my D. Jason Koskinen
// Last Update: March 5, 2006
#include "NumiDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"
#include "G4Cons.hh"
#include "G4Torus.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4PVReplica.hh"
#include "G4AssemblyVolume.hh"
#include "NumiDataInput.hh"
#include <math.h>

void NumiDetectorConstruction::ConstructHadronAbsorber()
{
  G4ThreeVector tunnelPos = G4ThreeVector(0,0,NumiData->TunnelLength/2.+NumiData->TunnelZ0);
  G4double in = .0254;
  G4float CShld_length = 72*in*m;
  
  //blue block dimensions, all lengths in Geant4 are half-lengths
  G4float bb_length=1.3208*m;
  G4float bb_width=1.3208*m;
  G4float bb_height=.6604*m;

  // Core Aluminum and Steel block dimensions
  G4float Al_width = 51*in*m;
  G4float Al_height = 51*in*m;
  G4float Al_length = 12*in*m;
  G4float Stl_width = 1.2954*m;
  G4float Stl_height = 1.2954*m;
  G4float Stl_length = .23114*m;

  // Concrete Shielding block dimensions
  G4float A_width = 1.5*12*in*m;
  G4float A_height = 7.5*12*in*m;
  G4float A_length = 3*12*in*m;
  G4float B_width = 3*12*in*m;
  G4float B_height = 7.5*12*in*m;
  G4float B_length = 3*12*in*m;
  G4float C_width = 3*12*in*m;
  G4float C_height = 3*12*in*m;
  G4float C_length = 6*12*in*m;
  G4float D_width = 1.5*12*in*m;
  G4float D_height = 3*12*in*m;
  G4float D_length = 6*12*in*m;
  G4float J_width = 12*12*in*m;
  G4float J_height = 1.5*12*in*m;
  G4float J_length = 3*12*in*m;

  // Only want to retrieve some information once, i.e. z position of the end of the decay pipe.
  G4float DP_end = NumiData->DecayPipeZ0+NumiData->DecayPipeLength;

  // rotation crap. Establishes the location of the concrete shielding and Absorber box in relation to a 3.3 degree rotation
  G4float D_hadr = 50.25*in*m+NumiData->HadrBox_length/2.0;
  G4float H_hadr = sqrt(pow(11.13*in*m,2)+pow(D_hadr,2));
  G4float a_hadr = atan(11.1325*in*m/D_hadr);
  G4float b_hadr = 3.3*pi/180-a_hadr;
  G4float z_hadr = H_hadr*cos(b_hadr)+DP_end;
  G4float y_hadr = H_hadr*sin(b_hadr);

  G4float D_shld = (50.25-36)*in*m;
  G4float H_shld = sqrt(pow(11.13*in*m,2)+pow(D_shld,2));
  G4float a_shld = atan(11.1325*in*m/D_shld);
  G4float b_shld = 3.3*pi/180-a_shld;
  G4float z_shld = H_shld*cos(b_shld)+DP_end;
  G4float y_shld = H_shld*sin(b_shld);
  
  G4Material* blu_material=GetMaterial(11);
  G4Material* al=GetMaterial(9);
  G4Material* concrete=GetMaterial(17);
  G4Material* var_Al=GetMaterial(21);
  G4Material* var_Stl=GetMaterial(22);
  char no[3];
  G4String volName;

  // This constructs a Hadron Absorber volume to put everything in. 
  //------------------------------------------------------------------------

  G4Box*  expHadrBox = new G4Box("expHadrBox",NumiData->HadrBox_width/2.,NumiData->HadrBox_height/2.0,NumiData->HadrBox_length/2.0);
  HadrBox_log = new G4LogicalVolume(expHadrBox,Air,"LVHadrBox",0,0,0);
  G4RotationMatrix *zrot=new G4RotationMatrix();
  zrot->rotateX(3.3*M_PI/180);
  G4ThreeVector pos=G4ThreeVector(0,y_hadr,z_hadr)-tunnelPos;
  G4String nam="HadronAbsorber";
  //(rotation,G4ThreeVector,name,logical volume,mother,boolean operators,copy number)
  HadrBox=new G4PVPlacement(zrot,pos,nam,HadrBox_log,pvTUNE,false,0);

  // The following code was introduced with the idea to include 
  // some variable density bits of the stl and al core in order
  // to simulate the cooling pipes. 
  //------------------------------------------------------------------------

  G4Box* hole1_solid = new G4Box("hole1",4.5*in*m/2.0,42.75*in*m/2.0,12*in*m/2.0);
  G4Box* hole2_solid = new G4Box("hole2",4.5*in*m/2.0,42.75*in*m/2.0,9.1*in*m/2.0);
  G4LogicalVolume* hole1_log = new G4LogicalVolume(hole1_solid,var_Al,"hole1",0,0,0);
  G4LogicalVolume* hole2_log = new G4LogicalVolume(hole2_solid,var_Stl,"hole2",0,0,0);
  G4ThreeVector holePosL = G4ThreeVector(-22*in*m,-1.875*in*m,0);
  G4ThreeVector holePosR = G4ThreeVector(22*in*m,-1.875*in*m,0);
  G4ThreeVector blockPos;

  G4float Alx_pos = 0;
  G4float Aly_pos = 0;
  G4float Alz_pos = 0;
  G4float core_offset = 12.75*in*m; // This is the offset between the beginning edge of the blublocks and the the front edge of the core.

  // The Absorber core positioning and placement
  //-----------------------------------------------------------------------

  for(G4int ii=0;ii<8;ii++){
    Alz_pos=ii*Al_length+Al_length/2.0-NumiData->HadrBox_length/2.0+core_offset; 
    sprintf(no,"%d",ii+1);
    volName="Al_sBLK"; 
    volName.append(no);
    BLK_solid[ii] = new G4Box(volName,Al_width/2.0,Al_height/2.0,Al_length/2.0);
    volName="Al_lvBLK"; 
    volName.append(no);
    BLK_log[ii] = new G4LogicalVolume(BLK_solid[ii],al,volName,0,0,0);
    volName="Al_BLK"; 
    volName.append(no);
    blockPos=G4ThreeVector(Alx_pos,Aly_pos,Alz_pos);
    new G4PVPlacement(0,holePosL,hole1_log,"hole",BLK_log[ii],false,0);
    new G4PVPlacement(0,holePosR,hole1_log,"hole",BLK_log[ii],false,0);
    new G4PVPlacement(0,blockPos,volName,BLK_log[ii],HadrBox,false,0);
  }	

  G4float Stlx_pos = 0;
  G4float Stly_pos = 0;
  G4float Stlz_pos = 0;
  //Stl core block positioning
  for(G4int ii=0;ii<10;ii++){
    Stlz_pos=ii*Stl_length+Stl_length/2.0-NumiData->HadrBox_length/2.0+8*Al_length+core_offset; 
    sprintf(no,"%d",ii+1);
    volName="Stl_sBLK"; 
    volName.append(no);
    BLK_solid[ii] = new G4Box(volName,Stl_width/2.0,Stl_height/2.0,Stl_length/2.0);
    volName="Stl_lvBLK"; 
    volName.append(no);
    BLK_log[ii] = new G4LogicalVolume(BLK_solid[ii],al,volName,0,0,0);
    volName="Stl_BLK"; 
    volName.append(no);
    blockPos=G4ThreeVector(Stlx_pos,Stly_pos,Stlz_pos);
    new G4PVPlacement(0,holePosL,hole2_log,"hole",BLK_log[ii],false,0);
    new G4PVPlacement(0,holePosR,hole2_log,"hole",BLK_log[ii],false,0);
    new G4PVPlacement(0,blockPos,volName,BLK_log[ii],HadrBox,false,0);
  }	

  // Base, side and top levels surrounding the Absorber Core
  //-------------------------------------------------------------------

  G4float bbx_pos = 0;
  G4float bby_pos = 0;
  G4float bbz_pos = 0;

  // Blu block positioning
  // Bottom two layers
  for (G4int ii=0;ii<24;ii++){
    bbx_pos=(ii%3)*bb_width-bb_width;
    bby_pos=-ii/12*bb_height-Al_height/2.0-bb_height/2.0-2*bb_height-in*m;
    bbz_pos=(ii%12)/3*bb_length-NumiData->HadrBox_length/2.0+bb_length/2.0; 
    sprintf(no,"%d",ii+1);
    volName="blu_sBLK"; 
    volName.append(no);
    BLK_solid[ii] = new G4Box(volName,bb_width/2.0,bb_height/2.0,bb_length/2.0);
    volName="blu_lvBLK"; 
    volName.append(no);
    BLK_log[ii] = new G4LogicalVolume(BLK_solid[ii],blu_material,volName,0,0,0);
    volName="blu_BLK"; 
    volName.append(no);
     blockPos=G4ThreeVector(bbx_pos,bby_pos,bbz_pos);
    new G4PVPlacement(0,blockPos,volName,BLK_log[ii],HadrBox,false,0);
  }
  //layer immediatetly under the core
  for (G4int ii=0;ii<8;ii++){
    bbx_pos=0;
    bby_pos=-ii/4*bb_height-Al_height/2.0-bb_height/2.0-in*m;
    bbz_pos=ii%4*bb_length-NumiData->HadrBox_length/2.0+bb_length/2.0; 
    sprintf(no,"%d",ii+25);
    volName="blu_sBLK"; 
    volName.append(no);
    BLK_solid[ii] = new G4Box(volName,bb_width/2.0,bb_height/2.0,bb_length/2.0);
    volName="blu_lvBLK"; 
    volName.append(no);
    BLK_log[ii] = new G4LogicalVolume(BLK_solid[ii],blu_material,volName,0,0,0);
    volName="blu_BLK"; 
    volName.append(no);
     blockPos=G4ThreeVector(bbx_pos,bby_pos,bbz_pos);
    new G4PVPlacement(0,blockPos,volName,BLK_log[ii],HadrBox,false,0);
    if(ii==0){
      volName="stl_slab1_s";
      volName.append(no);
      G4Box *stl_slab1_s = new G4Box(volName,bb_width/2.0,in/2.0*m,bb_length/2.0*4);
      volName="stl_slab1_lv";
      volName.append(no);
      G4LogicalVolume *stl_slab1_lv = new G4LogicalVolume(stl_slab1_s,blu_material,volName,0,0,0);
      volName="stl_slab1";
      G4ThreeVector slabpos1 = G4ThreeVector(0,-Al_height/2.0-in*m/2.0,-NumiData->HadrBox_length/2.0+2*bb_length);
      new G4PVPlacement(0,slabpos1,volName,stl_slab1_lv,HadrBox,false,0);

      volName="stl_slab2_s";
      volName.append(no);
      G4Box *stl_slab2_s = new G4Box(volName,bb_width/2.0,in/2.0*m,bb_length/2.0*4);
      volName="stl_slab2_lv";
      G4LogicalVolume *stl_slab2_lv = new G4LogicalVolume(stl_slab2_s,blu_material,volName,0,0,0);
      volName="stl_slab2";
      volName.append(no);
      G4ThreeVector slabpos2 = G4ThreeVector(-bb_width,-bb_height*2.0-in*m/2.0-Al_height/2,-NumiData->HadrBox_length/2.0+2*bb_length);
      new G4PVPlacement(0,slabpos2,volName,stl_slab2_lv,HadrBox,false,0);

      volName="stl_slab3_s";
      G4Box *stl_slab3_s = new G4Box(volName,bb_width/2.0,in/2.0*m,bb_length/2.0*4);
      volName="stl_slab3_lv";
      volName.append(no);
      G4LogicalVolume *stl_slab3_lv = new G4LogicalVolume(stl_slab3_s,blu_material,volName,0,0,0);
      volName="stl_slab3";
      G4ThreeVector slabpos3 = G4ThreeVector(bb_width,-bb_height*2.0-in*m/2.0-Al_height/2,-NumiData->HadrBox_length/2.0+2*bb_length);
      new G4PVPlacement(0,slabpos3,volName,stl_slab3_lv,HadrBox,false,0);

      volName="stl_slab4_s";
      volName.append(no);
      G4Box *stl_slab4_s = new G4Box(volName,3.9624/2.0*m,in/2.0*m,5.2832/2.0*m);
      volName="stl_slab4_lv";
      volName.append(no);
      G4LogicalVolume *stl_slab4_lv = new G4LogicalVolume(stl_slab4_s,blu_material,volName,0,0,0);
      volName="stl_slab4";
      G4ThreeVector slabpos4 = G4ThreeVector(0,Al_height/2.0+in*m*1.5,-NumiData->HadrBox_length/2.0+2*bb_length);
      new G4PVPlacement(0,slabpos4,volName,stl_slab4_lv,HadrBox,false,0);
    }
  }

  // Creates the steel slabs wedged between the core
  // and the surrounding blu blocks to reduce the gap size.
  // Dixon estimates that the steel slabs are ~1 the length
  // of the blu blocks.
  volName="stl_slabL_s";
  volName.append(no);
  G4Box *stl_slabL_s = new G4Box(volName,1/2.0*in/2.0*m,Al_width/2.0,1*bb_length/2.0);
  volName="stl_slabL_lv";
  volName.append(no);
  G4LogicalVolume *stl_slabL_lv = new G4LogicalVolume(stl_slabL_s,blu_material,volName,0,0,0);
  volName="stl_slabL";
  G4ThreeVector slabpos5 = G4ThreeVector(-Al_width/2-1/2.0*in/2.0*m, 0, -NumiData->HadrBox_length/2.0+core_offset+8*Al_length+10*Stl_length-1*bb_length/2.0);
  new G4PVPlacement(0,slabpos5,volName,stl_slabL_lv,HadrBox,false,0);

  volName="stl_slabR_s";
  volName.append(no);
  G4Box *stl_slabR_s = new G4Box(volName,1/2.0*in/2.0*m,Al_width/2.0,1*bb_length/2.0);
  volName="stl_slabR_lv";
  volName.append(no);
  G4LogicalVolume *stl_slabR_lv = new G4LogicalVolume(stl_slabR_s,blu_material,volName,0,0,0);
  volName="stl_slabR";
  G4ThreeVector slabpos6 = G4ThreeVector(Al_width/2+1/2.0*in/2.0*m, 0, -NumiData->HadrBox_length/2.0+core_offset+8*Al_length+10*Stl_length-1*bb_length/2.0);
  new G4PVPlacement(0,slabpos6,volName,stl_slabR_lv,HadrBox,false,0);

  //Layer left(beam-view) and under of the core, and above the steel slab
  for (G4int ii=0;ii<8;ii++){
    bbx_pos=-bb_width;
    bby_pos=ii/4*bb_height-Al_height/2.0-bb_height*1.5;
    bbz_pos=ii%4*bb_length-NumiData->HadrBox_length/2.0+bb_length/2.0;
    sprintf(no,"%d",ii+33);
    volName="blu_sBLK"; 
    volName.append(no);
    BLK_solid[ii] = new G4Box(volName,bb_width/2.0,bb_height/2.0,bb_length/2.0);
    volName="blu_lvBLK"; 
    volName.append(no);
    BLK_log[ii] = new G4LogicalVolume(BLK_solid[ii],blu_material,volName,0,0,0);
    volName="blu_BLK"; 
    volName.append(no);
     blockPos=G4ThreeVector(bbx_pos,bby_pos,bbz_pos);
    new G4PVPlacement(0,blockPos,volName,BLK_log[ii],HadrBox,false,0);
  }

  for (G4int ii=0;ii<8;ii++){
    bbx_pos=bb_width;
    bby_pos=ii/4*bb_height-Al_height/2.0-bb_height*1.5;
    bbz_pos=ii%4*bb_length-NumiData->HadrBox_length/2.0+bb_length/2.0; 
    sprintf(no,"%d",ii+41);
    volName="blu_sBLK"; 
    volName.append(no);
    BLK_solid[ii] = new G4Box(volName,bb_width/2.0,bb_height/2.0,bb_length/2.0);
    volName="blu_lvBLK"; 
    volName.append(no);
    BLK_log[ii] = new G4LogicalVolume(BLK_solid[ii],blu_material,volName,0,0,0);
    volName="blu_BLK"; 
    volName.append(no);
     blockPos=G4ThreeVector(bbx_pos,bby_pos,bbz_pos);
    new G4PVPlacement(0,blockPos,volName,BLK_log[ii],HadrBox,false,0);
  }

  // Blocks to the right of the core (Beam-view looking downstream)
  for (G4int ii=0;ii<8;ii++){
    bbx_pos=bb_width/2.0+bb_height/2.0+(ii%2)*bb_height;
    bby_pos=in/2.0*m;
    bbz_pos=(ii/2)*bb_width-NumiData->HadrBox_length/2.0+bb_width/2.0; 
    sprintf(no,"%d",ii+49);
    volName="blu_sBLK"; 
    volName.append(no);
    BLK_solid[ii] = new G4Box(volName,bb_height/2.0,bb_length/2.0,bb_width/2.0);
    volName="blu_lvBLK"; 
    volName.append(no);
    BLK_log[ii] = new G4LogicalVolume(BLK_solid[ii],blu_material,volName,0,0,0);
    volName="blu_BLK"; 
    volName.append(no);
     blockPos=G4ThreeVector(bbx_pos,bby_pos,bbz_pos);
    new G4PVPlacement(0,blockPos,volName,BLK_log[ii],HadrBox,false,0);
  }

  // Blocks to the left of the core
  for (G4int ii=0;ii<8;ii++){
    bbx_pos=-bb_width/2.0-bb_height/2.0-(ii%2)*bb_height;
    bby_pos=in/2.0*m;
    bbz_pos=(ii/2)*bb_width-NumiData->HadrBox_length/2.0+bb_width/2.0; 
    sprintf(no,"%d",ii+57);
    volName="blu_sBLK"; 
    volName.append(no);
    BLK_solid[ii] = new G4Box(volName,bb_height/2.0,bb_length/2.0,bb_width/2.0);
    volName="blu_lvBLK"; 
    volName.append(no);
    BLK_log[ii] = new G4LogicalVolume(BLK_solid[ii],blu_material,volName,0,0,0);
    volName="blu_BLK"; 
    volName.append(no);
     blockPos=G4ThreeVector(bbx_pos,bby_pos,bbz_pos);
    new G4PVPlacement(0,blockPos,volName,BLK_log[ii],HadrBox,false,0);
  }

  // Blocks above the Core
  for (G4int ii=0;ii<24;ii++){
    bbx_pos=ii%3*(-bb_width)+bb_width;
    bby_pos=ii/12*bb_height+(in*m*2+Al_height/2.0+bb_height/2.0);
    bbz_pos=(ii%12)/3*(bb_length)-NumiData->HadrBox_length/2.0+bb_length/2.0; 
    sprintf(no,"%d",ii+65);
    volName="blu_sBLK"; 
    volName.append(no);
    BLK_solid[ii] = new G4Box(volName,bb_width/2.0,bb_height/2.0,bb_length/2.0);
    volName="blu_lvBLK"; 
    volName.append(no);
    BLK_log[ii] = new G4LogicalVolume(BLK_solid[ii],blu_material,volName,0,0,0);
    volName="blu_BLK"; 
    volName.append(no);
     blockPos=G4ThreeVector(bbx_pos,bby_pos,bbz_pos);
    new G4PVPlacement(0,blockPos,volName,BLK_log[ii],HadrBox,false,0);
  }

  // Concrete blocks surrounding core
  //------------------------------------------------------------------------

  G4float conc_width=.9144*m;
  G4float conc_length=7*.9144*m;
  G4float conc_height=2*2.1336*m;
  BLK_solid[90]= new G4Box("conc_sBLK",conc_width/2.0,conc_height/2.0,conc_length/2.0);
  BLK_log[90]=new G4LogicalVolume(BLK_solid[90],concrete,"conc_lvBLK",0,0,0);
   blockPos=G4ThreeVector(-81*in*m-conc_width/2.0,-in*m-Al_height/2.0-4*bb_height+conc_height/2.0,-NumiData->HadrBox_length/2.0+conc_length/2.0);
  new G4PVPlacement(0,blockPos,"conc_BLK",BLK_log[90],HadrBox,false,0);

  BLK_solid[91]= new G4Box("conc_sBLK",conc_width/2.0,conc_height/2.0,conc_length/2.0-.9144/2.0*m);
  BLK_log[91]=new G4LogicalVolume(BLK_solid[91],concrete,"conc_lvBLK",0,0,0);
  blockPos=G4ThreeVector(81*in*m+conc_width/2.0,-in*m-Al_height/2.0-4*bb_height+conc_height/2.0,-NumiData->HadrBox_length/2.0+(conc_length-.9144*m)/2.0);
  new G4PVPlacement(0,blockPos,"conc_BLK",BLK_log[91],HadrBox,false,0);

  BLK_solid[92]= new G4Box("conc_sBLK",conc_length/2.0-2*.9144/2.0*m,conc_height/2.0,conc_width/2.0);
  BLK_log[92]=new G4LogicalVolume(BLK_solid[92],concrete,"conc_lvBLK",0,0,0);
  blockPos=G4ThreeVector(81*in*m+conc_width/2.0-(conc_length-2*.9144*m)/2.0,-in*m-Al_height/2.0-4*bb_height+conc_height/2.0,-NumiData->HadrBox_length/2.0+conc_length+10*in*m-conc_width/2.0);
  new G4PVPlacement(0,blockPos,"conc_BLK",BLK_log[92],HadrBox,false,0);

  BLK_solid[93]= new G4Box("topstl1_sBLK",3*bb_width/2.0,18*in*m/2.0,4*bb_width/2.0);
  BLK_log[93]=new G4LogicalVolume(BLK_solid[93],concrete,"topstl1_lvBLK",0,0,0);
  blockPos=G4ThreeVector(0,in*m*2+Al_height/2.0+bb_height*2.0+18*in*m/2.0,-NumiData->HadrBox_length/2.0+4*bb_width/2.0);
  new G4PVPlacement(0,blockPos,"topstl1_BLK",BLK_log[93],HadrBox,false,0);

  BLK_solid[94]= new G4Box("topstl2_sBLK",conc_width/2.0,18*in*m/2.0,conc_length/2.0);
  BLK_log[94]=new G4LogicalVolume(BLK_solid[94],concrete,"topstl2_lvBLK",0,0,0);
  blockPos=G4ThreeVector(-81*in*m-conc_width/2.0,-in*m-Al_height/2.0-4*bb_height+conc_height+18*in*m/2.0,-NumiData->HadrBox_length/2.0+conc_length/2.0);
  new G4PVPlacement(0,blockPos,"topstl2_BLK",BLK_log[94],HadrBox,false,0);

  BLK_solid[95]= new G4Box("topstl3_sBLK",conc_width/2.0,18*in*m/2.0,conc_length/2.0);
  BLK_log[95]=new G4LogicalVolume(BLK_solid[95],concrete,"topstl3_lvBLK",0,0,0);
  blockPos=G4ThreeVector(81*in*m+conc_width/2.0,-in*m-Al_height/2.0-4*bb_height+conc_height+18*in*m/2.0,-NumiData->HadrBox_length/2.0+conc_length/2.0);
  new G4PVPlacement(0,blockPos,"topstl3_BLK",BLK_log[95],HadrBox,false,0);

  G4cout << "Hadron Absorber Constructed" << G4endl;

  // Constructs a volume in which to put the pre-absorber concrete shielding.
  //-------------------------------------------------------------------------

  G4Box*  expShldBox = new G4Box("expShldBox",NumiData->HadrBox_width/2.,NumiData->HadrBox_height/2.0,CShld_length/2.0);
  ShldBox_log =new G4LogicalVolume(expShldBox,Air,"LVShldBox",0,0,0);
  G4ThreeVector shldpos=G4ThreeVector(0,y_shld,z_shld)-tunnelPos;
  nam="ConcShield";
  ShldBox=new G4PVPlacement(zrot,shldpos,nam,ShldBox_log,pvTUNE,false,0);
  
  CShld_solid[0] = new G4Box("CShld_sBLK0",(A_width+B_width)/2.0,A_height/2.0,2*A_length/2.0);
  CShld_log[0] = new G4LogicalVolume(CShld_solid[0],concrete,"CShld_lvBLK0",0,0,0);
  blockPos = G4ThreeVector(1/2.0*A_height+(A_width+B_width)/2.0,-NumiData->HadrBox_height/2.0+B_height/2.0,-CShld_length/2.0+B_length);
  new G4PVPlacement(0,blockPos,"CShld_BLK0",CShld_log[0],ShldBox,false,0);
  CShld_solid[1] = new G4Box("CShld_sBLK1",(A_width+B_width)/2.0,A_height/2.0,2*A_length/2.0);
  CShld_log[1] = new G4LogicalVolume(CShld_solid[1],concrete,"CShld_lvBLK1",0,0,0);
  blockPos = G4ThreeVector(-1/2.0*A_height-(A_width+B_width)/2.0,-NumiData->HadrBox_height/2.0+B_height/2.0,-CShld_length/2.0+B_length);
  new G4PVPlacement(0,blockPos,"CShld_BLK1",CShld_log[1],ShldBox,false,0);

  CShld_solid[2] = new G4Box("CShld_sBLK2",(A_width+B_width)/2.0,A_height/2.0,2*A_length/2.0);
  CShld_log[2] = new G4LogicalVolume(CShld_solid[2],concrete,"CShld_lvBLK2",0,0,0);
  blockPos = G4ThreeVector(-1/2.0*A_height-(A_width+B_width)/2.0,-NumiData->HadrBox_height/2.0+1.5*B_height,-CShld_length/2.0+B_length);
  new G4PVPlacement(0,blockPos,"CShld_BLK2",CShld_log[2],ShldBox,false,0);

  CShld_solid[3] = new G4Box("CShld_sBLK3",(A_width+B_width)/2.0,A_height/2.0,A_length/2.0);
  CShld_log[3] = new G4LogicalVolume(CShld_solid[3],concrete,"CShld_lvBLK3",0,0,0);
  blockPos = G4ThreeVector(1/2.0*A_height+(A_width+B_width)/2.0,-NumiData->HadrBox_height/2.0+1.5*B_height,-CShld_length/2.0+B_length/2.0);
  new G4PVPlacement(0,blockPos,"CShld_BLK3",CShld_log[3],ShldBox,false,0);

  CShld_solid[4] = new G4Box("CShld_sBLK4",B_width/2.0,B_height/2.0,B_length/2.0);
  CShld_log[4] = new G4LogicalVolume(CShld_solid[4],concrete,"CShld_lvBLK4",0,0,0);
  blockPos = G4ThreeVector(1/2.0*A_height+B_width/2.0,-NumiData->HadrBox_height/2.0+1.5*B_height,-CShld_length/2.0+1.5*B_length);
  new G4PVPlacement(0,blockPos,"CShld_BLK4",CShld_log[4],ShldBox,false,0);

  CShld_solid[5] = new G4Box("CShld_sBLK5",C_width/2.0,C_height/2.0,C_length/2.0);
  CShld_log[5] = new G4LogicalVolume(CShld_solid[5],concrete,"CShld_lvBLK5",0,0,0);
  blockPos = G4ThreeVector(1/2.0*D_width+1/2.0*C_width,-NumiData->HadrBox_height/2.0+C_height/2.0,-CShld_length/2.0+C_length/2.0);
  new G4PVPlacement(0,blockPos,"CShld_BLK5",CShld_log[5],ShldBox,false,0);

  CShld_solid[6] = new G4Box("CShld_sBLK6",C_width/2.0,C_height/2.0,C_length/2.0);
  CShld_log[6] = new G4LogicalVolume(CShld_solid[6],concrete,"CShld_lvBLK6",0,0,0);
  blockPos = G4ThreeVector(-1/2.0*D_width-1/2.0*C_width,-NumiData->HadrBox_height/2.0+C_height/2.0,-CShld_length/2.0+C_length/2.0);
  new G4PVPlacement(0,blockPos,"CShld_BLK6",CShld_log[6],ShldBox,false,0);

  CShld_solid[7] = new G4Box("CShld_sBLK7",C_width/2.0,C_height/2.0,C_length/2.0);
  CShld_log[7] = new G4LogicalVolume(CShld_solid[7],concrete,"CShld_lvBLK7",0,0,0);
  blockPos = G4ThreeVector(1/2.0*D_width+1/2.0*C_width,-NumiData->HadrBox_height/2.0+C_height+A_width+C_height/2.0,-CShld_length/2.0+C_length/2.0);
  new G4PVPlacement(0,blockPos,"CShld_BLK7",CShld_log[7],ShldBox,false,0);

  CShld_solid[8] = new G4Box("CShld_sBLK8",C_width/2.0,C_height/2.0,C_length/2.0);
  CShld_log[8] = new G4LogicalVolume(CShld_solid[8],concrete,"CShld_lvBLK8",0,0,0);
  blockPos = G4ThreeVector(-1/2.0*D_width-1/2.0*C_width,-NumiData->HadrBox_height/2.0+C_height+A_width+C_height/2.0,-CShld_length/2.0+C_length/2.0);
  new G4PVPlacement(0,blockPos,"CShld_BLK8",CShld_log[8],ShldBox,false,0);

  CShld_solid[9] = new G4Box("CShld_sBLK9",A_height/2.0,A_width/2.0,A_length);
  CShld_log[9] = new G4LogicalVolume(CShld_solid[9],concrete,"CShld_lvBLK9",0,0,0);
  blockPos = G4ThreeVector(0,-NumiData->HadrBox_height/2.0+C_height+A_width/2.0,-CShld_length/2.0+A_length);
  new G4PVPlacement(0,blockPos,"CShld_BLK9",CShld_log[9],ShldBox,false,0);

  CShld_solid[10] = new G4Box("CShld_sBLK10",D_width/2.0,D_height/2.0,D_length/2.0);
  CShld_log[10] = new G4LogicalVolume(CShld_solid[10],concrete,"CShld_lvBLK10",0,0,0);
  blockPos = G4ThreeVector(0,-NumiData->HadrBox_height/2.0+D_height/2.0,-CShld_length/2.0+D_length/2.0);
  new G4PVPlacement(0,blockPos,"CShld_BLK10",CShld_log[10],ShldBox,false,0);

  CShld_solid[11] = new G4Box("CShld_sBLK11",D_width/2.0,D_height/2.0,D_length/2.0);
  CShld_log[11] = new G4LogicalVolume(CShld_solid[11],concrete,"CShld_lvBLK11",0,0,0);
  blockPos = G4ThreeVector(0,-NumiData->HadrBox_height/2.0+D_height+A_width+D_height/2.0,-CShld_length/2.0+D_length/2.0);
  new G4PVPlacement(0,blockPos,"CShld_BLK11",CShld_log[11],ShldBox,false,0);

  CShld_solid[12] = new G4Box("CShld_sBLK12",J_width/2.0,2*J_height/2.0,2*J_length/2.0);
  CShld_log[12] = new G4LogicalVolume(CShld_solid[12],concrete,"CShld_lvBLK12",0,0,0);
  blockPos = G4ThreeVector(0,-NumiData->HadrBox_height/2.0+2*B_height+J_height+3*in*m,-CShld_length/2.0+D_length/2.0);
  new G4PVPlacement(0,blockPos,"CShld_BLK12",CShld_log[12],ShldBox,false,0);

  // This is to create a block of air near the end of the beam pipe to compensate for a smaller block position.
  CShld_solid[15] = new G4Box("air_sBLK15",1.5*12*in*m/2.0,1.5*12*in*m/2.0,1.5*12*in*m/2.0);
  CShld_log[15] = new G4LogicalVolume(CShld_solid[15],Air,"air_lvBLK15",0,0,0);
  blockPos = G4ThreeVector((A_width+B_width)/2.0-1.5*12*in*m/2.0,A_height/2.0-1.5*12*in*m/2.0,-2*A_length/2.0+1.5*12*in*m/2.0);
  new G4PVPlacement(0,blockPos,CShld_log[15],"air_BLK15",CShld_log[2],false,0);

  G4cout << "Decay Pipe end cap Shielding Constructed" << G4endl;

}
