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

void NumiDetectorConstruction::ConstructTarget()
{
  // Mother volume (TGAR) position
  G4ThreeVector target_hall_position=G4ThreeVector(0,0,NumiData->TargetAreaLength/2.+NumiData->TargetAreaZ0);
  
  //------------------------------------------------------------------------------------------------
  //Target
  //------------------------------------------------------------------------------------------------
  G4RotationMatrix rotation;
  G4ThreeVector translation; 
  G4String vol_name;
  char no[3];
  /*
  //Create mother volume for whole target
  G4Tubs* Target_volume1=new G4Tubs("Trgt_vol1",0.,122.*mm,220.*mm,0,360.*deg);
  G4Tubs* Target_volume2=new G4Tubs("Trgt_vol2",0.,16.*mm,500.*mm,0,360.*deg);
  translation=G4ThreeVector(0.*mm,0.*mm,700.*mm);
  G4UnionSolid* Target_volume=new G4UnionSolid("Target_volume",Target_volume1,Target_volume2,G4Transform3D(rotation,translation));
  TRGT_lv=new G4LogicalVolume(Target_volume,Air,"TRGT_lv",0,0,0);
  TRGT_lv->SetVisAttributes(G4VisAttributes::Invisible);

  G4ThreeVector target_vol_position =G4ThreeVector(0.,0.,-(220.-52.+32.2)*mm)+
    G4ThreeVector(NumiData->TargetX0,NumiData->TargetY0,NumiData->TargetZ0)-target_hall_position;
  TRGT=new G4PVPlacement(0,target_vol_position,"TRGT",TRGT_lv,TGAR,false,0);

  //Be Window
  G4Tubs* Window= new G4Tubs("Window",0.,21.5*mm,1.*mm,0,360.*deg);
  G4LogicalVolume* BeWindow_lv=new G4LogicalVolume(Window,Be,"BeWindow_lv",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,(220.-52.-24.-23.-2*127.5-18.-7.)*mm),"BEWN",BeWindow_lv,TRGT,false,0);

  //Canister front part
  G4Tubs* Canister_front1= new G4Tubs("Canister_front1",17.5*mm,106.*mm,9.*mm,0,360.*deg);
  G4Tubs* Canister_front_ind= new G4Tubs("Canister_ind",36.*mm,78.*mm,6.*mm,0,360.*deg);
  translation=G4ThreeVector(0.*mm,0.*mm,-10.*mm);
  G4SubtractionSolid* Canister_front=
    new G4SubtractionSolid("Canister_front",Canister_front1,Canister_front_ind,G4Transform3D(rotation,translation));
  G4Tubs* screw_hole=new G4Tubs("screw_hole",0,3.01*mm,50.*mm,0,360.*deg);
  for (G4int ii=0;ii<6;ii++){
    translation=G4ThreeVector(sin(ii*60.*deg)*29.*mm,cos(ii*60.*deg)*29.*mm,-45.*mm);
    Canister_front=new G4SubtractionSolid("Canister_front",Canister_front,screw_hole,G4Transform3D(rotation,translation));
  }

  translation=G4ThreeVector(0.*mm,0.*mm,-9.*mm);
  Canister_front=new G4SubtractionSolid("Canister_front",Canister_front,Window,G4Transform3D(rotation,translation));
  G4LogicalVolume* Canister_log=new G4LogicalVolume(Canister_front,Fe,"Canister_log",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,(220.-52.-24.-23.-2*127.5-9.)*mm),"CNST",Canister_log,TRGT,false,0); 

  //Front cover
  G4Tubs* Front_cover1=new G4Tubs("Front_cover",17.5*mm,36.*mm,6.*mm,0.,360.*deg);
  G4Tubs* Front_cover_ind=new G4Tubs("Front_cover",0.*mm,21.5*mm,5.*mm,0.,360.*deg);
  translation=G4ThreeVector(0.*mm,0.*mm,-2.5*mm);
  G4SubtractionSolid* Front_cover=new G4SubtractionSolid("Front_cover",Front_cover1,Front_cover_ind,G4Transform3D(rotation,translation));
  translation=G4ThreeVector(0.*mm,0.*mm,3.*mm);
  Front_cover=new G4SubtractionSolid("Front_cover",Front_cover,Window,G4Transform3D(rotation,translation));
  for (G4int ii=0;ii<6;ii++){
    translation=G4ThreeVector(sin(ii*60.*deg)*29.*mm,cos(ii*60.*deg)*29.*mm,-40.*mm);
    Front_cover=new G4SubtractionSolid("Front_cover",Front_cover,screw_hole,G4Transform3D(rotation,translation));
  }
  G4LogicalVolume* Front_cover_lv=new G4LogicalVolume(Front_cover,Fe,"Front_cover_lv",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,(220.-52.-24.-23.-2*127.5-18.-6.)*mm),"FRTC",Front_cover_lv,TRGT,false,0);

  //Canister Body
  G4Trap* Canister_body_left=new G4Trap("Canister_body_left",127.5*mm,127.5*mm,34.5*mm,83.5*mm,24.5*mm);
  G4Trap* Canister_body_right=new G4Trap("Canister_body_right",127.5*mm,127.5*mm,83.5*mm,51.*mm,26.5*mm);
  G4Box* Canister_body_mid=new G4Box("Canister_body_mid",37.*mm,83.5*mm,127.5*mm);
  rotation.rotateY(90.*deg);
  translation=G4ThreeVector(-61.4*mm,0.*mm,0.*mm);
  G4UnionSolid* Canister_body=new G4UnionSolid("Canister_body",Canister_body_mid,Canister_body_left,G4Transform3D(rotation,translation));
  translation=G4ThreeVector(63.4*mm,0.*mm,0.*mm);
  Canister_body=new G4UnionSolid("Canister_body",Canister_body,Canister_body_right,G4Transform3D(rotation,translation));
  
  G4Tubs* CB_hole=new G4Tubs("CB_hole",0.,77.5*mm,128.5*mm,0.,360.*deg);
  rotation=G4RotationMatrix(0.,0.,0.);
  translation=G4ThreeVector(-2.*mm,0.*mm,0.*mm);
  G4SubtractionSolid* CB=new G4SubtractionSolid("CB1",Canister_body,CB_hole,G4Transform3D(rotation,translation));
 
  G4Tubs* BudalM_hole=new G4Tubs("BudalM_hole",0.,32.5*mm,10.*cm,0.,360.*deg);
  rotation=G4RotationMatrix(0.*deg,90.*deg,90.*deg);
  translation=G4ThreeVector(75.*mm,0.*mm,50.*mm);
  CB=new G4SubtractionSolid("CB2",CB,BudalM_hole,G4Transform3D(rotation,translation));

  G4LogicalVolume* Canister_body_lv=new G4LogicalVolume(CB,Fe,"Canister_body_lv",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(2.*mm,0.,(220.-52.-24.-23.-127.5)*mm),"CBDY",Canister_body_lv,TRGT,false,0);
  
  // Canister end
  G4Tubs* Canister_end_tube= new G4Tubs("Canister_end_tube",77.5*mm,120.*mm,2.*mm,0.,360.*deg);           
  G4Cons* Canister_end_cone= new G4Cons("Canister_end_cone",77.5*mm,120.*mm,77.5*mm,111.*mm,9.5*mm,0.,360.*deg);
  rotation=G4RotationMatrix(0.,0.,0.); 
  translation=G4ThreeVector(0.,0.,11.5*mm); 
  G4UnionSolid* Canister_end1= new G4UnionSolid("Canister_end1",Canister_end_tube,Canister_end_cone,G4Transform3D(rotation,translation));
  G4Tubs* screw_hole_M8=new G4Tubs("screw_hole",0,4.01*mm,50.*mm,0,360.*deg);
  translation=G4ThreeVector(sin(22.5*deg)*94.*mm,cos(22.5*deg)*94.*mm,0.*mm);
  G4SubtractionSolid* Canister_end=new G4SubtractionSolid("Canister_end",Canister_end1,screw_hole_M8,G4Transform3D(rotation,translation));
  for (G4int ii=2;ii<17;ii++){
    translation=G4ThreeVector(sin(ii*22.5*deg)*94.*mm,cos(ii*22.5*deg)*94.*mm,0.*mm);
    Canister_end=new G4SubtractionSolid("Canister_end",Canister_end,screw_hole_M8,G4Transform3D(rotation,translation));
  }
  
  G4LogicalVolume* Canister_end_lv=new G4LogicalVolume(Canister_end,Fe,"Canister_end_lv",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,(220.-52.-24.-21.)*mm),"CEND",Canister_end_lv,TRGT,false,0);
  
  //Back cover
  G4Tubs* Back_cover1=new G4Tubs("Back_cover",30.*mm,106.*mm,12.*mm,0.,360.*deg);
  G4Tubs* Pipe_hole=new G4Tubs("Pipe_hole",0.*mm,5.*mm,50.*mm,0.,360.*deg);

  translation=G4ThreeVector(0.,60.*mm,0.);
  G4SubtractionSolid* Back_cover=new G4SubtractionSolid("Back_cover",Back_cover1,Pipe_hole,G4Transform3D(rotation,translation));
  translation=G4ThreeVector(0.,-60.*mm,0.);
  Back_cover=new G4SubtractionSolid("Back_cover",Back_cover,Pipe_hole,G4Transform3D(rotation,translation));


  G4Tubs* BC_Inner_hole=new G4Tubs("BC_Inner_hole",0.,74.*mm,8.*mm,0.,360.*deg);
  translation=G4ThreeVector(0.,0.,-12.5*mm);
  Back_cover=new G4SubtractionSolid("Back_cover",Back_cover,BC_Inner_hole,G4Transform3D(rotation,translation));

  for (G4int ii=1;ii<17;ii++){
    translation=G4ThreeVector(sin(ii*22.5*deg)*94.*mm,cos(ii*22.5*deg)*94.*mm,0.*mm);
    Back_cover=new G4SubtractionSolid("Back_cover",Back_cover,screw_hole_M8,G4Transform3D(rotation,translation));
  }
  G4Tubs* screw_hole_M4=new G4Tubs("screw_hole_M4",0.,4.01*mm,10.*mm,0.,360.*deg);
  for (G4int ii=0;ii<4;ii++){
    translation=G4ThreeVector(sin(ii*90.*deg)*35.*mm,cos(ii*90.*deg)*35.*mm,12.*mm);
    Back_cover=new G4SubtractionSolid("Back_cover",Back_cover,screw_hole_M4,G4Transform3D(rotation,translation));
  }
  
  G4LogicalVolume* Back_cover_lv=new G4LogicalVolume(Back_cover,Fe,"Back_cover_lv",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,(220.-52.-12.)*mm),"ENDC",Back_cover_lv,TRGT,false,0);
  
  //screws M6
  G4double sside=6.*mm;
  G4Trap* screw_M6_1=new G4Trap("screw_M6_1",1.*mm,1.*mm,sside/2,sside,0.86603*sside/2.);
  G4Trap* screw_M6_2=new G4Trap("screw_M6_2",1.*mm,1.*mm,sside,sside/2,0.86603*sside/2.);
  translation=G4ThreeVector(0*mm,0.*mm,0.86603*sside);
  rotation=G4RotationMatrix(0.,0.,0.); 
  G4UnionSolid* head_M6=new G4UnionSolid("head_M6",screw_M6_1,screw_M6_2,G4Transform3D(rotation,translation));

  G4Tubs* screw_M6_3=new G4Tubs("scrw_M6_3",0.,3.*mm,11.5*mm,0,360.*deg);
  rotation.rotateX(0.*deg);  rotation.rotateY(90.*deg);  rotation.rotateZ(0.*deg); 
  translation=G4ThreeVector(-0.86603*sside/2,0.,-11.5*mm);
  G4UnionSolid* screw=new G4UnionSolid("screw",screw_M6_3,head_M6,G4Transform3D(rotation,translation));
  G4LogicalVolume* screw_M6_log=new G4LogicalVolume(screw,Fe,"screw_log",0,0,0);
    rotation=G4RotationMatrix(0.,0.,0.); 
  G4AssemblyVolume* front_screws=new G4AssemblyVolume();
  G4Transform3D trans3d;
  for (G4int ii=0;ii<6;ii++){
    translation=G4ThreeVector(sin(ii*60.*deg)*29.*mm,cos(ii*60.*deg)*29.*mm,0);
    trans3d=G4Transform3D(rotation,translation);
    front_screws->AddPlacedVolume(screw_M6_log,trans3d);
  }
  translation=G4ThreeVector(0.*mm,0.*mm,(220.-52.-24.-23.-2*127.5-18.-6.+5.5-2)*mm);
  trans3d=G4Transform3D(rotation,translation);
  front_screws->MakeImprint(TRGT_lv,trans3d,1);
  
  //screws M8
  sside=8.*mm;
  G4Trap* screw_M8_1=new G4Trap("screw_M8_1",1.*mm,1.*mm,sside/2,sside,0.86603*sside/2.);
  G4Trap* screw_M8_2=new G4Trap("screw_M8_2",1.*mm,1.*mm,sside,sside/2,0.86603*sside/2.);
  translation=G4ThreeVector(0*mm,0.*mm,0.86603*sside);
  rotation=G4RotationMatrix(0.,0.,0.); 
  G4UnionSolid* head_M8=new G4UnionSolid("head_M8",screw_M8_1,screw_M8_2,G4Transform3D(rotation,translation));

  G4Tubs* screw_M8_3=new G4Tubs("scrw_M8_3",0.,4.*mm,11.5*mm,0,360.*deg);
  rotation.rotateX(0.*deg);  rotation.rotateY(90.*deg);  rotation.rotateZ(0.*deg); 
  translation=G4ThreeVector(-0.86603*sside/2,0.,11.5*mm);
  G4UnionSolid* screw_M8=new G4UnionSolid("screw",screw_M8_3,head_M8,G4Transform3D(rotation,translation));
  G4LogicalVolume* screw_M8_log=new G4LogicalVolume(screw_M8,Fe,"screw_M8_log",0,0,0);
  rotation=G4RotationMatrix(0.,0.,0.); 
  G4AssemblyVolume* end_screws=new G4AssemblyVolume();
  for (G4int ii=1;ii<17;ii++){
    translation=G4ThreeVector(sin(ii*22.5*deg)*94.*mm,cos(ii*22.5*deg)*94.*mm,0);
    trans3d=G4Transform3D(rotation,translation);
    end_screws->AddPlacedVolume(screw_M8_log,trans3d);
  }
  translation=G4ThreeVector(0.*mm,0.*mm,(220.-52.-11.5+2)*mm);
  trans3d=G4Transform3D(rotation,translation);
  end_screws->MakeImprint(TRGT_lv,trans3d,1);
*/
  // Carbon Target
  for (G4int ii=0;ii<NumiData->TargetSegmentNo;ii++){
    
    G4double TGT_w=NumiData->TargetSWidth/2.;
    G4double TGT_h=NumiData->TargetSHeight/2.;
    G4double TGT_l=NumiData->TargetSLength/2.; //these are 0.5*true dimension

    //G4ThreeVector target_position=G4ThreeVector(0.,0.,TGT_l+NumiData->TargetZ0+(NumiData->TargetSegmentPitch+TGT_l*2)*ii)
    //-target_hall_position;
    G4VSolid* TGT1_solid;
    if (NumiData->TargetEndRounded) {
      TGT_l=TGT_l-NumiData->TargetSWidth/2.;
      G4Box* TGT2_solid=new G4Box("TGT2_solid",TGT_w,TGT_h,TGT_l);
      G4Tubs* TGTRE_solid=new G4Tubs("TGTRE_solid",0.,TGT_w,TGT_h,0.,360.*deg); 
      rotation=G4RotationMatrix(0,0,0); rotation.rotateX(90.*deg);
      translation=G4ThreeVector(0,0,-TGT_l);
      TGT1_solid=new G4UnionSolid("TGT1_solid",TGT2_solid,TGTRE_solid,G4Transform3D(rotation,translation));
      translation=G4ThreeVector(0,0,TGT_l);
      TGT1_solid=new G4UnionSolid("TGT1_solid",TGT1_solid,TGTRE_solid,G4Transform3D(rotation,translation));
    }
    else{
      TGT1_solid=new G4Box("TGT1_solid",TGT_w,TGT_h,TGT_l);
    }
    G4Tubs* CPG_solid=new G4Tubs("CPG_solid",0.,NumiData->TargetCPGRadius,2.*TGT_l,0.,360.*deg);
    rotation=G4RotationMatrix(0,0,0);
    translation=G4ThreeVector(0,NumiData->TargetCPGPosition,0);
    G4SubtractionSolid* TGT_solid=new G4SubtractionSolid("TGT_solid",TGT1_solid,CPG_solid,G4Transform3D(rotation,translation));
    translation=G4ThreeVector(0,-NumiData->TargetCPGPosition,0);
    TGT_solid=new G4SubtractionSolid("TGT_solid",TGT_solid,CPG_solid,G4Transform3D(rotation,translation));
    G4LogicalVolume* TGT1_log=new G4LogicalVolume(TGT_solid,Target,"TGT1_log",0,0,0);

    rotation=G4RotationMatrix(0,0,0);
    rotation.rotateX(atan(NumiData->TargetDydz));
    rotation.rotateY(atan(NumiData->TargetDxdz));
    // with this translation rotation axis is at the begining of the volume (x0,y0,z0) and not at its center
    translation=G4ThreeVector(-(sin(atan(NumiData->TargetDxdz)))*cos(atan(NumiData->TargetDydz))*NumiData->TargetSLength/2.,(sin(atan(NumiData->TargetDydz)))*NumiData->TargetSLength/2.,(1-cos(atan(NumiData->TargetDxdz))*cos(atan(NumiData->TargetDydz)))*NumiData->TargetSLength/2.);
    G4ThreeVector target_position=G4ThreeVector(NumiData->TargetX0,NumiData->TargetY0,NumiData->TargetZ0+NumiData->TargetSLength/2.+(NumiData->TargetSegmentPitch+NumiData->TargetSLength)*ii)-target_hall_position-translation;
    new G4PVPlacement(G4Transform3D(rotation,target_position),"TGT1",TGT1_log,TGAR,false,0);

  }
  //Container
  for (G4int ii=0;ii<NumiData->NContainerN;ii++){

    vol_name=NumiData->CTubeVolName[ii];

    G4VSolid* Cont_solid;
    Cont_solid=new G4Tubs("Cont_solid",NumiData->CTubeRin[ii],NumiData->CTubeRout[ii],NumiData->CTubeLength[ii]/2.,0.,360.*deg);
    G4ThreeVector Cont_position=G4ThreeVector(0.,0.,NumiData->CTubeZ0[ii]+NumiData->CTubeLength[ii]/2.)-target_hall_position;
    
    if (NumiData->CTubeHoleIndex[ii]==1){
      G4Tubs* CP_solid=new G4Tubs("CP_solid",0.,NumiData->CPipeRadiusOut[0],2*NumiData->CTubeLength[ii],0.,360.*deg);
      rotation=G4RotationMatrix(0,0,0);
      translation=G4ThreeVector(0,NumiData->CPipeY0[0],0);
      Cont_solid = new G4SubtractionSolid("Cont_solid",Cont_solid,CP_solid,G4Transform3D(rotation,translation));
      translation=G4ThreeVector(0,-NumiData->CPipeY0[0],0);
      Cont_solid = new G4SubtractionSolid("Cont_solid",Cont_solid,CP_solid,G4Transform3D(rotation,translation));
    
      // clip off edges so that target fits inside
      G4double TGT_w=NumiData->TargetSWidth/2.;
      G4double TGT_h=NumiData->TargetSHeight/2.;
      G4double TGT_l=NumiData->TargetSLength/2.; //these are 0.5*true dimension

      G4Box* TGT_box_solid=new G4Box("TGT1_solid",1.1*TGT_w,1.025*TGT_h,TGT_l);
      Cont_solid = new G4SubtractionSolid("edge",Cont_solid,TGT_box_solid);
    }
    G4LogicalVolume* ContTube_lv=new G4LogicalVolume(Cont_solid,GetMaterial(NumiData->CTubeGeantMat[ii]),"CONT_lv",0,0,0);
    CNT[ii]=new G4PVPlacement(0,Cont_position,vol_name,ContTube_lv,TGAR,false,0);
  }

  //Cooling pipes
  for (G4int ii=0;ii<NumiData->NCPipeN;ii++){
    sprintf(no,"%d",ii+1);
    vol_name="CPI";
    vol_name.append(no);

    if (NumiData->CPipeCurvRad[ii]==0){
      G4Tubs* CP_solid=new G4Tubs("CP_solid",NumiData->CPipeRadiusIn[ii],NumiData->CPipeRadiusOut[ii],NumiData->CPipeLength[ii]/2.,0.,360.*deg);
      G4Tubs* CP_water_solid=new 
	G4Tubs("CP_water_solid",0.,NumiData->CPipeRadiusOut[ii]-NumiData->CPipeWallThick[ii],NumiData->CPipeLength[ii]/2.,0.,360.*deg);
      G4LogicalVolume* CPipe_lv=new G4LogicalVolume(CP_solid,GetMaterial(NumiData->CPGeantMat[ii]),"CPIP_lv",0,0,0);
      
      rotation=G4RotationMatrix(0,0,0);
      rotation.rotateX(atan(NumiData->CPipeDYDZ[ii]));
      rotation.rotateY(atan(NumiData->CPipeDXDZ[ii]));
      translation=G4ThreeVector(-(sin(atan(NumiData->CPipeDXDZ[ii])))*cos(atan(NumiData->CPipeDYDZ[ii]))*NumiData->CPipeLength[ii]/2.,(sin(atan(NumiData->CPipeDYDZ[ii])))*NumiData->CPipeLength[ii]/2.,(1-cos(atan(NumiData->CPipeDXDZ[ii]))*cos(atan(NumiData->CPipeDYDZ[ii])))*NumiData->CPipeLength[ii]/2.);
      G4ThreeVector CPipe_position=G4ThreeVector(NumiData->CPipeX0[ii],NumiData->CPipeY0[ii],NumiData->CPipeZ0[ii]+NumiData->CPipeLength[ii]/2.)-target_hall_position-translation;

      CPIP[ii]=new G4PVPlacement(G4Transform3D(rotation,CPipe_position),vol_name,CPipe_lv,TGAR,false,0);
 
      if (NumiData->CPipeFilledWater[ii]) {
	G4Tubs* CP_water_solid=new 
	G4Tubs("CP_water_solid",NumiData->CPipeRadiusIn[ii]+NumiData->CPipeWallThick[ii],NumiData->CPipeRadiusOut[ii]-NumiData->CPipeWallThick[ii],NumiData->CPipeLength[ii]/2.,0.,360.*deg);
	G4LogicalVolume* CPipeW_lv=new G4LogicalVolume(CP_water_solid,Water,"CPIPW_lv",0,0,0);
	G4VisAttributes * WaterPipeAtt=new G4VisAttributes(G4Colour(0.,0.,1.));
	CPipeW_lv->SetVisAttributes(WaterPipeAtt);
	new G4PVPlacement(0,0,"CPIPW",CPipeW_lv,CPIP[ii],false,0);
      }
    }
    if (NumiData->CPipeCurvRad[ii]!=0){
      G4double deltaphi=0;
      if((-NumiData->CPipeOpenAng[ii]+NumiData->CPipeCloseAng[ii])>NumiData->CPipeOpenAng[ii]){
	deltaphi=-2*NumiData->CPipeOpenAng[ii]+NumiData->CPipeCloseAng[ii]; //something is strange with G4Torus
      }
      else {
	deltaphi=360.*deg-2*NumiData->CPipeOpenAng[ii]+NumiData->CPipeCloseAng[ii]; //something is strange with G4Torus
      }

      G4Torus* CPipe_tor=new G4Torus("CPipe_tor",NumiData->CPipeRadiusOut[ii]-NumiData->CPipeWallThick[ii],NumiData->CPipeRadiusOut[ii],NumiData->CPipeCurvRad[ii],NumiData->CPipeOpenAng[ii],deltaphi);
      G4LogicalVolume* CPipeTor_lv=new G4LogicalVolume(CPipe_tor,GetMaterial(NumiData->CPGeantMat[ii]),"CPPipeTor_lv",0,0,0);

      rotation=G4RotationMatrix(0,0,0);
      rotation.rotateX(atan(NumiData->CPipeDYDZ[ii]));
      rotation.rotateY(atan(NumiData->CPipeDXDZ[ii]));
      translation=G4ThreeVector(-(sin(atan(NumiData->CPipeDXDZ[ii])))*cos(atan(NumiData->CPipeDYDZ[ii]))*NumiData->CPipeLength[ii]/2.,(sin(atan(NumiData->CPipeDYDZ[ii])))*NumiData->CPipeLength[ii]/2.,(1-cos(atan(NumiData->CPipeDXDZ[ii]))*cos(atan(NumiData->CPipeDYDZ[ii])))*NumiData->CPipeLength[ii]/2.);
      G4ThreeVector CPipe_position=G4ThreeVector(NumiData->CPipeX0[ii],NumiData->CPipeY0[ii],NumiData->CPipeZ0[ii]+NumiData->CPipeLength[ii]/2.)-target_hall_position-translation;

      CPIP[ii]=new G4PVPlacement(G4Transform3D(rotation,CPipe_position),vol_name,CPipeTor_lv,TGAR,false,0);

      if (NumiData->CPipeFilledWater[ii]) {
	G4Torus* CPipeW_tor=new G4Torus("CPipe_tor",NumiData->CPipeRadiusIn[ii]+NumiData->CPipeWallThick[ii],NumiData->CPipeRadiusOut[ii]-NumiData->CPipeWallThick[ii],NumiData->CPipeCurvRad[ii],NumiData->CPipeOpenAng[ii],deltaphi);
	G4LogicalVolume* CPipeTorW_lv=new G4LogicalVolume(CPipeW_tor,Water,"CPPipeTorW_lv",0,0,0);
	G4VisAttributes * WaterPipeAtt=new G4VisAttributes(G4Colour(0.,0.,1.));
	CPipeTorW_lv->SetVisAttributes(WaterPipeAtt);
	new G4PVPlacement(0,0,"CPIPW",CPipeTorW_lv,CPIP[ii],false,0);
      }
    }

  }
  /*
  //Cooling pipes
  G4Torus* pipe_tor=new G4Torus("pipe_tor",0.,3.*mm,24.*mm,0.,180.*deg);
  G4Tubs* pipe_target=new G4Tubs("pipe_target",0.,3.*mm,TGT_l+163.5*mm,0,360.*deg);
  G4Tubs* pipe= new G4Tubs("pipe",0,3.*mm,157.*mm,0,360.*deg);
  rotation=G4RotationMatrix(0.*deg,90.*deg,0.);  
  translation=G4ThreeVector(24.*mm,-(TGT_l+163.49*mm),0.);
  G4UnionSolid* pipe1= new G4UnionSolid("pipe1",pipe_tor,pipe_target,G4Transform3D(rotation,translation));
  rotation=G4RotationMatrix(0.*deg,90.*deg,0.);  
  translation=G4ThreeVector(-24.*mm,-156.99*mm,0.);
  pipe1= new G4UnionSolid("pipe1",pipe1,pipe,G4Transform3D(rotation,translation));
  G4LogicalVolume* pipe1_lv=new G4LogicalVolume(pipe1,Fe,"pipe1_lv",0,0,0);
 
  //Water
  G4Torus* pipe_tor_w=new G4Torus("pipe_tor_w",0.,2.5*mm,24.*mm,0.,180.*deg);
  G4Tubs* pipe_target_w=new G4Tubs("pipe_target_w",0.,2.5*mm,TGT_l+163.5*mm,0,360.*deg);
  G4Tubs* pipe_w= new G4Tubs("pipe",0,2.5*mm,157.*mm,0,360.*deg);
  rotation=G4RotationMatrix(0.*deg,90.*deg,0.);  
  translation=G4ThreeVector(24.*mm,-(TGT_l+163.49*mm),0.);
  G4UnionSolid* pipe1_w= new G4UnionSolid("pipe1_w",pipe_tor_w,pipe_target_w,G4Transform3D(rotation,translation));
  rotation=G4RotationMatrix(0.*deg,90.*deg,0.);  
  translation=G4ThreeVector(-24.*mm,-156.99*mm,0.);
  pipe1_w= new G4UnionSolid("pipe1_w",pipe1_w,pipe_w,G4Transform3D(rotation,translation));
  G4LogicalVolume* pipe1_w_lv=new G4LogicalVolume(pipe1_w,Water,"pipe1_w_lv",0,0,0);
  G4VisAttributes * WaterPipeAtt=new G4VisAttributes(G4Colour(0.,0.,1.));
  pipe1_w_lv->SetVisAttributes(WaterPipeAtt);
  
  rotation=G4RotationMatrix(0.*deg,90.*deg,90.*deg);  
  translation=G4ThreeVector(0.*mm,36.*mm,-(156.99+157.-220.)*mm);
  G4VPhysicalVolume* CPI1=new G4PVPlacement(G4Transform3D(rotation,translation),"CPI1",pipe1_lv,TRGT,false,0);

  rotation=G4RotationMatrix(0.*deg,0.*deg,0.*deg);  
  translation=G4ThreeVector(0.*mm,0.*mm,0.*mm);
  new G4PVPlacement(G4Transform3D(rotation,translation),"CPW1",pipe1_w_lv,CPI1,false,0);

  rotation=G4RotationMatrix(0.*deg,90.*deg,-90.*deg);  
  translation=G4ThreeVector(0.*mm,-36.*mm,-(156.99+157.-220.)*mm);
  new G4PVPlacement(G4Transform3D(rotation,translation),"CPI2",pipe1_lv,TRGT,false,0);
  */
  G4cout << "Target Constructed" << G4endl;
}
