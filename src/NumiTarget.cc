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
#include "NumiDataInput.hh"

void NumiDetectorConstruction::ConstructTarget()
{


  // Mother volume (TGAR) position
  G4ThreeVector target_hall_position=G4ThreeVector(0,0,NumiData->TargetAreaLength/2.+NumiData->TargetAreaZ0);
  G4ThreeVector targetPosition=G4ThreeVector(NumiData->TargetX0,NumiData->TargetY0,NumiData->TargetZ0);
  
  //Create Target Region
  G4Region* TargetRegion = new G4Region("Numi Target");

  //------------------------------------------------------------------------------------------------
  //Target
  //------------------------------------------------------------------------------------------------
  G4RotationMatrix rotation;
  G4ThreeVector translation; 
  G4String Sname,PVname,LVname,vol_name;
  char no[3];
 
  // Carbon Target
  
  //First create one segment
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
  G4LogicalVolume* LVTargetFin=new G4LogicalVolume(TGT_solid,Target,"LVTargetFin",0,0,0);
  TargetRegion->AddRootLogicalVolume(LVTargetFin);
  //Now create TargetSegmentNo of target fins and place them
  for (G4int ii=0;ii<NumiData->TargetSegmentNo;ii++){
    rotation=G4RotationMatrix(0,0,0);
    rotation.rotateX(atan(NumiData->TargetDydz));
    rotation.rotateY(atan(NumiData->TargetDxdz));
    // with this translation rotation axis is at the begining of the volume (x0,y0,z0) and not at its center
    translation=G4ThreeVector(-(sin(atan(NumiData->TargetDxdz)))*cos(atan(NumiData->TargetDydz))*NumiData->TargetSLength/2.,(sin(atan(NumiData->TargetDydz)))*NumiData->TargetSLength/2.,(1-cos(atan(NumiData->TargetDxdz))*cos(atan(NumiData->TargetDydz)))*NumiData->TargetSLength/2.);
    G4ThreeVector targetSegmentPosition=targetPosition+G4ThreeVector(0.,0.,NumiData->TargetSLength/2.+(NumiData->TargetSegmentPitch+NumiData->TargetSLength)*ii)-target_hall_position-translation;
    new G4PVPlacement(G4Transform3D(rotation,targetSegmentPosition),"TGT1",LVTargetFin,TGAR,false,0);
  }

  // Budal Monitor
  rotation=G4RotationMatrix(0,0,0);
  rotation.rotateX(atan(NumiData->BudalDydz));
  rotation.rotateY(atan(NumiData->BudalDxdz));
  rotation.rotateZ(90.*deg);
  // with this translation rotation axis is at the begining of the volume (x0,y0,z0) and not at its center
  translation=G4ThreeVector(-(sin(atan(NumiData->BudalDxdz)))*cos(atan(NumiData->BudalDydz))*NumiData->TargetSLength/2.,(sin(atan(NumiData->BudalDydz)))*NumiData->TargetSLength/2.,(1-cos(atan(NumiData->BudalDxdz))*cos(atan(NumiData->BudalDydz)))*NumiData->TargetSLength/2.);
  G4ThreeVector budalMonitorPosition=targetPosition+G4ThreeVector(0.,0.,NumiData->BudalZ0+NumiData->TargetSLength/2.)-target_hall_position-translation;
  G4LogicalVolume* LVBudalMonitor=new G4LogicalVolume(TGT1_solid,Target,"LVBudalMonitor",0,0,0);
  TargetRegion->AddRootLogicalVolume(LVBudalMonitor);
  new G4PVPlacement(G4Transform3D(rotation,budalMonitorPosition),"BudalMonitor",LVBudalMonitor,TGAR,false,0);

  //Container
  for (G4int ii=0;ii<NumiData->NContainerN;ii++){

    vol_name=NumiData->CTubeVolName[ii];
    G4VSolid* ContSolid;
    ContSolid=new G4Tubs("ContSolid",NumiData->CTubeRin[ii],NumiData->CTubeRout[ii],NumiData->CTubeLength[ii]/2.,0.,360.*deg);
    G4ThreeVector ContPosition=G4ThreeVector(0.,0.,NumiData->CTubeZ0[ii]+NumiData->CTubeLength[ii]/2.)-target_hall_position+targetPosition;    G4LogicalVolume* LVContTube=new G4LogicalVolume(ContSolid,GetMaterial(NumiData->CTubeGeantMat[ii]),vol_name.append("LV"),0,0,0);
    TargetRegion->AddRootLogicalVolume(LVContTube);
    CNT[ii]=new G4PVPlacement(0,ContPosition,vol_name,LVContTube,TGAR,false,0);
  }

  //Rings that hold target and cooling pipes
  for (G4int ii=0;ii<NumiData->NTgtRingN;ii++){
    vol_name=NumiData->TgtRingVolName[ii];

    G4VSolid* STgtRing;
    STgtRing = new G4Tubs("STgtRing",NumiData->TgtRingRin[ii],NumiData->TgtRingRout[ii],NumiData->TgtRingLength[ii]/2.,0.,360.*deg);
    G4ThreeVector TgtRingPosition=G4ThreeVector(0.,0.,NumiData->TgtRingZ0[ii]+NumiData->TgtRingLength[ii]/2.)-target_hall_position+targetPosition;
    //make holes for pipes
    G4Tubs* CP_solid=new G4Tubs("CP_solid",0.,NumiData->CPipeRadiusOut[0],2*NumiData->TgtRingLength[ii],0.,360.*deg);
    rotation=G4RotationMatrix(0,0,0);
    translation=G4ThreeVector(0,NumiData->CPipeY0[0],0);
    STgtRing = new G4SubtractionSolid("STgtRing",STgtRing,CP_solid,G4Transform3D(rotation,translation));
    translation=G4ThreeVector(0,-NumiData->CPipeY0[0],0);
    STgtRing = new G4SubtractionSolid("RingSolid",STgtRing,CP_solid,G4Transform3D(rotation,translation));
    // clip off edges so that target fits inside
    G4double TGT_w=NumiData->TargetSWidth/2.;
    G4double TGT_h=NumiData->TargetSHeight/2.;
    G4double TGT_l=NumiData->TargetSLength/2.; //these are 0.5*true dimension
    G4Box* TGT_box_solid=new G4Box("TgtBox",1.1*TGT_w,1.025*TGT_h,TGT_l);
    STgtRing = new G4SubtractionSolid("RingSolid",STgtRing,TGT_box_solid);
    
    G4LogicalVolume* LVTgtRing=new G4LogicalVolume(STgtRing,GetMaterial(NumiData->TgtRingGeantMaterial[ii]),vol_name.append("LV"),0,0,0);
    TargetRegion->AddRootLogicalVolume(LVTgtRing);
    new G4PVPlacement(0,TgtRingPosition,vol_name,LVTgtRing,TGAR,false,0);
  }

  //Cooling pipes
  for (G4int ii=0;ii<NumiData->NCPipeN;ii++){
    sprintf(no,"%d",ii+1);
    vol_name="PV";
    vol_name=vol_name.append(NumiData->CPipeVolName[ii]);
    Sname="S";
    Sname=Sname.append(NumiData->CPipeVolName[ii]);
    LVname="LV";
    LVname=LVname.append(NumiData->CPipeVolName[ii]);
    PVname="";
    PVname=PVname.append(NumiData->CPipeVolName[ii]);
 
    if (NumiData->CPipeCurvRad[ii]==0){
      G4Tubs* CP_solid=new G4Tubs(Sname,NumiData->CPipeRadiusIn[ii],NumiData->CPipeRadiusOut[ii],NumiData->CPipeLength[ii]/2.,0.,360.*deg);
      LVCPipe[ii]=new G4LogicalVolume(CP_solid,GetMaterial(NumiData->CPGeantMat[ii]),LVname,0,0,0);
      TargetRegion->AddRootLogicalVolume(LVCPipe[ii]);

      // Position the pipe
      rotation=G4RotationMatrix(0,0,0);
      rotation.rotateX(atan(NumiData->CPipeDYDZ[ii]));
      rotation.rotateY(atan(NumiData->CPipeDXDZ[ii]));
      translation=G4ThreeVector(-(sin(atan(NumiData->CPipeDXDZ[ii])))*cos(atan(NumiData->CPipeDYDZ[ii]))*NumiData->CPipeLength[ii]/2.,(sin(atan(NumiData->CPipeDYDZ[ii])))*NumiData->CPipeLength[ii]/2.,(1-cos(atan(NumiData->CPipeDXDZ[ii]))*cos(atan(NumiData->CPipeDYDZ[ii])))*NumiData->CPipeLength[ii]/2.);
      G4ThreeVector CPipe_position=G4ThreeVector(NumiData->CPipeX0[ii],NumiData->CPipeY0[ii],NumiData->CPipeZ0[ii]+NumiData->CPipeLength[ii]/2.)-target_hall_position-translation+targetPosition;

      PVCPipe[ii]=new G4PVPlacement(G4Transform3D(rotation,CPipe_position),PVname,LVCPipe[ii],TGAR,false,0);
 
      if (NumiData->CPipeFilledWater[ii]) {
	G4Tubs* CP_water_solid;
	Sname=Sname.append("_water");
	LVname=LVname.append("_water");
	PVname=PVname.append("_water");
	if (NumiData->CPipeRadiusIn[ii]==0) {
	  CP_water_solid=new 
	    G4Tubs(Sname,0.,NumiData->CPipeRadiusOut[ii]-NumiData->CPipeWallThick[ii],NumiData->CPipeLength[ii]/2.,0.,360.*deg);
	}
	else {
	  CP_water_solid=new 
	    G4Tubs(Sname,NumiData->CPipeRadiusIn[ii]+NumiData->CPipeWallThick[ii],NumiData->CPipeRadiusOut[ii]-NumiData->CPipeWallThick[ii],NumiData->CPipeLength[ii]/2.,0.,360.*deg);
	}
	LVCPipeW[ii]=new G4LogicalVolume(CP_water_solid,Water,LVname,0,0,0);
	G4VisAttributes * WaterPipeAtt=new G4VisAttributes(G4Colour(0.,0.,1.));
	LVCPipeW[ii]->SetVisAttributes(WaterPipeAtt);
	new G4PVPlacement(0,0,PVname,LVCPipeW[ii],PVCPipe[ii],false,0);
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
      G4Torus* CPipe_tor=new G4Torus(Sname,0.,NumiData->CPipeRadiusOut[ii],NumiData->CPipeCurvRad[ii],NumiData->CPipeOpenAng[ii],deltaphi);
      LVCPipe[ii]=new G4LogicalVolume(CPipe_tor,GetMaterial(NumiData->CPGeantMat[ii]),LVname,0,0,0);
      TargetRegion->AddRootLogicalVolume(LVCPipe[ii]);

      // Position the pipe  
      rotation=G4RotationMatrix(0,0,0);
      rotation.rotateX(atan(NumiData->CPipeDYDZ[ii]));
      rotation.rotateY(atan(NumiData->CPipeDXDZ[ii]));
      translation=G4ThreeVector(-(sin(atan(NumiData->CPipeDXDZ[ii])))*cos(atan(NumiData->CPipeDYDZ[ii]))*NumiData->CPipeLength[ii]/2.,(sin(atan(NumiData->CPipeDYDZ[ii])))*NumiData->CPipeLength[ii]/2.,(1-cos(atan(NumiData->CPipeDXDZ[ii]))*cos(atan(NumiData->CPipeDYDZ[ii])))*NumiData->CPipeLength[ii]/2.);
      G4ThreeVector CPipe_position=G4ThreeVector(NumiData->CPipeX0[ii],NumiData->CPipeY0[ii],NumiData->CPipeZ0[ii]+NumiData->CPipeLength[ii]/2.)-target_hall_position-translation+targetPosition;
     
      PVCPipe[ii]=new G4PVPlacement(G4Transform3D(rotation,CPipe_position),PVname,LVCPipe[ii],TGAR,false,0);
      
      if (NumiData->CPipeFilledWater[ii]) {
	Sname=Sname.append("_water");
	LVname=LVname.append("_water");
	PVname=PVname.append("_water");
	G4Torus* CPipeW_tor=new G4Torus(Sname,0.,NumiData->CPipeRadiusOut[ii]-NumiData->CPipeWallThick[ii],NumiData->CPipeCurvRad[ii],NumiData->CPipeOpenAng[ii],deltaphi);
	LVCPipeW[ii]=new G4LogicalVolume(CPipeW_tor,Water,LVname,0,0,0);
	G4VisAttributes * WaterPipeAtt=new G4VisAttributes(G4Colour(0.,0.,1.));
	LVCPipeW[ii]->SetVisAttributes(WaterPipeAtt);
	new G4PVPlacement(0,0,PVname,LVCPipeW[ii],PVCPipe[ii],false,0);
      }
    }
  }
  
  G4cout << "Target Constructed" << G4endl;
}
