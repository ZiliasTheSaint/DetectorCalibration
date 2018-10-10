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
//
// $Id: DetectorConstruction.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Dan Fulea, JUL 2013, National Institute of Public Health Bucharest, Cluj-Napoca Regional Centre, Romania 

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4RunManager.hh"


DetectorConstruction::DetectorConstruction()
: 
 solidWorldD(0),logicWorldD(0),physicalWorldD(0),

 solidUpperAir(0),logicUpperAir(0),physicalUpperAir(0),

 //solidDetector(0),
 logicDetector(0),physicalDetector(0),physicalHole(0),hole_diameter(0),hole_total_height(0),logicHole(0),
 fsolidWindow(0),fLogicWindow(0),fphysicalWindow(0),
 fsolidVacuum_top(0),flogicVacuum_top(0),fphysicalVacuum_top(0),
 fsolidVacuum_surface(0),flogicVacuum_surface(0),fphysicalVacuum_surface(0),
 fsolidVacuum_bottom(0),flogicVacuum_bottom(0),fphysicalVacuum_bottom(0),
 fsolidMonture_bottom(0),flogicMonture_bottom(0),fphysicalMonture_bottom(0),
 fsolidMonture_surface(0),flogicMonture_surface(0),fphysicalMonture_surface(0),

 fTargetMaterial(0),fWindowMaterial(0),fmonture_mat(0),fHoleMaterial(0),
 fgap_mat(0),fworld_material(0),

 factive_diameter(0),factive_height(0),
 fmonture_thickness(0),ftotal_diameter(0),ftotal_height(0),
 fwindow_thickness(0),fendCapToDetector(0),

 fUpperAir_height(0),fFrontalBeam_radius(0),fFrontalBeam_angle(0),

 logicSourceTopCylinder(0),solidSourceTopCylinder(0),physicalSourceTopCylinder(0),
 logicSourceEnvelope(0),solidSourceEnvelope(0),physicalSourceEnvelope(0),
 fsource_material(0),fsource_envelope_material(0),fsource_external_diameter(0),fsource_total_height(0),fsource_envelope_thickness(0),
 
 fsource_internal_diameter(0),fsource_upper_height(0),
 solidEnvelopeBottomCylinder(0),
     logicEnvelopeBottomCylinder(0),    
     physicalEnvelopeBottomCylinder(0), 
	 solidEnvelopeSurfaceCylinder(0),
     logicEnvelopeSurfaceCylinder(0),    
     physicalEnvelopeSurfaceCylinder(0), 
	 solidSourceSurfaceCylinder(0),
     logicSourceSurfaceCylinder(0),    
     physicalSourceSurfaceCylinder(0), 

 fCheckOverlaps(false)
{  
  // Option to switch on/off checking of volumes overlaps
  //fCheckOverlaps = true;  

  DefineMaterials();
  DefaultConstruction();
  detectorMessenger = new DetectorMessenger(this);
}
 
DetectorConstruction::~DetectorConstruction()
{ delete detectorMessenger;}

G4VPhysicalVolume* DetectorConstruction::Construct()
{	
	//ConstructDetector();
	return ConstructDetector();
}

void DetectorConstruction::DefineMaterials()
{ 
	G4String symbol;             //a=mass of a mole;
    G4double a, z, density;      //z=mean number of protons;  
    G4int iz, n;                 //iz=number of protons  in an isotope; 
                             // n=number of nucleons in an isotope;
	G4int ncomponents, natoms;
    G4double abundance, fractionmass;
	//---------------------
	G4NistManager* man = G4NistManager::Instance();
    G4Material* CO2_std = man->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
    G4Material* Ar_std = man->FindOrBuildMaterial("G4_Ar");
	//----------------------------
	G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
	G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
	G4Element* Ar  = new G4Element("Argon"  ,symbol="Ar" , z= 18., a= 39.95*g/mole);
	
	G4double pressure=600*0.00133322368*bar;
	
	G4Material* CO2 = 
    new G4Material("CarbonicGas", density= 1.842*mg/cm3, ncomponents=2, kStateGas, STP_Temperature, pressure);
    CO2->AddElement(C, natoms=1);
    CO2->AddElement(O, natoms=2);

	G4Material* Ar_gas=new G4Material("Ar_gas", density= 1.67*kg/m3, ncomponents=1, kStateGas, STP_Temperature, pressure);
	Ar_gas->AddElement(Ar, natoms=1);

	G4double rhoCO2=1.842*mg/cm3;
	G4double rhoAr=1.67*kg/m3;
	G4double fco2=0.2;
	G4double far=0.8;//fraction by volume in Ar-Co2 mixture.

	density=fco2*rhoCO2+far*rhoAr;//because is fraction by volume!!

	//now fraction by weight
	G4double wco2=fco2*rhoCO2/density;
	G4double war=far*rhoAr/density;

	G4Material* ArCO2_80_20 = 
		new G4Material("ArCO2_80_20",density, ncomponents=2, kStateGas, STP_Temperature, pressure);
		ArCO2_80_20->AddMaterial(Ar_gas, fractionmass=war);
		ArCO2_80_20->AddMaterial(CO2 , fractionmass=wco2);//20.*perCent);
    //-----------------
    G4Material* H = man->FindOrBuildMaterial("G4_H");
	G4Material* Carbon = man->FindOrBuildMaterial("G4_C");
	G4Material* N = man->FindOrBuildMaterial("G4_N");
	G4Material* Oxygen = man->FindOrBuildMaterial("G4_O");
	G4Material* Na = man->FindOrBuildMaterial("G4_Na");
	G4Material* Al = man->FindOrBuildMaterial("G4_Al");
	G4Material* Si = man->FindOrBuildMaterial("G4_Si");
	G4Material* K = man->FindOrBuildMaterial("G4_K");
	G4Material* Ca = man->FindOrBuildMaterial("G4_Ca");
	G4Material* Fe = man->FindOrBuildMaterial("G4_Fe");
    G4Material* soil_typicalloam_seltzer=
		new G4Material("soil_typicalloam_seltzer",1.5*g/cm3,ncomponents=10);
	    soil_typicalloam_seltzer->AddMaterial(H, fractionmass=0.028081);//PEGS4data->RHOZ
		soil_typicalloam_seltzer->AddMaterial(Carbon, fractionmass=0.144308);
		soil_typicalloam_seltzer->AddMaterial(N, fractionmass=1.0E-4);
		soil_typicalloam_seltzer->AddMaterial(Oxygen, fractionmass=0.496406);
		soil_typicalloam_seltzer->AddMaterial(Na, fractionmass=0.00816);
		soil_typicalloam_seltzer->AddMaterial(Al, fractionmass=0.089285);
		soil_typicalloam_seltzer->AddMaterial(Si, fractionmass=0.21315);
		soil_typicalloam_seltzer->AddMaterial(K, fractionmass=0.005562);
		soil_typicalloam_seltzer->AddMaterial(Ca, fractionmass=0.005366);
		soil_typicalloam_seltzer->AddMaterial(Fe, fractionmass=0.009582);

/*
d = 1.5*g/cm3;
  moon =  new G4Material("moon",d,6);
  moon -> AddElement(elO,0.45);
  moon -> AddElement(elMg,0.05);
  moon -> AddElement(elAl,0.13);
  moon -> AddElement(elSi,0.21);
  moon -> AddElement(elCa,0.10);
  moon -> AddElement(elFe,0.06); 

  // wood
  G4Material* wood = new G4Material
    (name="wood", density=0.9*g/cm3, ncomponents=3);
  wood->AddElement(H , 4);
  wood->AddElement(O , 1);
  wood->AddElement(C , 2);

  G4Material* Graphite = 
    new G4Material("Graphite", density= 1.7*g/cm3, ncomponents=1);
  Graphite->AddElement(C, fractionmass=1.);

  G4Material* CH = new G4Material("Plastic", density= 1.04*g/cm3, ncomponents=2);
  CH->AddElement(C, natoms=1);
  CH->AddElement(H, natoms=1);

  G4Material* H2O = new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78*eV);

  G4Material* steam = 
    new G4Material("WaterSteam", density= 1.0*mg/cm3, ncomponents=1);
  steam->AddMaterial(H2O, fractionmass=1.);
  steam->GetIonisation()->SetMeanExcitationEnergy(71.6*eV);  
*/
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;

}
void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());

  //PrintDetectorParameters();//already printed!!
  if (sourceType==0 || sourceType==1) {
	G4cout //<< "\n------------------------------------------------------------"
         << "\n---> The source is " << GetSourceType() << " at distance: "
         << fUpperAir_height << " mm in: " << fworld_material->GetName()		 	 
         << "\n------------------------------------------------------------\n";
  }

  if (sourceType==2) {
	G4cout << "\n------------------------------------------------------------"
         << "\n---> The source is " << GetSourceType() << " ; diameter: "
         << fsource_external_diameter << " mm; height: " <<fsource_total_height
		 <<" mm; composition: "<< fsource_material->GetName()
		 <<"\n"
		 <<" source envelope made of: "<<fsource_envelope_material->GetName()<<" ; thickness [mm] = "<<fsource_envelope_thickness
         << "\n------------------------------------------------------------\n";
  }

  if (sourceType==3) {
	G4cout << "\n------------------------------------------------------------"
         << "\n---> The source is " << GetSourceType() << " ; external diameter: "
         << fsource_external_diameter << " mm; total height (including bottom envelope): " <<fsource_total_height
		 <<" mm; composition: "<< fsource_material->GetName()
		 <<"\n"
		 <<" internal diameter: "<<fsource_internal_diameter<< " mm; upperheight (including envelope): "<<fsource_upper_height<<" mm;"
		 <<"\n"
		 <<" source envelope made of: "<<fsource_envelope_material->GetName()<<" ; thickness [mm] = "<<fsource_envelope_thickness
         << "\n------------------------------------------------------------\n";
  }
}

void DetectorConstruction::DefaultConstruction()
{
  sourceType=0;//default!!
	// Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();  
  
  //Detector coupled with frontal beam parallel with z-axis at detector face!
  //Detector monture data. Its shape is an open tube at the top. At bottom side, it's a very thin compact tube (disc shape) and is part of the monture.
  ftotal_diameter=7.3*cm;  
  ftotal_height=7.3*cm;
  fmonture_thickness=0.05*cm;
  fmonture_mat = nist->FindOrBuildMaterial("G4_Al");

  //Active detector data. Active detector is a tube with center placed on the symmetry axis (z-axis) of the detector monture tube. 
  //Also. the active detector symmetry axis coincides with the monture tube symmetry axis.
  factive_diameter=6.3*cm;
  factive_height=6.3*cm; 
  fTargetMaterial  = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");//, fromIsotopes);
  
  fHoleMaterial = nist->FindOrBuildMaterial("G4_Au");
  hole_diameter = 10.0;//5.85 *mm;
  hole_total_height = 40.0;//14.1 *mm;
  isCoaxial =false;
  //Window and gaps. The window is a very thin compact tube (disc shape) placed at detector monture entrance. In most cases, it is made of same material as the monture.
  //Between window and active detector entrance there is gap, usualy named endCapToDetector filled with vacuum.
  fwindow_thickness=0.05*cm;  
  fWindowMaterial  = nist->FindOrBuildMaterial("G4_Al");//, fromIsotopes);
  
  fendCapToDetector=0.5*cm;
  fgap_mat = nist->FindOrBuildMaterial("G4_Galactic");//vacuum

  fworld_material = nist->FindOrBuildMaterial("G4_AIR");//does not really matter here.

  //-------------------
  fUpperAir_height=3.8*cm;//0.0001*cm;//3.8*cm;//radioactive source placed on top side of detector at this distance
  fFrontalBeam_radius=0.1*cm;//source is parallel frontal beam
  fFrontalBeam_angle=0.0;//frontal beam angle to world z-axis

  fsource_material= nist->FindOrBuildMaterial("G4_WATER");
  fsource_envelope_material= nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
  fsource_external_diameter=7.5*cm;
  fsource_total_height=3.4*cm;
  fsource_envelope_thickness=0.1*cm;

  fsource_internal_diameter=8.55*cm;
  fsource_upper_height=3.8*cm;
}

G4VPhysicalVolume* DetectorConstruction::BuildDetectorWithPointOrBeamSource()
{
//we have:
  //active detector tube
  //one empty tube surrounding active detector made of vacuum
  //two cylinders (on top and on bottom of detector) made of vacuum.
  //this is surrounded by monture (also a tube) which is also presented at bottom (another tube).
  //on top we have the window tube.

  //Construct geometry in Geant4 "language".

  //The world is a surrounding cylinder.
  G4double world_diameter = 1.0*(ftotal_diameter);//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  G4double world_Z  = 1.0*(ftotal_height+fUpperAir_height);//@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  solidWorldD=0;logicWorldD=0;physicalWorldD=0;
  //world solid  
  solidWorldD=new G4Tubs("World_solid", 0.0, 0.5*world_diameter, 0.5*world_Z, 0.0, 2.0*pi);
  
  //world logic
  logicWorldD = new G4LogicalVolume(
						solidWorldD,          //its solid
                        fworld_material,           //its material
                       "World_logic");            //its name
  //world physical
  physicalWorldD = new G4PVPlacement(0,                     //no rotation
                  G4ThreeVector(),       //at (0,0,0) CENTER PLACED AT THESE COORDINATES!
                  logicWorldD,            //its logical volume
                      "World_phys",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      fCheckOverlaps);        //overlaps checking

  
  
  //---some initializations:
  G4double deltah=world_Z;
  G4double zoffset=0.5*deltah;//redundant
  G4ThreeVector pos = G4ThreeVector(0, 0, -zoffset);//redundant
  G4double deltah_save=0.0;
  //---end initialisation

  //------------------Zero level. Source to detector layer
  solidUpperAir=0;logicUpperAir=0;physicalUpperAir=0;
  deltah=deltah-fUpperAir_height;
  zoffset=0.5*deltah;
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.
  solidUpperAir=new G4Tubs("SourceTop_solid", 0, 0.5*ftotal_diameter, 0.5*fUpperAir_height, 0.0, 2.0*pi);
  logicUpperAir =new G4LogicalVolume(solidUpperAir, fworld_material, "SourceTop_logic");                       
  physicalUpperAir=new G4PVPlacement(0, pos, logicUpperAir, "SourceTop_phys", logicWorldD, false, 0, fCheckOverlaps);          

  deltah=deltah-fUpperAir_height;//next zofsset take into account the last half layer

  //-----------------First layer. Detector window. The top
  fsolidWindow=0;fLogicWindow=0;fphysicalWindow=0;
  deltah=deltah-fwindow_thickness;
  zoffset=0.5*deltah;
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.
  fsolidWindow=new G4Tubs("Window_solid", 0, 0.5*ftotal_diameter, 0.5*fwindow_thickness, 0.0, 2.0*pi);
  fLogicWindow =new G4LogicalVolume(fsolidWindow, fWindowMaterial, "Window_logic");                       
  fphysicalWindow=new G4PVPlacement(0, pos, fLogicWindow, "Window_phys", logicWorldD, false, 0, fCheckOverlaps);          

  deltah=deltah-fwindow_thickness;//next zofsset take into account the last half layer

  G4VisAttributes* VisAtt1= new G4VisAttributes(G4Colour(0.3,0.8,0.1));//green!
  VisAtt1->SetVisibility(true);
  VisAtt1->SetForceSolid(true);
  //VisAtt1->SetForceWireframe(true);
  fLogicWindow->SetVisAttributes(VisAtt1);

  //===============
  deltah_save=deltah;//@@@@@@@@@@@@@@save it here!!!
  //===================
  //---------------Next layer. The vacuum at the top of active detector 

  fsolidVacuum_top=0;flogicVacuum_top=0;fphysicalVacuum_top=0;
  deltah=deltah-fendCapToDetector;
  zoffset=0.5*deltah;
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.
  G4double rmax=0.5*ftotal_diameter-fmonture_thickness;//until the monture!
  
  fsolidVacuum_top=new G4Tubs("Vacuum_top_solid", 0, rmax, 0.5*fendCapToDetector, 0.0, 2.0*pi);
  flogicVacuum_top =new G4LogicalVolume(fsolidVacuum_top, fgap_mat, "Vacuum_top_logic");          
  fphysicalVacuum_top=new G4PVPlacement(0, pos, flogicVacuum_top, "Vacuum_top_phys", logicWorldD, false, 0, fCheckOverlaps); 
					
  deltah=deltah-fendCapToDetector;//next zofsset take into account the last half layer

  G4VisAttributes* VisAtt2= new G4VisAttributes(G4Colour(0.2,0.3,0.8));//blue!
  VisAtt2->SetVisibility(true);
  //VisAtt2->SetForceSolid(true);
  flogicVacuum_top->SetVisAttributes(VisAtt2);

  //---------------Next layer. The active detector
  G4Tubs* solidDetector=0;logicDetector=0;physicalDetector=0;

  deltah=deltah-factive_height;//deltah=deltah-active_height;
  zoffset=0.5*deltah;
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.  
  solidDetector=new G4Tubs("Detector_solid", 0.0, 0.5*factive_diameter, 0.5*factive_height, 0.0, 2.2*pi);                
  //logicDetector = new G4LogicalVolume(solidDetector, fTargetMaterial, "Detector_logic");  
  //physicalDetector = new G4PVPlacement(0, pos, logicDetector, "Detector_phys", logicWorldD, false, 0, fCheckOverlaps); 
  
  deltah=deltah-factive_height;//next zofsset take into account the last half layer
  //=========================COAXIAL========================
  if(isCoaxial){
  G4Tubs* solidHole = new G4Tubs("hole_solid", 0.0, 0.5*hole_diameter, 0.5*hole_total_height, 0.0, 2.0*pi); 
  G4double zHolePlace=-0.5*factive_height+0.5*hole_total_height;
  G4SubtractionSolid* coaxDet = new G4SubtractionSolid("coaxDet_solid", solidDetector, solidHole
	  ,0,//rot matrix
	  G4ThreeVector(0.0 *cm, 0.0 *cm, -zHolePlace)//YESS, minus here means to back side!!!- with - means + and radiation comes from - =>back side!!
	  );
  
 logicDetector = new G4LogicalVolume(coaxDet, fTargetMaterial, "Detector_logic");  
 physicalDetector = new G4PVPlacement(0, pos, logicDetector, "Detector_phys", logicWorldD, false, 0, fCheckOverlaps); 

 G4ThreeVector pos2 = G4ThreeVector(0, 0, -zoffset-zHolePlace);//zoffset gradually increase so it is increase here!!!
 logicHole = new G4LogicalVolume(solidHole, fHoleMaterial, "Hole_logic"); 
 physicalHole = new G4PVPlacement(0, pos2, logicHole, "Hole_phys", logicWorldD, false, 0, fCheckOverlaps); 

 G4VisAttributes* VisAtt33= new G4VisAttributes(G4Colour(1.0,1.0,0.0));//yellow!
 VisAtt33->SetVisibility(true);
 VisAtt33->SetForceSolid(true);//false);//true);
 logicHole->SetVisAttributes(VisAtt33); 
  } else{
	logicDetector = new G4LogicalVolume(solidDetector, fTargetMaterial, "Detector_logic");  
    physicalDetector = new G4PVPlacement(0, pos, logicDetector, "Detector_phys", logicWorldD, false, 0, fCheckOverlaps);
  }
  //==============================================================
  G4VisAttributes* VisAtt3= new G4VisAttributes(G4Colour(0.8,0.2,0.3));//red!
  VisAtt3->SetVisibility(true);
  VisAtt3->SetForceSolid(true);
  //detector_logic->SetVisAttributes(VisAtt3);
  logicDetector->SetVisAttributes(VisAtt3);

  //----------also we have (at same layer) the vacuum suroundings
  fsolidVacuum_surface=0;flogicVacuum_surface=0;fphysicalVacuum_surface=0;

  G4double rmin=0.5*factive_diameter;//inner radius
  rmax=0.5*ftotal_diameter-fmonture_thickness;//until the monture!  
  fsolidVacuum_surface=new G4Tubs("Vacuum_surface_solid", rmin, rmax, 0.5*factive_height, 0.0, 2.0*pi);
  flogicVacuum_surface= new G4LogicalVolume(fsolidVacuum_surface, fgap_mat, "Vacuum_surface_logic");
  fphysicalVacuum_surface= new G4PVPlacement(0, pos, flogicVacuum_surface, "Vacuum_surface_phys", logicWorldD, false, 0, fCheckOverlaps); 
					
  //deltah is already updated earlier!
  flogicVacuum_surface->SetVisAttributes(VisAtt2);

  //---------------Next layer=remainder vacuum--------------
  fsolidVacuum_bottom=0;flogicVacuum_bottom=0;fphysicalVacuum_bottom=0;

  G4double remainder=ftotal_height-fwindow_thickness-fendCapToDetector-factive_height-fmonture_thickness;
  deltah=deltah-remainder;
  zoffset=0.5*deltah;
  rmax=0.5*ftotal_diameter-fmonture_thickness;//until the monture!
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.  
  fsolidVacuum_bottom=new G4Tubs("Vacuum_bottom_solid", 0.0, rmax, 0.5*remainder, 0.0, 2.0*pi);
  flogicVacuum_bottom = new G4LogicalVolume(fsolidVacuum_bottom, fgap_mat, "Vacuum_bottom_logic"); 
  fphysicalVacuum_bottom=new G4PVPlacement(0, pos, flogicVacuum_bottom, "Vacuum_bottom_phys", logicWorldD, false, 0, fCheckOverlaps);   
					
  deltah=deltah-remainder;//next zofsset take into account the last half layer
  
  flogicVacuum_bottom->SetVisAttributes(VisAtt2);

  //-------------Next layer. Detector monture. The bottom=last layer
  fsolidMonture_bottom=0;flogicMonture_bottom=0;fphysicalMonture_bottom=0;
  deltah=deltah-fmonture_thickness;
  zoffset=0.5*deltah;
  //zoffset=0.5*ftotal_height-0.5*fmonture_thickness;
  pos = G4ThreeVector(0, 0, -zoffset);//pos = G4ThreeVector(0, 0, zoffset);//its position ...at bottom, toward pozitive z-axis
  fsolidMonture_bottom=new G4Tubs("Monture_bottom_solid", 0, 0.5*ftotal_diameter, 0.5*fmonture_thickness, 0.0, 2.0*pi);
  flogicMonture_bottom = new G4LogicalVolume(fsolidMonture_bottom, fmonture_mat, "Monture_bottom_logic");   
  fphysicalMonture_bottom=new G4PVPlacement(0, pos, flogicMonture_bottom, "Monture_bottom_phys", logicWorldD, false, 0, fCheckOverlaps); 

  deltah=deltah-fmonture_thickness;
  //----------------------------------------------------------------------
  //---------------------------------------------------------------------
  //Detector monture. The outer surface.
  fsolidMonture_surface=0;flogicMonture_surface=0;fphysicalMonture_surface=0;

  rmin=0.5*ftotal_diameter-fmonture_thickness;//inner radius
  G4double halfheight=0.5*(ftotal_height-fmonture_thickness-fwindow_thickness);
  deltah_save=deltah_save-2.0*halfheight;
  zoffset=0.5*deltah_save;
  //zoffset=0.5*(fwindow_thickness-fmonture_thickness);//to match with world at top!
  pos = G4ThreeVector(0, 0, -zoffset);//pos = G4ThreeVector(0, 0, zoffset);//its position displaced by true value!
  fsolidMonture_surface=new G4Tubs("Monture_surface_solid", rmin, 0.5*ftotal_diameter, halfheight, 0.0, 2.0*pi);
  flogicMonture_surface =new G4LogicalVolume(fsolidMonture_surface, fmonture_mat, "Monture_surface_logic");
  fphysicalMonture_surface=new G4PVPlacement(0, pos, flogicMonture_surface, "Monture_surface_phys", logicWorldD, false, 0, fCheckOverlaps);    
  //---------------------------------------------
  
  PrintDetectorParameters();
  
  //logicWorldD->SetVisAttributes (G4VisAttributes::Invisible);

  //always return the physical World  
  return physicalWorldD;//world_phys;
}

G4VPhysicalVolume* DetectorConstruction::BuildDetectorWithSarpaganSource()
{
  //The world is a surrounding cylinder.  Top and exterior source envelope is neglected!
  G4double world_diameter = 1.0*(std::max(ftotal_diameter,fsource_external_diameter));//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  G4double world_Z  = 1.0*(ftotal_height+fsource_total_height);//@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  solidWorldD=0;logicWorldD=0;physicalWorldD=0;
  //world solid  
  solidWorldD=new G4Tubs("World_solid", 0.0, 0.5*world_diameter, 0.5*world_Z, 0.0, 2.0*pi);
  
  //world logic
  logicWorldD = new G4LogicalVolume(
						solidWorldD,          //its solid
                        fworld_material,           //its material
                       "World_logic");            //its name
  //world physical
  physicalWorldD = new G4PVPlacement(0,                     //no rotation
                  G4ThreeVector(),       //at (0,0,0) CENTER PLACED AT THESE COORDINATES!
                  logicWorldD,            //its logical volume
                      "World_phys",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      fCheckOverlaps);        //overlaps checking

  
  
  //---some initializations:
  G4double deltah=world_Z;
  G4double zoffset=0.5*deltah;//redundant
  G4ThreeVector pos = G4ThreeVector(0, 0, -zoffset);//redundant
  G4double deltah_save=0.0;
  //---end initialisation

  //------------------Zero level. Source to detector layer
  solidSourceTopCylinder=0;logicSourceTopCylinder=0;physicalSourceTopCylinder=0;
  deltah=deltah-(fsource_total_height-fsource_envelope_thickness);
  zoffset=0.5*deltah;
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.
  solidSourceTopCylinder=new G4Tubs("Sarpagan_solid", 0, 0.5*fsource_external_diameter,
	  0.5*(fsource_total_height-fsource_envelope_thickness), 0.0, 2.0*pi);
  logicSourceTopCylinder =new G4LogicalVolume(solidSourceTopCylinder, fsource_material, "Sarpagan_logic");                       
  physicalSourceTopCylinder=new G4PVPlacement(0, pos, logicSourceTopCylinder, "Sarpagan_phys", logicWorldD, false, 0, fCheckOverlaps);          

  deltah=deltah-(fsource_total_height-fsource_envelope_thickness);//next zofsset take into account the last half layer

  G4VisAttributes* VisAtt0= new G4VisAttributes(G4Colour(0.45,0.25,0.0));//brown
  VisAtt0->SetVisibility(true);
  VisAtt0->SetForceSolid(true);
  logicSourceTopCylinder->SetVisAttributes(VisAtt0);
  //----------------Bottom envelope plastic Top and surface plastic can be neglected.
  solidSourceEnvelope=0;logicSourceEnvelope=0;physicalSourceEnvelope=0;
  deltah=deltah-fsource_envelope_thickness;//fUpperAir_height;
  zoffset=0.5*deltah;
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.
  solidSourceEnvelope=new G4Tubs("Envelope_solid", 0, 0.5*fsource_external_diameter, 0.5*fsource_envelope_thickness, 0.0, 2.0*pi);
  logicSourceEnvelope =new G4LogicalVolume(solidSourceEnvelope, fsource_envelope_material, "Envelope_logic");                       
  physicalSourceEnvelope=new G4PVPlacement(0, pos, logicSourceEnvelope, "Envelope_phys", logicWorldD, false, 0, fCheckOverlaps);          

  deltah=deltah-fsource_envelope_thickness;//next zofsset take into account the last half layer
  //-----------------First layer. Detector window. The top
  fsolidWindow=0;fLogicWindow=0;fphysicalWindow=0;
  deltah=deltah-fwindow_thickness;
  zoffset=0.5*deltah;
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.
  fsolidWindow=new G4Tubs("Window_solid", 0, 0.5*ftotal_diameter, 0.5*fwindow_thickness, 0.0, 2.0*pi);
  fLogicWindow =new G4LogicalVolume(fsolidWindow, fWindowMaterial, "Window_logic");                       
  fphysicalWindow=new G4PVPlacement(0, pos, fLogicWindow, "Window_phys", logicWorldD, false, 0, fCheckOverlaps);          

  deltah=deltah-fwindow_thickness;//next zofsset take into account the last half layer

  G4VisAttributes* VisAtt1= new G4VisAttributes(G4Colour(0.3,0.8,0.1));//green!
  VisAtt1->SetVisibility(true);
  VisAtt1->SetForceSolid(true);
  //VisAtt1->SetForceWireframe(true);
  fLogicWindow->SetVisAttributes(VisAtt1);

  //===============
  deltah_save=deltah;//@@@@@@@@@@@@@@save it here!!!
  //===================
  //---------------Next layer. The vacuum at the top of active detector 

  fsolidVacuum_top=0;flogicVacuum_top=0;fphysicalVacuum_top=0;
  deltah=deltah-fendCapToDetector;
  zoffset=0.5*deltah;
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.
  G4double rmax=0.5*ftotal_diameter-fmonture_thickness;//until the monture!
  
  fsolidVacuum_top=new G4Tubs("Vacuum_top_solid", 0, rmax, 0.5*fendCapToDetector, 0.0, 2.0*pi);
  flogicVacuum_top =new G4LogicalVolume(fsolidVacuum_top, fgap_mat, "Vacuum_top_logic");          
  fphysicalVacuum_top=new G4PVPlacement(0, pos, flogicVacuum_top, "Vacuum_top_phys", logicWorldD, false, 0, fCheckOverlaps); 
					
  deltah=deltah-fendCapToDetector;//next zofsset take into account the last half layer

  G4VisAttributes* VisAtt2= new G4VisAttributes(G4Colour(0.2,0.3,0.8));//blue!
  VisAtt2->SetVisibility(true);
  //VisAtt2->SetForceSolid(true);
  flogicVacuum_top->SetVisAttributes(VisAtt2);

  //---------------Next layer. The active detector
  G4Tubs* solidDetector=0;logicDetector=0;physicalDetector=0;

  deltah=deltah-factive_height;//deltah=deltah-active_height;
  zoffset=0.5*deltah;
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.  

  solidDetector=new G4Tubs("Detector_solid", 0.0, 0.5*factive_diameter, 0.5*factive_height, 0.0, 2.0*pi);                
  //logicDetector = new G4LogicalVolume(solidDetector, fTargetMaterial, "Detector_logic");  
  //physicalDetector = new G4PVPlacement(0, pos, logicDetector, "Detector_phys", logicWorldD, false, 0, fCheckOverlaps); 
  
  deltah=deltah-factive_height;//next zofsset take into account the last half layer
    
  //=========================COAXIAL========================
  if(isCoaxial){
  G4Tubs* solidHole = new G4Tubs("hole_solid", 0.0, 0.5*hole_diameter, 0.5*hole_total_height, 0.0, 2.0*pi); 
  G4double zHolePlace=-0.5*factive_height+0.5*hole_total_height;
  G4SubtractionSolid* coaxDet = new G4SubtractionSolid("coaxDet_solid", solidDetector, solidHole
	  ,0,//rot matrix
	  G4ThreeVector(0.0 *cm, 0.0 *cm, -zHolePlace)//YESS, minus here means to back side!!!- with - means + and radiation comes from - =>back side!!
	  );
  
 logicDetector = new G4LogicalVolume(coaxDet, fTargetMaterial, "Detector_logic");  
 physicalDetector = new G4PVPlacement(0, pos, logicDetector, "Detector_phys", logicWorldD, false, 0, fCheckOverlaps); 

 G4ThreeVector pos2 = G4ThreeVector(0, 0, -zoffset-zHolePlace);//zoffset gradually increase so it is increase here!!!
 logicHole = new G4LogicalVolume(solidHole, fHoleMaterial, "Hole_logic"); 
 physicalHole = new G4PVPlacement(0, pos2, logicHole, "Hole_phys", logicWorldD, false, 0, fCheckOverlaps); 

 G4VisAttributes* VisAtt33= new G4VisAttributes(G4Colour(1.0,1.0,0.0));//yellow!
 VisAtt33->SetVisibility(true);
 VisAtt33->SetForceSolid(true);//false);//true);
 logicHole->SetVisAttributes(VisAtt33); 
  } else{
	logicDetector = new G4LogicalVolume(solidDetector, fTargetMaterial, "Detector_logic");  
    physicalDetector = new G4PVPlacement(0, pos, logicDetector, "Detector_phys", logicWorldD, false, 0, fCheckOverlaps);
  }
  //==============================================================
 G4VisAttributes* VisAtt3= new G4VisAttributes(G4Colour(0.8,0.2,0.3));//red!
 VisAtt3->SetVisibility(true);
 VisAtt3->SetForceSolid(true);//false);//true);
 logicDetector->SetVisAttributes(VisAtt3);

  //----------also we have (at same layer) the vacuum suroundings
  fsolidVacuum_surface=0;flogicVacuum_surface=0;fphysicalVacuum_surface=0;

  G4double rmin=0.5*factive_diameter;//inner radius
  rmax=0.5*ftotal_diameter-fmonture_thickness;//until the monture!  
  fsolidVacuum_surface=new G4Tubs("Vacuum_surface_solid", rmin, rmax, 0.5*factive_height, 0.0, 2.0*pi);
  flogicVacuum_surface= new G4LogicalVolume(fsolidVacuum_surface, fgap_mat, "Vacuum_surface_logic");
  fphysicalVacuum_surface= new G4PVPlacement(0, pos, flogicVacuum_surface, "Vacuum_surface_phys", logicWorldD, false, 0, fCheckOverlaps); 
					
  //deltah is already updated earlier!
  flogicVacuum_surface->SetVisAttributes(VisAtt2);

  //---------------Next layer=remainder vacuum--------------
  fsolidVacuum_bottom=0;flogicVacuum_bottom=0;fphysicalVacuum_bottom=0;

  G4double remainder=ftotal_height-fwindow_thickness-fendCapToDetector-factive_height-fmonture_thickness;
  deltah=deltah-remainder;
  zoffset=0.5*deltah;
  rmax=0.5*ftotal_diameter-fmonture_thickness;//until the monture!
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.  
  fsolidVacuum_bottom=new G4Tubs("Vacuum_bottom_solid", 0.0, rmax, 0.5*remainder, 0.0, 2.0*pi);
  flogicVacuum_bottom = new G4LogicalVolume(fsolidVacuum_bottom, fgap_mat, "Vacuum_bottom_logic"); 
  fphysicalVacuum_bottom=new G4PVPlacement(0, pos, flogicVacuum_bottom, "Vacuum_bottom_phys", logicWorldD, false, 0, fCheckOverlaps);   
					
  deltah=deltah-remainder;//next zofsset take into account the last half layer
  
  flogicVacuum_bottom->SetVisAttributes(VisAtt2);

  //-------------Next layer. Detector monture. The bottom=last layer
  fsolidMonture_bottom=0;flogicMonture_bottom=0;fphysicalMonture_bottom=0;
  deltah=deltah-fmonture_thickness;
  zoffset=0.5*deltah;
  //zoffset=0.5*ftotal_height-0.5*fmonture_thickness;
  pos = G4ThreeVector(0, 0, -zoffset);//pos = G4ThreeVector(0, 0, zoffset);//its position ...at bottom, toward pozitive z-axis
  fsolidMonture_bottom=new G4Tubs("Monture_bottom_solid", 0, 0.5*ftotal_diameter, 0.5*fmonture_thickness, 0.0, 2.0*pi);
  flogicMonture_bottom = new G4LogicalVolume(fsolidMonture_bottom, fmonture_mat, "Monture_bottom_logic");   
  fphysicalMonture_bottom=new G4PVPlacement(0, pos, flogicMonture_bottom, "Monture_bottom_phys", logicWorldD, false, 0, fCheckOverlaps); 

  deltah=deltah-fmonture_thickness;
  //----------------------------------------------------------------------
  //---------------------------------------------------------------------
  //Detector monture. The outer surface.
  fsolidMonture_surface=0;flogicMonture_surface=0;fphysicalMonture_surface=0;

  rmin=0.5*ftotal_diameter-fmonture_thickness;//inner radius
  G4double halfheight=0.5*(ftotal_height-fmonture_thickness-fwindow_thickness);
  deltah_save=deltah_save-2.0*halfheight;
  zoffset=0.5*deltah_save;
  //zoffset=0.5*(fwindow_thickness-fmonture_thickness);//to match with world at top!
  pos = G4ThreeVector(0, 0, -zoffset);//pos = G4ThreeVector(0, 0, zoffset);//its position displaced by true value!
  fsolidMonture_surface=new G4Tubs("Monture_surface_solid", rmin, 0.5*ftotal_diameter, halfheight, 0.0, 2.0*pi);
  flogicMonture_surface =new G4LogicalVolume(fsolidMonture_surface, fmonture_mat, "Monture_surface_logic");
  fphysicalMonture_surface=new G4PVPlacement(0, pos, flogicMonture_surface, "Monture_surface_phys", logicWorldD, false, 0, fCheckOverlaps);    
  //---------------------------------------------
  
  PrintDetectorParameters();
  
  //logicWorldD->SetVisAttributes (G4VisAttributes::Invisible);

  //always return the physical World  
  return physicalWorldD;//world_phys;
}

G4VPhysicalVolume* DetectorConstruction::BuildDetectorWithMarinelliSource()
{
  //The world is a surrounding cylinder.  Top and exterior source envelope is neglected!
  G4double world_diameter = 1.0*(fsource_external_diameter);//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  G4double world_Z  = 1.0*(std::max(ftotal_height+fsource_upper_height,fsource_total_height));//@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  solidWorldD=0;logicWorldD=0;physicalWorldD=0;
  //world solid  
  solidWorldD=new G4Tubs("World_solid", 0.0, 0.5*world_diameter, 0.5*world_Z, 0.0, 2.0*pi);
  
  //world logic
  logicWorldD = new G4LogicalVolume(
						solidWorldD,          //its solid
                        fworld_material,           //its material
                       "World_logic");            //its name
  //world physical
  physicalWorldD = new G4PVPlacement(0,                     //no rotation
                  G4ThreeVector(),       //at (0,0,0) CENTER PLACED AT THESE COORDINATES!
                  logicWorldD,            //its logical volume
                      "World_phys",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      fCheckOverlaps);        //overlaps checking

  
  
  //---some initializations:
  G4double deltah=world_Z;
  G4double zoffset=0.5*deltah;//redundant
  G4ThreeVector pos = G4ThreeVector(0, 0, -zoffset);//redundant
  G4double deltah_save=0.0;
  //---end initialisation

  //------------------Zero level. Source to detector layer
  solidSourceTopCylinder=0;logicSourceTopCylinder=0;physicalSourceTopCylinder=0;
  deltah=deltah-(fsource_upper_height-fsource_envelope_thickness);
  zoffset=0.5*deltah;
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.
  solidSourceTopCylinder=new G4Tubs("Marinelli_top_solid", 0, 0.5*fsource_external_diameter,
	  0.5*(fsource_upper_height-fsource_envelope_thickness), 0.0, 2.0*pi);
  logicSourceTopCylinder =new G4LogicalVolume(solidSourceTopCylinder, fsource_material, "Marinelli_top_logic");                       
  physicalSourceTopCylinder=new G4PVPlacement(0, pos, logicSourceTopCylinder, "Marinelli_top_phys", logicWorldD, false, 0, fCheckOverlaps);          

  deltah=deltah-(fsource_upper_height-fsource_envelope_thickness);//next zofsset take into account the last half layer

  G4VisAttributes* VisAtt0= new G4VisAttributes(G4Colour(0.45,0.25,0.0));//brown
  VisAtt0->SetVisibility(true);
  //VisAtt0->SetForceSolid(true);
  logicSourceTopCylinder->SetVisAttributes(VisAtt0);
  //====================@@@@@@@@@@@@@@@@@MARINELLI OUTER SURFACE!!!
  solidEnvelopeSurfaceCylinder=0;logicEnvelopeSurfaceCylinder=0;physicalEnvelopeSurfaceCylinder=0;
  G4double deltah2=deltah;
    
  G4double remainder2=fsource_total_height-fsource_upper_height-fsource_envelope_thickness;
  deltah2=deltah2-remainder2;
  zoffset=0.5*deltah2;
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.
  G4double rmin2=0.5*fsource_internal_diameter;//inner radius
  G4double rmax2=0.5*fsource_internal_diameter+fsource_envelope_thickness;
  solidEnvelopeSurfaceCylinder=new G4Tubs("Envelope_surface_solid", rmin2, rmax2,
	  0.5*remainder2, 0.0, 2.0*pi);
  logicEnvelopeSurfaceCylinder =new G4LogicalVolume(solidEnvelopeSurfaceCylinder, fsource_envelope_material, "Envelope_surface_logic");                       
  physicalEnvelopeSurfaceCylinder=new G4PVPlacement(0, pos, logicEnvelopeSurfaceCylinder, "Envelope_surface_phys", logicWorldD, false, 0, fCheckOverlaps);          

  solidSourceSurfaceCylinder=0;logicSourceSurfaceCylinder=0;physicalSourceSurfaceCylinder=0;
  rmin2=0.5*fsource_internal_diameter+fsource_envelope_thickness;//inner radius
  rmax2=0.5*fsource_external_diameter;
  solidSourceSurfaceCylinder=new G4Tubs("Source_surface_solid", rmin2, rmax2,
	  0.5*remainder2, 0.0, 2.0*pi);
  logicSourceSurfaceCylinder =new G4LogicalVolume(solidSourceSurfaceCylinder, fsource_material, "Source_surface_logic");                       
  physicalSourceSurfaceCylinder=new G4PVPlacement(0, pos, logicSourceSurfaceCylinder, "Source_surface_phys", logicWorldD, false, 0, fCheckOverlaps);          

  deltah2=deltah2-remainder2;
  logicSourceSurfaceCylinder->SetVisAttributes(VisAtt0);

  solidEnvelopeBottomCylinder=0;logicEnvelopeBottomCylinder=0;physicalEnvelopeBottomCylinder=0;
  deltah2=deltah2-fsource_envelope_thickness;
  zoffset=0.5*deltah2;
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.
  rmin2=0.5*fsource_internal_diameter;//inner radius
  rmax2=0.5*fsource_external_diameter;
  solidEnvelopeBottomCylinder=new G4Tubs("Envelope_bottom_solid", rmin2, rmax2,
	  0.5*fsource_envelope_thickness, 0.0, 2.0*pi);
  logicEnvelopeBottomCylinder =new G4LogicalVolume(solidEnvelopeBottomCylinder, fsource_envelope_material, "Envelope_bottom_logic");                       
  physicalEnvelopeBottomCylinder=new G4PVPlacement(0, pos, logicEnvelopeBottomCylinder, "Envelope_bottom_phys", logicWorldD, false, 0, fCheckOverlaps);          

  //----------------top envelope
  solidSourceEnvelope=0;logicSourceEnvelope=0;physicalSourceEnvelope=0;
  deltah=deltah-fsource_envelope_thickness;
  zoffset=0.5*deltah;
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.
  solidSourceEnvelope=new G4Tubs("Envelope_top_solid", 0, 0.5*ftotal_diameter, 0.5*fsource_envelope_thickness, 0.0, 2.0*pi);
  logicSourceEnvelope =new G4LogicalVolume(solidSourceEnvelope, fsource_envelope_material, "Envelope_top_logic");                       
  physicalSourceEnvelope=new G4PVPlacement(0, pos, logicSourceEnvelope, "Envelope_top_phys", logicWorldD, false, 0, fCheckOverlaps);          

  deltah=deltah-fsource_envelope_thickness;//next zofsset take into account the last half layer
  //==============================================================================
  
  //-----------------First layer. Detector window. The top
  fsolidWindow=0;fLogicWindow=0;fphysicalWindow=0;
  deltah=deltah-fwindow_thickness;
  zoffset=0.5*deltah;
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.
  fsolidWindow=new G4Tubs("Window_solid", 0, 0.5*ftotal_diameter, 0.5*fwindow_thickness, 0.0, 2.0*pi);
  fLogicWindow =new G4LogicalVolume(fsolidWindow, fWindowMaterial, "Window_logic");                       
  fphysicalWindow=new G4PVPlacement(0, pos, fLogicWindow, "Window_phys", logicWorldD, false, 0, fCheckOverlaps);          

  deltah=deltah-fwindow_thickness;//next zofsset take into account the last half layer

  G4VisAttributes* VisAtt1= new G4VisAttributes(G4Colour(0.3,0.8,0.1));//green!
  VisAtt1->SetVisibility(true);
  VisAtt1->SetForceSolid(true);
  //VisAtt1->SetForceWireframe(true);
  fLogicWindow->SetVisAttributes(VisAtt1);

  //===============
  deltah_save=deltah;//@@@@@@@@@@@@@@save it here!!!
  //===================
  //---------------Next layer. The vacuum at the top of active detector 

  fsolidVacuum_top=0;flogicVacuum_top=0;fphysicalVacuum_top=0;
  deltah=deltah-fendCapToDetector;
  zoffset=0.5*deltah;
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.
  G4double rmax=0.5*ftotal_diameter-fmonture_thickness;//until the monture!
  
  fsolidVacuum_top=new G4Tubs("Vacuum_top_solid", 0, rmax, 0.5*fendCapToDetector, 0.0, 2.0*pi);
  flogicVacuum_top =new G4LogicalVolume(fsolidVacuum_top, fgap_mat, "Vacuum_top_logic");          
  fphysicalVacuum_top=new G4PVPlacement(0, pos, flogicVacuum_top, "Vacuum_top_phys", logicWorldD, false, 0, fCheckOverlaps); 
					
  deltah=deltah-fendCapToDetector;//next zofsset take into account the last half layer

  G4VisAttributes* VisAtt2= new G4VisAttributes(G4Colour(0.2,0.3,0.8));//blue!
  VisAtt2->SetVisibility(true);
  //VisAtt2->SetForceSolid(true);
  flogicVacuum_top->SetVisAttributes(VisAtt2);

  //---------------Next layer. The active detector
  G4Tubs* solidDetector=0;logicDetector=0;physicalDetector=0;

  deltah=deltah-factive_height;//deltah=deltah-active_height;
  zoffset=0.5*deltah;
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.  
  solidDetector=new G4Tubs("Detector_solid", 0.0, 0.5*factive_diameter, 0.5*factive_height, 0.0, 2.2*pi);                
  //logicDetector = new G4LogicalVolume(solidDetector, fTargetMaterial, "Detector_logic");  
  //physicalDetector = new G4PVPlacement(0, pos, logicDetector, "Detector_phys", logicWorldD, false, 0, fCheckOverlaps); 
  
  deltah=deltah-factive_height;//next zofsset take into account the last half layer
  //=========================COAXIAL========================
  if(isCoaxial){
  G4Tubs* solidHole = new G4Tubs("hole_solid", 0.0, 0.5*hole_diameter, 0.5*hole_total_height, 0.0, 2.0*pi); 
  G4double zHolePlace=-0.5*factive_height+0.5*hole_total_height;
  G4SubtractionSolid* coaxDet = new G4SubtractionSolid("coaxDet_solid", solidDetector, solidHole
	  ,0,//rot matrix
	  G4ThreeVector(0.0 *cm, 0.0 *cm, -zHolePlace)//YESS, minus here means to back side!!!- with - means + and radiation comes from - =>back side!!
	  );
  
 logicDetector = new G4LogicalVolume(coaxDet, fTargetMaterial, "Detector_logic");  
 physicalDetector = new G4PVPlacement(0, pos, logicDetector, "Detector_phys", logicWorldD, false, 0, fCheckOverlaps); 

 G4ThreeVector pos2 = G4ThreeVector(0, 0, -zoffset-zHolePlace);//zoffset gradually increase so it is increase here!!!
 logicHole = new G4LogicalVolume(solidHole, fHoleMaterial, "Hole_logic"); 
 physicalHole = new G4PVPlacement(0, pos2, logicHole, "Hole_phys", logicWorldD, false, 0, fCheckOverlaps); 

 G4VisAttributes* VisAtt33= new G4VisAttributes(G4Colour(1.0,1.0,0.0));//yellow!
 VisAtt33->SetVisibility(true);
 VisAtt33->SetForceSolid(true);//false);//true);
 logicHole->SetVisAttributes(VisAtt33); 
  } else{
	logicDetector = new G4LogicalVolume(solidDetector, fTargetMaterial, "Detector_logic");  
    physicalDetector = new G4PVPlacement(0, pos, logicDetector, "Detector_phys", logicWorldD, false, 0, fCheckOverlaps);
  }
  //==============================================================
  G4VisAttributes* VisAtt3= new G4VisAttributes(G4Colour(0.8,0.2,0.3));//red!
  VisAtt3->SetVisibility(true);
  VisAtt3->SetForceSolid(true);
  //detector_logic->SetVisAttributes(VisAtt3);
  logicDetector->SetVisAttributes(VisAtt3);

  //----------also we have (at same layer) the vacuum suroundings
  fsolidVacuum_surface=0;flogicVacuum_surface=0;fphysicalVacuum_surface=0;

  G4double rmin=0.5*factive_diameter;//inner radius
  rmax=0.5*ftotal_diameter-fmonture_thickness;//until the monture!  
  fsolidVacuum_surface=new G4Tubs("Vacuum_surface_solid", rmin, rmax, 0.5*factive_height, 0.0, 2.0*pi);
  flogicVacuum_surface= new G4LogicalVolume(fsolidVacuum_surface, fgap_mat, "Vacuum_surface_logic");
  fphysicalVacuum_surface= new G4PVPlacement(0, pos, flogicVacuum_surface, "Vacuum_surface_phys", logicWorldD, false, 0, fCheckOverlaps); 
					
  //deltah is already updated earlier!
  flogicVacuum_surface->SetVisAttributes(VisAtt2);

  //---------------Next layer=remainder vacuum--------------
  fsolidVacuum_bottom=0;flogicVacuum_bottom=0;fphysicalVacuum_bottom=0;

  G4double remainder=ftotal_height-fwindow_thickness-fendCapToDetector-factive_height-fmonture_thickness;
  deltah=deltah-remainder;
  zoffset=0.5*deltah;
  rmax=0.5*ftotal_diameter-fmonture_thickness;//until the monture!
  pos = G4ThreeVector(0, 0, -zoffset);//its position ...at top, toward negative z-axis.  
  fsolidVacuum_bottom=new G4Tubs("Vacuum_bottom_solid", 0.0, rmax, 0.5*remainder, 0.0, 2.0*pi);
  flogicVacuum_bottom = new G4LogicalVolume(fsolidVacuum_bottom, fgap_mat, "Vacuum_bottom_logic"); 
  fphysicalVacuum_bottom=new G4PVPlacement(0, pos, flogicVacuum_bottom, "Vacuum_bottom_phys", logicWorldD, false, 0, fCheckOverlaps);   
					
  deltah=deltah-remainder;//next zofsset take into account the last half layer
  
  flogicVacuum_bottom->SetVisAttributes(VisAtt2);

  //-------------Next layer. Detector monture. The bottom=last layer
  fsolidMonture_bottom=0;flogicMonture_bottom=0;fphysicalMonture_bottom=0;
  deltah=deltah-fmonture_thickness;
  zoffset=0.5*deltah;
  //zoffset=0.5*ftotal_height-0.5*fmonture_thickness;
  pos = G4ThreeVector(0, 0, -zoffset);//pos = G4ThreeVector(0, 0, zoffset);//its position ...at bottom, toward pozitive z-axis
  fsolidMonture_bottom=new G4Tubs("Monture_bottom_solid", 0, 0.5*ftotal_diameter, 0.5*fmonture_thickness, 0.0, 2.0*pi);
  flogicMonture_bottom = new G4LogicalVolume(fsolidMonture_bottom, fmonture_mat, "Monture_bottom_logic");   
  fphysicalMonture_bottom=new G4PVPlacement(0, pos, flogicMonture_bottom, "Monture_bottom_phys", logicWorldD, false, 0, fCheckOverlaps); 

  deltah=deltah-fmonture_thickness;
  //----------------------------------------------------------------------
  //---------------------------------------------------------------------
  //Detector monture. The outer surface.
  fsolidMonture_surface=0;flogicMonture_surface=0;fphysicalMonture_surface=0;

  rmin=0.5*ftotal_diameter-fmonture_thickness;//inner radius
  G4double halfheight=0.5*(ftotal_height-fmonture_thickness-fwindow_thickness);
  deltah_save=deltah_save-2.0*halfheight;
  zoffset=0.5*deltah_save;
  //zoffset=0.5*(fwindow_thickness-fmonture_thickness);//to match with world at top!
  pos = G4ThreeVector(0, 0, -zoffset);//pos = G4ThreeVector(0, 0, zoffset);//its position displaced by true value!
  fsolidMonture_surface=new G4Tubs("Monture_surface_solid", rmin, 0.5*ftotal_diameter, halfheight, 0.0, 2.0*pi);
  flogicMonture_surface =new G4LogicalVolume(fsolidMonture_surface, fmonture_mat, "Monture_surface_logic");
  fphysicalMonture_surface=new G4PVPlacement(0, pos, flogicMonture_surface, "Monture_surface_phys", logicWorldD, false, 0, fCheckOverlaps);    
  //---------------------------------------------
  
  PrintDetectorParameters();
  
  //logicWorldD->SetVisAttributes (G4VisAttributes::Invisible);

  //always return the physical World  
  return physicalWorldD;//world_phys;
}

G4VPhysicalVolume* DetectorConstruction::ConstructDetector()
{
	// Clean old geometry, if any
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  if (sourceType==0 || sourceType==1) return BuildDetectorWithPointOrBeamSource();
  else if (sourceType==2) return BuildDetectorWithSarpaganSource();
  else if (sourceType==3) return BuildDetectorWithMarinelliSource();
  else
	  return BuildDetectorWithPointOrBeamSource();//default  
}

void DetectorConstruction::PrintDetectorParameters()
{
  G4cout << "\n------------------------------------------------------------"
         << "\n---> The detector is " << factive_diameter << " mm active diameter and "
         << factive_height << " mm active height of crystal: " << fTargetMaterial->GetName()
		 << "\n"
         << " window made of "<< fWindowMaterial->GetName()<<" ; thickness [mm] = "<<fwindow_thickness
 		 << "\n"
         << " monture made of "<< fmonture_mat->GetName()<<" ; thickness [mm] = "<<fmonture_thickness
 		 << "\n"
         << " gap made of "<< fgap_mat->GetName()<<" ; endCapToDet [mm] = "<<fendCapToDetector
		 << "\n"
         << " Total diameter [mm] = "<< ftotal_diameter<<" ; total height [mm] = "<<ftotal_height;		 
         //<< "\n------------------------------------------------------------\n";
  if (isCoaxial){
	  G4cout << "\n"
		  <<" Coaxial hole is filled with "<<fHoleMaterial->GetName()<<" ; diameter [mm] "<<hole_diameter<<" ; height [mm] "<<hole_total_height;
  }
  G4cout << "\n------------------------------------------------------------\n";
}

void DetectorConstruction::SetTargetMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();
  G4bool fromIsotopes = false;

  G4Material* pttoMaterial = 
              nistManager->FindOrBuildMaterial(materialName, fromIsotopes);

  if (fTargetMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fTargetMaterial = pttoMaterial;
        if (logicDetector) logicDetector->SetMaterial(fTargetMaterial);
        G4cout << "\n----> The target is made of " << materialName << G4endl;
     } else {
        G4cout << "\n-->  WARNING from SetTargetMaterial : "
               << materialName << " not found" << G4endl;
     }
  }
}

void DetectorConstruction::SetHoleMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();
  G4bool fromIsotopes = false;

  G4Material* pttoMaterial = 
              nistManager->FindOrBuildMaterial(materialName, fromIsotopes);

  if (fHoleMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fHoleMaterial = pttoMaterial;
        if (logicHole) logicHole->SetMaterial(fHoleMaterial);
        G4cout << "\n----> The hole is made of " << materialName << G4endl;
     } else {
        G4cout << "\n-->  WARNING from SetHoleMaterial : "
               << materialName << " not found" << G4endl;
     }
  }
}

void DetectorConstruction::SetWindowMaterial(G4String materialName)
{	
  G4NistManager* nistManager = G4NistManager::Instance();
  G4bool fromIsotopes = false;

  G4Material* pttoMaterial = 
              nistManager->FindOrBuildMaterial(materialName, fromIsotopes);

  if (fWindowMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fWindowMaterial = pttoMaterial;
        if (fLogicWindow) fLogicWindow->SetMaterial(fWindowMaterial);
        G4cout << "\n----> The window is made of " << materialName << G4endl;
     } else {
        G4cout << "\n-->  WARNING from SetWindowMaterial : "
               << materialName << " not found" << G4endl;
     }
  }
}

void DetectorConstruction::SetActive_height(G4double value)
{
  factive_height=value;
}

void DetectorConstruction::SetActive_diameter(G4double value)
{
  factive_diameter=value;
}

void DetectorConstruction::SetHole_height(G4double value)
{
  hole_total_height=value;
}

void DetectorConstruction::SetHole_diameter(G4double value)
{
  hole_diameter=value;
}
//------------
void DetectorConstruction::SetWorldMaterial (G4String materialName)
{
	G4NistManager* nistManager = G4NistManager::Instance();
  G4bool fromIsotopes = false;

  G4Material* pttoMaterial = 
              nistManager->FindOrBuildMaterial(materialName, fromIsotopes);

  if (fworld_material != pttoMaterial) {
     if ( pttoMaterial ) {
        fworld_material = pttoMaterial;
        if (logicWorldD) logicWorldD->SetMaterial(fworld_material);
        G4cout << "\n----> The world is made of " << materialName << G4endl;
     } else {
        G4cout << "\n-->  WARNING from SetWorldMaterial : "
               << materialName << " not found" << G4endl;
     }
  }
}
void DetectorConstruction::SetGapMaterial (G4String materialName)
{

	//G4Material* pttoMaterialqq = G4Material::GetMaterial(materialName);
  //if (pttoMaterialqq) GapMaterial = pttoMaterialqq;
  //================
		G4NistManager* nistManager = G4NistManager::Instance();
  G4bool fromIsotopes = false;

  G4Material* pttoMaterial = 
              nistManager->FindOrBuildMaterial(materialName, fromIsotopes);

  if (fgap_mat != pttoMaterial) {
     if ( pttoMaterial ) {
        fgap_mat = pttoMaterial;
        if (flogicVacuum_bottom) flogicVacuum_bottom->SetMaterial(fgap_mat);
		if (flogicVacuum_surface) flogicVacuum_surface->SetMaterial(fgap_mat);
		if (flogicVacuum_top) flogicVacuum_top->SetMaterial(fgap_mat);
        G4cout << "\n----> The gap is made of " << materialName << G4endl;
     } else {
        G4cout << "\n-->  WARNING from SetGapMaterial : "
               << materialName << " not found" << G4endl;
     }
  }
}
void DetectorConstruction::SetMontureMaterial (G4String materialName)
{
		G4NistManager* nistManager = G4NistManager::Instance();
  G4bool fromIsotopes = false;

  G4Material* pttoMaterial = 
              nistManager->FindOrBuildMaterial(materialName, fromIsotopes);

  if (fmonture_mat != pttoMaterial) {
     if ( pttoMaterial ) {
        fmonture_mat = pttoMaterial;
        if (flogicMonture_bottom) flogicMonture_bottom->SetMaterial(fmonture_mat);
		if (flogicMonture_surface) flogicMonture_surface->SetMaterial(fmonture_mat);
        G4cout << "\n----> The monture is made of " << materialName << G4endl;
     } else {
        G4cout << "\n-->  WARNING from SetMontureMaterial : "
               << materialName << " not found" << G4endl;
     }
  }
}
void DetectorConstruction::SetTotal_diameter (G4double value)
{
	ftotal_diameter=value;
}
void DetectorConstruction::SetTotal_height (G4double value)
{
	ftotal_height=value;
}
void DetectorConstruction::SetEndCapToDetector (G4double value)
{
	fendCapToDetector=value;
}
void DetectorConstruction::SetWindow_thickness (G4double value)
{
	fwindow_thickness=value;
}
void DetectorConstruction::SetMonture_thickness (G4double value)
{
	fmonture_thickness=value;
}

void DetectorConstruction::SetPointSourceOrFrontalBeamToDetectorDistance (G4double value)
{
	fUpperAir_height=value;
}

void DetectorConstruction::SetFrontalBeamRadius (G4double value)
{
	fFrontalBeam_radius=value;
}

void DetectorConstruction::SetFrontalBeamAngle (G4double value)
{
	fFrontalBeam_angle=value;
}

void DetectorConstruction::SetSourceMaterial (G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();
  G4bool fromIsotopes = false;

  G4Material* pttoMaterial = 
              nistManager->FindOrBuildMaterial(materialName, fromIsotopes);

  if (fsource_material != pttoMaterial) {
     if ( pttoMaterial ) {
        fsource_material = pttoMaterial;
        if (logicSourceTopCylinder) logicSourceTopCylinder->SetMaterial(fsource_material);
        G4cout << "\n----> The source is made of " << materialName << G4endl;
     } else {
        G4cout << "\n-->  WARNING from SetSourceMaterial : "
               << materialName << " not found" << G4endl;
     }
  }
}

void DetectorConstruction::SetSourceEnvelopeMaterial (G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();
  G4bool fromIsotopes = false;

  G4Material* pttoMaterial = 
              nistManager->FindOrBuildMaterial(materialName, fromIsotopes);

  if (fsource_envelope_material != pttoMaterial) {
     if ( pttoMaterial ) {
        fsource_envelope_material = pttoMaterial;
        if (logicSourceEnvelope) logicSourceEnvelope->SetMaterial(fsource_envelope_material);
        G4cout << "\n----> The source envelope is made of " << materialName << G4endl;
     } else {
        G4cout << "\n-->  WARNING from SetSourceEnvelopeMaterial : "
               << materialName << " not found" << G4endl;
     }
  }
}

void DetectorConstruction::SetSourceExternalDiameter (G4double value)
{
	fsource_external_diameter=value;
}

void DetectorConstruction::SetSourceTotalHeight (G4double value)
{
	fsource_total_height=value;
}

void DetectorConstruction::SetSourceEnvelopeThickness (G4double value)
{
	fsource_envelope_thickness=value;
}

void DetectorConstruction::SetSourceInternalDiameter (G4double value)
{
	fsource_internal_diameter=value;
}

void DetectorConstruction::SetSourceUpperHeight (G4double value)
{
	fsource_upper_height=value;
}
