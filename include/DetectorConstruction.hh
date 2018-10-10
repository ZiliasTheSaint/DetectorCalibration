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
// $Id: DetectorConstruction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Dan Fulea, JUL 2013, National Institute of Public Health Bucharest, Cluj-Napoca Regional Centre, Romania

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Tubs.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class DetectorMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
     
     G4VPhysicalVolume* Construct();
     void UpdateGeometry();
	 //--------------------------------------------
	 void SetSourceType(G4int val){sourceType=val;};
	 //------------------------------------------------
	 void SetTargetMaterial (G4String);
     void SetWindowMaterial(G4String);
     void SetActive_diameter (G4double);
	 void SetActive_height (G4double);
	 void SetWorldMaterial (G4String);
	 void SetGapMaterial (G4String);
	 void SetMontureMaterial (G4String);
	 void SetTotal_diameter (G4double);
	 void SetTotal_height (G4double);
	 void SetEndCapToDetector (G4double);
	 void SetWindow_thickness (G4double);
	 void SetMonture_thickness (G4double);
	 //-----------------------------------------------
	 void SetPointSourceOrFrontalBeamToDetectorDistance(G4double);
	 void SetFrontalBeamRadius(G4double);
	 void SetFrontalBeamAngle(G4double);
	 //-----------------------------
	 void SetSourceMaterial(G4String);
	 void SetSourceEnvelopeMaterial(G4String);
	 void SetSourceTotalHeight (G4double);
	 void SetSourceExternalDiameter (G4double);
	 void SetSourceEnvelopeThickness (G4double);

	 void SetSourceInternalDiameter (G4double);
	 void SetSourceUpperHeight (G4double);

	 void SetHoleMaterial(G4String);
	 void SetHole_diameter (G4double);
	 void SetHole_height (G4double);
	 void SetIfCoaxial(G4String val){
		 if (val=="on")
			 isCoaxial=true;
		 else
			 isCoaxial=false;
	 };
 //public://-----------------------
	 G4int GetSourceTypeCode(){return sourceType;};
	 G4String GetSourceType(){
	  G4String s="";
	  if(sourceType==0) s="frontal beam";
	  else if (sourceType==1) s="point source";
	  else if (sourceType==2) s="Sarpagan";
	  else if (sourceType==3) s="Marinelli baker";
	  else s="unidentified; It was set frontal beam as default ";
	  return s;
     };
	 void PrintDetectorParameters(); 
	 G4double GetPointSourceOrFrontalBeamToDetectorDistance()           {return fUpperAir_height;};
	 G4double GetWorldSizeZ()           {return ftotal_height+fUpperAir_height;};
	 G4double GetWorldSizeRadius()           {return 0.5*ftotal_diameter;};
	 G4double GetFrontalBeamRadius()           {return fFrontalBeam_radius;};
	 G4double GetFrontalBeamAngle()           {return fFrontalBeam_angle;};

	 G4double GetDetectorHeight(){return ftotal_height;};
	 G4double GetDetectorDiameter(){return ftotal_diameter;};
	 G4double GetSourceExternalDiameter(){return fsource_external_diameter;};
	 G4double GetSourceTotalHeight(){return fsource_total_height;};
	 G4double GetSourceEnvelopeThickness(){return fsource_envelope_thickness;};

	 G4double GetSourceInternalDiameter(){return fsource_internal_diameter;};
	 G4double GetSourceUpperHeight(){return fsource_upper_height;};
	 ///The scoring volume. If necessary add another:
	 const G4VPhysicalVolume* GetDetector()  {return physicalDetector;};//{return physiAbsorber;};//@@@@@@@@@
	 G4LogicalVolume* GetDetector_logical(){return logicDetector;};//used for computing mass and dose
            
  private:
     G4int sourceType;

	 G4Tubs*            solidWorldD;    //pointer to the solid World 
     G4LogicalVolume*   logicWorldD;    //pointer to the logical World
     G4VPhysicalVolume* physicalWorldD;    //pointer to the physical World

	 //G4Tubs*            solidDetector; //pointer to the solid active detector
     G4LogicalVolume*   logicDetector; //pointer to the logical active detector
     G4VPhysicalVolume* physicalDetector; //pointer to the physical active detector
	 G4LogicalVolume* logicHole;
	 G4VPhysicalVolume* physicalHole;
    //G4LogicalVolume*   fLogicTarget;     // pointer to the logical Target, i.e. the active detector
    //G4LogicalVolume**  fLogicChamber;    // pointer to the logical Chamber. Note ** only if fLogicchamber is an array!
											 //!!!!..see example B2
	//no needs for envelope, everything is contained in world mother
	 G4Tubs*            solidUpperAir;//source-detector media    
     G4LogicalVolume*   logicUpperAir;    
     G4VPhysicalVolume* physicalUpperAir;    

	 G4Tubs*            fsolidWindow;
	 G4LogicalVolume*   fLogicWindow;
	 G4VPhysicalVolume* fphysicalWindow;

	 G4Tubs*            fsolidVacuum_top;
	 G4LogicalVolume*   flogicVacuum_top;
	 G4VPhysicalVolume* fphysicalVacuum_top;

	 G4Tubs*            fsolidVacuum_surface;
	 G4LogicalVolume*   flogicVacuum_surface;
	 G4VPhysicalVolume* fphysicalVacuum_surface;

	 G4Tubs*            fsolidVacuum_bottom;
	 G4LogicalVolume*   flogicVacuum_bottom;
	 G4VPhysicalVolume* fphysicalVacuum_bottom;

	 G4Tubs*            fsolidMonture_bottom;
	 G4LogicalVolume*   flogicMonture_bottom;
	 G4VPhysicalVolume* fphysicalMonture_bottom;

	 G4Tubs*            fsolidMonture_surface;
	 G4LogicalVolume*   flogicMonture_surface;
	 G4VPhysicalVolume* fphysicalMonture_surface;

     G4Material*        fTargetMaterial;  // pointer to the target  material
     G4Material*        fWindowMaterial; // pointer to the window material
	 G4Material*        fHoleMaterial;//pointer to the hole material...contact made usualy of Au (coaxial semiconductor such as HpGe)

	 G4double hole_diameter;
	 G4double hole_total_height;
	 G4bool isCoaxial;
	
	 G4Material* fmonture_mat;
	 G4Material* fgap_mat;
	 G4Material* fworld_material;
	 
	 G4Tubs*            solidSourceTopCylinder;
     G4LogicalVolume*   logicSourceTopCylinder;    
     G4VPhysicalVolume* physicalSourceTopCylinder;   
	 G4Tubs*            solidSourceEnvelope;
     G4LogicalVolume*   logicSourceEnvelope;    
     G4VPhysicalVolume* physicalSourceEnvelope;   
	 G4Material* fsource_material;
	 G4Material* fsource_envelope_material;
	 G4double fsource_external_diameter;
	 G4double fsource_total_height;
	 G4double fsource_envelope_thickness;
	 G4double fsource_internal_diameter;
	 G4double fsource_upper_height;

	 G4Tubs*            solidEnvelopeBottomCylinder;
     G4LogicalVolume*   logicEnvelopeBottomCylinder;    
     G4VPhysicalVolume* physicalEnvelopeBottomCylinder; 
	 G4Tubs*            solidEnvelopeSurfaceCylinder;
     G4LogicalVolume*   logicEnvelopeSurfaceCylinder;    
     G4VPhysicalVolume* physicalEnvelopeSurfaceCylinder; 
	 G4Tubs*            solidSourceSurfaceCylinder;
     G4LogicalVolume*   logicSourceSurfaceCylinder;    
     G4VPhysicalVolume* physicalSourceSurfaceCylinder; 

	 G4double factive_diameter;
	 G4double factive_height;
	 G4double ftotal_diameter;
	 G4double ftotal_height;
	 G4double fmonture_thickness;
	 G4double fwindow_thickness;
	 G4double fendCapToDetector;

	 G4double fUpperAir_height;
	 G4double fFrontalBeam_radius;
	 G4double fFrontalBeam_angle;

	 G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps 
	 
     DetectorMessenger* detectorMessenger;  //pointer to the Messenger
      
  private:
    
	 void DefineMaterials();
	 void DefaultConstruction();
	 G4VPhysicalVolume* BuildDetectorWithPointOrBeamSource();
	 G4VPhysicalVolume* BuildDetectorWithSarpaganSource();
	 G4VPhysicalVolume* BuildDetectorWithMarinelliSource();
	 G4VPhysicalVolume* ConstructDetector();
};



//inline void DetectorConstruction::ComputeCalorParameters()
//{
  // Compute derived parameters of the calorimeter
     /*LayerThickness = AbsorberThickness + GapThickness;
     CalorThickness = NbOfLayers*LayerThickness;
     
     WorldSizeX = 1.2*CalorThickness; WorldSizeYZ = 1.2*CalorSizeYZ;*/
//}



#endif

