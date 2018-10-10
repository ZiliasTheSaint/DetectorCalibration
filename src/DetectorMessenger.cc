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
// $Id: DetectorMessenger.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Dan Fulea, JUL 2013, National Institute of Public Health Bucharest, Cluj-Napoca Regional Centre, Romania

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

DetectorMessenger::DetectorMessenger(
                                           DetectorConstruction* Det)
:Detector(Det)
{ 
  fB2Directory = new G4UIdirectory("/gdet/");
  fB2Directory->SetGuidance("UI commands specific to this example.");

  fDetDirectory = new G4UIdirectory("/gdet/det/");
  fDetDirectory->SetGuidance("Detector construction control");

  fTargMatCmd = new G4UIcmdWithAString("/gdet/det/setTargetMaterial",this);
  fTargMatCmd->SetGuidance("Select Material of the active detector.");
  fTargMatCmd->SetParameterName("choice",false);
  fTargMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fWindowMatCmd = new G4UIcmdWithAString("/gdet/det/setWindowMaterial",this);
  fWindowMatCmd->SetGuidance("Select Material of the Window.");
  fWindowMatCmd->SetParameterName("choice",false);
  fWindowMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  factive_diameterCmd = new G4UIcmdWithADoubleAndUnit("/gdet/det/active_diameter",this);
  factive_diameterCmd->SetGuidance("Define active diameter.");
  factive_diameterCmd->SetParameterName("Size",false);
  factive_diameterCmd->SetRange("Size>0.");
  factive_diameterCmd->SetUnitCategory("Length");
  factive_diameterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  //factive_diameterCmd->AvailableForStates(G4State_Idle);

  factive_heightCmd = new G4UIcmdWithADoubleAndUnit("/gdet/det/active_height",this);
  factive_heightCmd->SetGuidance("Define active height.");
  factive_heightCmd->SetParameterName("Size",false);
  factive_heightCmd->SetRange("Size>0.");
  factive_heightCmd->SetUnitCategory("Length");
  factive_heightCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  //factive_heightCmd->AvailableForStates(G4State_Idle);
  //-------------------------  
  fworld_matCmd = new G4UIcmdWithAString("/gdet/det/setWorldMaterial",this);
  fworld_matCmd->SetGuidance("Select Material of the world.");
  fworld_matCmd->SetParameterName("choice",false);
  fworld_matCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
  fgap_matCmd = new G4UIcmdWithAString("/gdet/det/setGapMaterial",this);
  fgap_matCmd->SetGuidance("Select Material of the gap.");
  fgap_matCmd->SetParameterName("choice",false);
  fgap_matCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
  fmonture_matCmd = new G4UIcmdWithAString("/gdet/det/setMontureMaterial",this);
  fmonture_matCmd->SetGuidance("Select Material of the monture.");
  fmonture_matCmd->SetParameterName("choice",false);
  fmonture_matCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ftotal_diameterCmd = new G4UIcmdWithADoubleAndUnit("/gdet/det/total_diameter",this);
  ftotal_diameterCmd->SetGuidance("Define total diameter.");
  ftotal_diameterCmd->SetParameterName("Size",false);
  ftotal_diameterCmd->SetRange("Size>0.");
  ftotal_diameterCmd->SetUnitCategory("Length");
  ftotal_diameterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
  ftotal_heightCmd = new G4UIcmdWithADoubleAndUnit("/gdet/det/total_height",this);
  ftotal_heightCmd->SetGuidance("Define total height.");
  ftotal_heightCmd->SetParameterName("Size",false);
  ftotal_heightCmd->SetRange("Size>0.");
  ftotal_heightCmd->SetUnitCategory("Length");
  ftotal_heightCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  	
  fendCapToDetectorCmd = new G4UIcmdWithADoubleAndUnit("/gdet/det/endCapToDetector",this);
  fendCapToDetectorCmd->SetGuidance("Define end cap to detector distance (upper gap).");
  fendCapToDetectorCmd->SetParameterName("Size",false);
  fendCapToDetectorCmd->SetRange("Size>0.");
  fendCapToDetectorCmd->SetUnitCategory("Length");
  fendCapToDetectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  	
  fwindow_thicknessCmd = new G4UIcmdWithADoubleAndUnit("/gdet/det/window_thickness",this);
  fwindow_thicknessCmd->SetGuidance("Define window thickness.");
  fwindow_thicknessCmd->SetParameterName("Size",false);
  fwindow_thicknessCmd->SetRange("Size>0.");
  fwindow_thicknessCmd->SetUnitCategory("Length");
  fwindow_thicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
  fmonture_thicknessCmd = new G4UIcmdWithADoubleAndUnit("/gdet/det/monture_thickness",this);
  fmonture_thicknessCmd->SetGuidance("Define monture thickness.");
  fmonture_thicknessCmd->SetParameterName("Size",false);
  fmonture_thicknessCmd->SetRange("Size>0.");
  fmonture_thicknessCmd->SetUnitCategory("Length");
  fmonture_thicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  //-------------------------
  fsourceTopCmd = new G4UIcmdWithADoubleAndUnit("/gdet/det/pointSourceOrFrontalBeamToDetectorDistance",this);
  fsourceTopCmd->SetGuidance("Define point source (or frontal beam) to detector distance.");
  fsourceTopCmd->SetParameterName("Size",false);
  fsourceTopCmd->SetRange("Size>0.");
  fsourceTopCmd->SetUnitCategory("Length");
  fsourceTopCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ffrontalBeamRadiusCmd = new G4UIcmdWithADoubleAndUnit("/gdet/det/frontalBeamRadius",this);
  ffrontalBeamRadiusCmd->SetGuidance("Define frontal beam radius.");
  ffrontalBeamRadiusCmd->SetParameterName("Size",false);
  ffrontalBeamRadiusCmd->SetRange("Size>0.");
  ffrontalBeamRadiusCmd->SetUnitCategory("Length");
  ffrontalBeamRadiusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ffrontalBeamAngleCmd = new G4UIcmdWithADoubleAndUnit("/gdet/det/frontalBeamAngle",this);
  ffrontalBeamAngleCmd->SetGuidance("Define frontal beam angle.");
  ffrontalBeamAngleCmd->SetParameterName("Angle",false);
  ffrontalBeamAngleCmd->SetUnitCategory("Angle");
  ffrontalBeamAngleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  //-------------
  sourceTypeCmd = new G4UIcmdWithAnInteger("/gdet/gun/source",this);
  sourceTypeCmd->SetGuidance("Source types.");
  sourceTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  //---------------

  fsource_matCmd = new G4UIcmdWithAString("/gdet/det/setSourceMaterial",this);
  fsource_matCmd->SetGuidance("Select source composition.");
  fsource_matCmd->SetParameterName("choice",false);
  fsource_matCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fsource_envelope_matCmd = new G4UIcmdWithAString("/gdet/det/setSourceEnvelopeMaterial",this);
  fsource_envelope_matCmd->SetGuidance("Select source envelope composition.");
  fsource_envelope_matCmd->SetParameterName("choice",false);
  fsource_envelope_matCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fsource_envelope_thicknessCmd = new G4UIcmdWithADoubleAndUnit("/gdet/det/source_envelope_thickness",this);
  fsource_envelope_thicknessCmd->SetGuidance("Define source envelope thickness.");
  fsource_envelope_thicknessCmd->SetParameterName("Size",false);
  fsource_envelope_thicknessCmd->SetRange("Size>0.");
  fsource_envelope_thicknessCmd->SetUnitCategory("Length");
  fsource_envelope_thicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fsource_total_heightCmd = new G4UIcmdWithADoubleAndUnit("/gdet/det/source_total_height",this);
  fsource_total_heightCmd->SetGuidance("Define source total height.");
  fsource_total_heightCmd->SetParameterName("Size",false);
  fsource_total_heightCmd->SetRange("Size>0.");
  fsource_total_heightCmd->SetUnitCategory("Length");
  fsource_total_heightCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fsource_external_diameterCmd = new G4UIcmdWithADoubleAndUnit("/gdet/det/source_external_diameter",this);
  fsource_external_diameterCmd->SetGuidance("Define source external diameter.");
  fsource_external_diameterCmd->SetParameterName("Size",false);
  fsource_external_diameterCmd->SetRange("Size>0.");
  fsource_external_diameterCmd->SetUnitCategory("Length");
  fsource_external_diameterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fsource_internal_diameterCmd = new G4UIcmdWithADoubleAndUnit("/gdet/det/source_internal_diameter",this);
  fsource_internal_diameterCmd->SetGuidance("Define source internal diameter. Marinelli baker only");
  fsource_internal_diameterCmd->SetParameterName("Size",false);
  fsource_internal_diameterCmd->SetRange("Size>0.");
  fsource_internal_diameterCmd->SetUnitCategory("Length");
  fsource_internal_diameterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fsource_upper_heightCmd = new G4UIcmdWithADoubleAndUnit("/gdet/det/source_upper_height",this);
  fsource_upper_heightCmd->SetGuidance("Define source upper height. Marinelli baker only");
  fsource_upper_heightCmd->SetParameterName("Size",false);
  fsource_upper_heightCmd->SetRange("Size>0.");
  fsource_upper_heightCmd->SetUnitCategory("Length");
  fsource_upper_heightCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  isCoaxialCmd = new G4UIcmdWithAString("/gdet/det/isCoaxial",this);
  isCoaxialCmd->SetGuidance("Set on if semiconductor detector such as Ge, off otherwise");
  isCoaxialCmd->SetGuidance("  Choice : on, off(default)");
  isCoaxialCmd->SetParameterName("choice",true);
  isCoaxialCmd->SetDefaultValue("off");
  isCoaxialCmd->SetCandidates("on off");
  isCoaxialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  holeMatCmd = new G4UIcmdWithAString("/gdet/det/setHoleMaterial",this);
  holeMatCmd->SetGuidance("Select Material of the coaxial detector hole.");
  holeMatCmd->SetParameterName("choice",false);
  holeMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  hole_diameterCmd = new G4UIcmdWithADoubleAndUnit("/gdet/det/hole_diameter",this);
  hole_diameterCmd->SetGuidance("Define hole diameter.");
  hole_diameterCmd->SetParameterName("Size",false);
  hole_diameterCmd->SetRange("Size>0.");
  hole_diameterCmd->SetUnitCategory("Length");
  hole_diameterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  hole_heightCmd = new G4UIcmdWithADoubleAndUnit("/gdet/det/hole_height",this);
  hole_heightCmd->SetGuidance("Define hole height.");
  hole_heightCmd->SetParameterName("Size",false);
  hole_heightCmd->SetRange("Size>0.");
  hole_heightCmd->SetUnitCategory("Length");
  hole_heightCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/gdet/det/update",this);
  UpdateCmd->SetGuidance("Update detector geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
  
}



DetectorMessenger::~DetectorMessenger()
{  
	delete fTargMatCmd;
	delete fWindowMatCmd;
	delete factive_diameterCmd;
	delete factive_heightCmd;
  
	delete fB2Directory;
	delete fDetDirectory;

	delete fworld_matCmd;
	delete fgap_matCmd;
	delete fmonture_matCmd;
    delete ftotal_diameterCmd;
    delete ftotal_heightCmd;
    delete fendCapToDetectorCmd;
    delete fwindow_thicknessCmd;
	delete fmonture_thicknessCmd;

	delete fsourceTopCmd;
	delete ffrontalBeamRadiusCmd;
	delete ffrontalBeamAngleCmd;

	delete fsource_matCmd;
	delete fsource_envelope_matCmd;
	delete fsource_total_heightCmd;
	delete fsource_external_diameterCmd;
	delete fsource_envelope_thicknessCmd;

	delete fsource_internal_diameterCmd;
	delete fsource_upper_heightCmd;

	delete hole_diameterCmd;
    delete hole_heightCmd;
	delete holeMatCmd;
	delete isCoaxialCmd;

	delete UpdateCmd;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{   
  if( command == fTargMatCmd )
   { Detector->SetTargetMaterial(newValue);}

  if( command == fWindowMatCmd )
   { Detector->SetWindowMaterial(newValue);}

  if( command == fsource_matCmd )
   { Detector->SetSourceMaterial(newValue);}

  if( command == fsource_envelope_matCmd )
   { Detector->SetSourceEnvelopeMaterial(newValue);}

  if( command == factive_diameterCmd ) {
    Detector
      ->SetActive_diameter(factive_diameterCmd->GetNewDoubleValue(newValue));
  }

  if( command == factive_heightCmd ) {
    Detector
      ->SetActive_height(factive_heightCmd->GetNewDoubleValue(newValue));
  }  
  //-------------------------------------
  
  if( command == fworld_matCmd )
   { Detector->SetWorldMaterial(newValue);}
	
  if( command == fgap_matCmd )
   { Detector->SetGapMaterial(newValue);
  }
	
  if( command == fmonture_matCmd )
   { Detector->SetMontureMaterial(newValue);}
    
  if( command == ftotal_diameterCmd ) {
    Detector
      ->SetTotal_diameter(ftotal_diameterCmd->GetNewDoubleValue(newValue));
  }
    
  if( command == ftotal_heightCmd ) {
    Detector
      ->SetTotal_height(ftotal_heightCmd->GetNewDoubleValue(newValue));
  }
    
  if( command == fendCapToDetectorCmd ) {
    Detector
      ->SetEndCapToDetector(fendCapToDetectorCmd->GetNewDoubleValue(newValue));
  }
    
  if( command == fwindow_thicknessCmd ) {
    Detector
      ->SetWindow_thickness(fwindow_thicknessCmd->GetNewDoubleValue(newValue));
  }
	
  if( command == fmonture_thicknessCmd ) {
    Detector
      ->SetMonture_thickness(fmonture_thicknessCmd->GetNewDoubleValue(newValue));
  }
  //---------
  if( command == fsourceTopCmd ) {
    Detector
      ->SetPointSourceOrFrontalBeamToDetectorDistance(fsourceTopCmd->GetNewDoubleValue(newValue));
  }

  if( command == ffrontalBeamRadiusCmd ) {
    Detector
      ->SetFrontalBeamRadius(ffrontalBeamRadiusCmd->GetNewDoubleValue(newValue));
  }

  if( command == ffrontalBeamAngleCmd ) {
    Detector
      ->SetFrontalBeamAngle(ffrontalBeamAngleCmd->GetNewDoubleValue(newValue));
  }
  //-----------------
  if( command == fsource_total_heightCmd ) {
    Detector
      ->SetSourceTotalHeight(fsource_total_heightCmd->GetNewDoubleValue(newValue));
  }
  if( command == fsource_external_diameterCmd ) {
    Detector
      ->SetSourceExternalDiameter(fsource_external_diameterCmd->GetNewDoubleValue(newValue));
  }
  if( command == fsource_envelope_thicknessCmd ) {
    Detector
      ->SetSourceEnvelopeThickness(fsource_envelope_thicknessCmd->GetNewDoubleValue(newValue));
  }
    
  if( command == sourceTypeCmd )
   { Detector->SetSourceType(sourceTypeCmd->GetNewIntValue(newValue));}

  if( command == fsource_internal_diameterCmd ) {
    Detector
      ->SetSourceInternalDiameter(fsource_internal_diameterCmd->GetNewDoubleValue(newValue));
  }

  if( command == fsource_upper_heightCmd ) {
    Detector
      ->SetSourceUpperHeight(fsource_upper_heightCmd->GetNewDoubleValue(newValue));
  }

  if( command == fsource_envelope_thicknessCmd ) {
    Detector
      ->SetSourceEnvelopeThickness(fsource_envelope_thicknessCmd->GetNewDoubleValue(newValue));
  }
  //=========
  if( command == isCoaxialCmd )
   { Detector->SetIfCoaxial(newValue);
   }
  if( command == holeMatCmd )
   { Detector->SetHoleMaterial(newValue);}
   
  if( command == hole_diameterCmd ) {
    Detector
      ->SetHole_diameter(hole_diameterCmd->GetNewDoubleValue(newValue));
  }

  if( command == hole_heightCmd ) {
    Detector
      ->SetHole_height(hole_heightCmd->GetNewDoubleValue(newValue));
  } 
  //--------
  if( command == UpdateCmd )
   { Detector->UpdateGeometry(); }
  //--------
}

