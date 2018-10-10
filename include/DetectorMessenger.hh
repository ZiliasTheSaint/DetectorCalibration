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
// $Id: DetectorMessenger.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Dan Fulea, JUL 2013, National Institute of Public Health Bucharest, Cluj-Napoca Regional Centre, Romania

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

class DetectorMessenger: public G4UImessenger
{
  public:
    DetectorMessenger(DetectorConstruction* );
   ~DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    DetectorConstruction* Detector;
    
	G4UIdirectory*           fB2Directory;
    G4UIdirectory*           fDetDirectory;

    G4UIcmdWithAString*      fTargMatCmd;//active detector
    G4UIcmdWithAString*      fWindowMatCmd;//window material

    G4UIcmdWithADoubleAndUnit* factive_diameterCmd;
    G4UIcmdWithADoubleAndUnit* factive_heightCmd;
    G4UIcmdWithADoubleAndUnit* hole_diameterCmd;
    G4UIcmdWithADoubleAndUnit* hole_heightCmd;
	G4UIcmdWithAString*      holeMatCmd;
	G4UIcmdWithAString*      isCoaxialCmd;

	G4UIcmdWithAString* fworld_matCmd;
	G4UIcmdWithAString* fgap_matCmd;
	G4UIcmdWithAString* fmonture_matCmd;
    G4UIcmdWithADoubleAndUnit* ftotal_diameterCmd;
    G4UIcmdWithADoubleAndUnit* ftotal_heightCmd;
    G4UIcmdWithADoubleAndUnit* fendCapToDetectorCmd;
    G4UIcmdWithADoubleAndUnit* fwindow_thicknessCmd;
	G4UIcmdWithADoubleAndUnit* fmonture_thicknessCmd;

	G4UIcmdWithADoubleAndUnit* fsourceTopCmd;
	G4UIcmdWithADoubleAndUnit* ffrontalBeamRadiusCmd;
	G4UIcmdWithADoubleAndUnit* ffrontalBeamAngleCmd;

	G4UIcmdWithAnInteger*     sourceTypeCmd;//link to source (i.e. PrimaryGeneratorAction)... it is the only way since Detector is initialized first!

	G4UIcmdWithAString* fsource_matCmd;
	G4UIcmdWithAString* fsource_envelope_matCmd;
	G4UIcmdWithADoubleAndUnit* fsource_total_heightCmd;
	G4UIcmdWithADoubleAndUnit* fsource_external_diameterCmd;
	G4UIcmdWithADoubleAndUnit* fsource_envelope_thicknessCmd;

	G4UIcmdWithADoubleAndUnit* fsource_internal_diameterCmd;
	G4UIcmdWithADoubleAndUnit* fsource_upper_heightCmd;

	G4UIcmdWithoutParameter*   UpdateCmd;
};

#endif

