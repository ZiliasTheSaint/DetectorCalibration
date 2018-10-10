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
// $Id: EventAction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Dan Fulea, JUL 2013, National Institute of Public Health Bucharest, Cluj-Napoca Regional Centre, Romania

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class RunAction;
class EventActionMessenger;
class PrimaryGeneratorAction;

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  virtual ~EventAction();

  void  BeginOfEventAction(const G4Event*);
  void  EndOfEventAction(const G4Event*);
    
  void AddEnergyInDetector(G4double de, G4double dl) {EnergyDet += de; TrackLAbs += dl;};
  void AddEnergyInDetector(G4double de) {EnergyDet += de;};
  void AddRealEnergy(G4double de){phener +=de;};

  //if second volume is required, the next line must be uncommented and modified accordingly
  //void AddEnergyInVolume2(G4double de, G4double dl) {EnergyGap += de; TrackLGap += dl;};
                         
  void SetPrintModulo(G4int    val)  {printModulo = val;};
  void SetEnergyThresholdForCounting (G4double val){energyThresholdForCounting=val;};
    
private:
   RunAction*  runAct;

   G4double  EnergyDet;//store energy deposition in volume=detector volume
   G4double  TrackLAbs;//store total track lenght (sum of all tracks; one track is itself a sum of subsequent steps 
   //which are required by physics) of charged particles in volume=detector volume
   G4double  phener;
   G4double  energyThresholdForCounting;//a virtual index to handle alpha-beta measurements.
   //If plateau is not an ideal straight (flat) line, then measure a known source at about end of plateau voltage.
   //Calculate efficiency. At this stage, the efficiency should match with the Monte-Carlo simulated efficiency,
   //because all quanta that reach the detector (without impediment) is supposed to be recorded (any small EDEP=count).
   //Then, set working voltage at your leisure...typically at half of plateau. Calculate again the efficiency.
   //Adjust accordingly energyThresholdForCounting to match. In this way, one can measure any source in any 
   //geometry and composition by knowing energyThresholdForCounting at working voltage!!

   G4int     printModulo;
                             
   EventActionMessenger*  eventMessenger;
   PrimaryGeneratorAction* gun;
};

#endif

    
