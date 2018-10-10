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
// $Id: EventActionMessenger.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Dan Fulea, JUL 2013, National Institute of Public Health Bucharest, Cluj-Napoca Regional Centre, Romania

#ifndef EventActionMessenger_h
#define EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class EventAction;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;


class EventActionMessenger: public G4UImessenger
{
public:
  EventActionMessenger(EventAction*);
  virtual ~EventActionMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  EventAction*          eventAction;
  G4UIdirectory*        eventDir;   
  G4UIcmdWithAnInteger* PrintCmd;
  G4UIcmdWithADoubleAndUnit* energyThresholdForCountingCmd;
};



#endif
