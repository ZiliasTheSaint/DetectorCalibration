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
// $Id: EventActionMessenger.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Dan Fulea, JUL 2013, National Institute of Public Health Bucharest, Cluj-Napoca Regional Centre, Romania

#include "EventActionMessenger.hh"

#include "EventAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

EventActionMessenger::EventActionMessenger(EventAction* EvAct)
:eventAction(EvAct)
{
  eventDir = new G4UIdirectory("/gdet/event/");
  eventDir->SetGuidance("event control");
   
  PrintCmd = new G4UIcmdWithAnInteger("/gdet/event/printModulo",this);
  PrintCmd->SetGuidance("Print events modulo n");
  PrintCmd->SetParameterName("EventNb",false);
  PrintCmd->SetRange("EventNb>0");

  energyThresholdForCountingCmd=new G4UIcmdWithADoubleAndUnit("/gdet/event/energyThresholdForCounting",this);
  energyThresholdForCountingCmd->SetGuidance("Define energy threshold for counting.");
  energyThresholdForCountingCmd->SetParameterName("Energy",false);
  energyThresholdForCountingCmd->SetRange("Energy>=0.");
  energyThresholdForCountingCmd->SetUnitCategory("Energy");
  energyThresholdForCountingCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

EventActionMessenger::~EventActionMessenger()
{
  delete PrintCmd;
  delete eventDir;
  delete energyThresholdForCountingCmd;
}

void EventActionMessenger::SetNewValue(
                                        G4UIcommand* command,G4String newValue)
{ 
  if(command == PrintCmd)
    {eventAction->SetPrintModulo(PrintCmd->GetNewIntValue(newValue));}

  if( command == energyThresholdForCountingCmd ) {
    eventAction
      ->SetEnergyThresholdForCounting(energyThresholdForCountingCmd->GetNewDoubleValue(newValue));
  }
}
