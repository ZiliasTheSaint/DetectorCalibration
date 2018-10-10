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
// $Id: EventAction.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Dan Fulea, JUL 2013, National Institute of Public Health Bucharest, Cluj-Napoca Regional Centre, Romania

#include "EventAction.hh"

#include "RunAction.hh"
#include "EventActionMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleGun.hh"
#include "Randomize.hh"
#include <iomanip>


EventAction::EventAction()
{
  gun = (PrimaryGeneratorAction*)
                G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  runAct = (RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();
  eventMessenger = new EventActionMessenger(this);
  printModulo = 1000;
  energyThresholdForCounting=0.;
}

EventAction::~EventAction()
{
  delete eventMessenger;
}

void EventAction::BeginOfEventAction(const G4Event* evt)
{  
  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0) { 
    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
    //CLHEP::HepRandom::showEngineStatus();
}
 
 // initialisation per event
 TrackLAbs = 0.;
 EnergyDet=0.;

 phener=0.;
}

void EventAction::EndOfEventAction(const G4Event* evt)
{
  //init for efficiency 
   G4double weight=gun->GetParticleWeight();
   G4double deltaE=gun->GetDeltaE();
   if (phener>energyThresholdForCounting)//0.)
    runAct->addPcumPerEvent(weight);

   const G4ParticleGun* particleGun =gun->GetParticleGun();
   G4double particleEnergy = particleGun->GetParticleEnergy();

   //G4double restMass = particleGun->GetParticleDefinition()->GetPDGMass();
   //G4double kineticEnergy=particleEnergy-restMass;
   //particleEnergy IS ALREADY A KINETIC ENERGY!!!!
   if (phener>=particleEnergy-deltaE && phener<=particleEnergy+deltaE)
   //if (phener>=kineticEnergy-deltaE && phener<=kineticEnergy+deltaE)
   {
	   runAct->addPPeakPerEvent(weight);
   }

   if (phener>=particleEnergy-2.0*deltaE && phener<particleEnergy-deltaE)
   //if (phener>=kineticEnergy-2.0*deltaE && phener<kineticEnergy-deltaE)
   {
	   runAct->addPBkgPerEvent(weight);
   }
  //////////////////
  //accumulates statistic
  //fill per Event require 2 volumes to score. This is just for further modification of application
  runAct->fillPerEvent(EnergyDet, 0.0, TrackLAbs, 0.0);
  
  //print per event (modulo n)
  //
  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0) {
    G4cout << "---> End of event: " << evtNb << G4endl;	

    G4cout
		//first scoring volume:
       << "   Detector: total energy: " << std::setw(7)
                                        << G4BestUnit(EnergyDet,"Energy")
       << "       total track length: " << std::setw(7)
                                        << G4BestUnit(TrackLAbs,"Length")
       << G4endl;
      /*
	  //second scoring volume:
	   << "        Gap: total energy: " << std::setw(7)
                                        << G4BestUnit(EnergyGap,"Energy")
       << "       total track length: " << std::setw(7)
                                        << G4BestUnit(TrackLGap,"Length")
       << G4endl;*/
	  
  }
}  
