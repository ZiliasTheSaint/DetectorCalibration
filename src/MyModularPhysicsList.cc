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
// Dan Fulea, JUL 2013, National Institute of Public Health Bucharest, Cluj-Napoca Regional Centre, Romania
#include "globals.hh"

#include "MyModularPhysicsList.hh"

#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmStandardPhysics.hh"

#include "G4EmStandardPhysics_option3.hh"
#include "G4EmExtraPhysics.hh"

#include "G4IonBinaryCascadePhysics.hh"
#include "G4QStoppingPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"

#include "HadronPhysicsQGSP_BIC.hh"
#include "HadronPhysicsQGSP_BERT_HP.hh"

#include "G4HadronElasticPhysicsHP.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4QStoppingPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "EmStandardPhysics_option4.hh"

#include "G4HadronElasticPhysicsXS.hh"
#include "G4HadronInelasticQBBC.hh"
/*
Based on recommended QGSP_BERT_HP for "low energies", i.e. up to 10GeV
from 15GeV to 100TeV will use FTFP_BERT model!!!

Modified..overall good=>based on "brand new" QBBC list.
*/
MyModularPhysicsList::MyModularPhysicsList() 
: G4VModularPhysicsList()
{
  //SetVerboseLevel(1);
  G4int ver=1;

  defaultCutValue = 0.7*CLHEP::mm;  

  // EM Physics
  //RegisterPhysics(new G4EmStandardPhysics());//standard
  //option 3 (if geant4.9.5) or 4 (if geant4.9.6) will be use for allowing precise simulation at low and intermediate energies!!
  //RegisterPhysics( new G4EmStandardPhysics_option3(ver) );//in geant4.9.6 will use _option4!!
  //RegisterPhysics( new EmStandardPhysics_option4(ver) );//just a replica here!
  //RegisterPhysics( new G4EmStandardPhysics(ver) );
  RegisterPhysics( new G4EmStandardPhysics_option4(ver) );

  // Synchroton Radiation & GN Physics (GammaNuclear Physics)
  RegisterPhysics( new G4EmExtraPhysics(ver) );

  // Decays
  RegisterPhysics(new G4DecayPhysics(ver));

  // Hadron Elastic scattering
  //RegisterPhysics( new G4HadronElasticPhysicsHP(ver) );
  RegisterPhysics( new G4HadronElasticPhysicsXS(ver) );

  // Hadron Physics
  //RegisterPhysics( new HadronPhysicsQGSP_BERT_HP(ver));//>15GeV use FTFP_BERT
  RegisterPhysics( new G4HadronInelasticQBBC(ver));

  // Stopping Physics
  RegisterPhysics( new G4StoppingPhysics(ver));//NA In geant4.9.5
  //RegisterPhysics( new G4QStoppingPhysics(ver));

  // Ion Physics
  RegisterPhysics( new G4IonPhysics(ver));

  //Neutron tracking cut
  RegisterPhysics( new G4NeutronTrackingCut(ver));//only if FTFP in use!..or QBBC!

  // Radioactive decay
  RegisterPhysics(new G4RadioactiveDecayPhysics(ver));  
  
}


MyModularPhysicsList::~MyModularPhysicsList()
{ 
}


void MyModularPhysicsList::SetCuts()
{
  //G4VUserPhysicsList::SetCuts();
	SetCutsWithDefault();
}  
