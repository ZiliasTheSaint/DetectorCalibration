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
// $Id: PrimaryGeneratorAction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Dan Fulea, JUL 2013, National Institute of Public Health Bucharest, Cluj-Napoca Regional Centre, Romania

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class DetectorConstruction;
class PrimaryGeneratorMessenger;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();    
  virtual ~PrimaryGeneratorAction();

  void GeneratePrimaries(G4Event*);
  void SetRndmFlag(G4String val) { rndmFlag = val;};

  // method to access particle gun
  //const G4ParticleGun* GetParticleGun() const { return particleGun; }
  G4ParticleGun* GetParticleGun() const { return particleGun; }
  G4double GetDeltaE() { return fdeltaE;};
  G4double GetParticleWeight() { return fweight;};  
  void SetDeltaE(G4double val){fdeltaE=val;};
  
private:
  G4ParticleGun*           particleGun;	 //pointer a to G4  class
  DetectorConstruction*    Detector;     //pointer to the geometry
    
  PrimaryGeneratorMessenger* gunMessenger;   //messenger of this class
  G4String                   rndmFlag;	     //flag for a rndm emmission, no variance reduction

  //=============================
  G4double fdeltaE;//delta in energy to define the peak!
  G4double fweight;//particle weight used in variance reduction
  G4double xin;G4double yin;G4double zin;//initial coordinates
  G4double uin;G4double vin;G4double win;//initial directional cosine
  void ProcessFrontalBeam();
  void ProcessPointSource();
  void ProcessSarpagan();
  void ProcessMarinelliSup(G4double,G4double,G4double);
  void ProcessMarinelliMiddle(G4double,G4double,G4double);
  void ProcessMarinelliInf(G4double,G4double,G4double);
};

#endif


