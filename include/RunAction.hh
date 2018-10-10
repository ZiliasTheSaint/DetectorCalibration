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
// $Id: RunAction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Dan Fulea, JUL 2013, National Institute of Public Health Bucharest, Cluj-Napoca Regional Centre, Romania

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class DetectorConstruction;
class PrimaryGeneratorAction;

class RunAction : public G4UserRunAction
{
public:
  RunAction();
  virtual ~RunAction();

  void BeginOfRunAction(const G4Run*);
  void   EndOfRunAction(const G4Run*);
  //default, two volume are used for scoring, named generically Abs and Gap.  
  void fillPerEvent(G4double, G4double, G4double, G4double);
  //-------------------------------------------------------
  G4double basicAnalyze(G4double x, G4double x2, G4int n){
	  G4double temp=x;G4double temp2=x2;G4int nevents=n;
	  temp=temp/nevents;temp2=temp2/nevents;
	  temp2 = temp2 - temp*temp;
	  if (nevents>1) temp2=temp2/(nevents-1.0);
	  if (temp2 >0.)
	   temp2 = sqrt(temp2); 
	  //else temp2 = 99.99;//never!
	  //=========percent
	  if (temp!=0.0){
	   temp2 = std::min(100.0*temp2/temp,99.9);
      } //else temp2 = 99.9;//no score means no score not necessarly an error!
	  return temp2;
  }
  //-----------------------------------------------------------
  void addPcumPerEvent(G4double w){pcum=pcum+w;pcum2=pcum2+w*w;};
  void addPPeakPerEvent(G4double w){ppeak=ppeak+w;ppeak2=ppeak2+w*w;};
  void addPBkgPerEvent(G4double w){pbkg=pbkg+w;pbkg2=pbkg2+w*w;};

private:
  DetectorConstruction*    Detector;     //pointer to the geometry
  PrimaryGeneratorAction* gun; //pointer to source
  //default, two volume are used for scoring, named generically Abs and Gap. 
  G4double sumEAbs, sum2EAbs;
  G4double sumEGap, sum2EGap;
  //for instance:
  //E reffers to Energy in those volumes while L reffers to total tracks of charged particles in those volumes.  
  G4double sumLAbs, sum2LAbs;
  G4double sumLGap, sum2LGap;    
  //=============
  G4double pcum, pcum2;//cumulated pulses
  G4double ppeak, ppeak2;//peak gross pulses
  G4double pbkg, pbkg2;//bkg pulses
};

#endif

