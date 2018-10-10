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
// $Id: RunAction.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// Dan Fulea, JUL 2013, National Institute of Public Health Bucharest, Cluj-Napoca Regional Centre, Romania

#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleGun.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

RunAction::RunAction()
{
  // add new units for dose 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);        
  Detector = (DetectorConstruction*)
             G4RunManager::GetRunManager()->GetUserDetectorConstruction();  
  gun = (PrimaryGeneratorAction*)
                G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
}

RunAction::~RunAction()
{}

void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4ParticleGun* particleGun =gun->GetParticleGun();
  G4String particleName = particleGun->GetParticleDefinition()->GetParticleName();
  G4double particleEnergy = particleGun->GetParticleEnergy();
  G4double deltaE=gun->GetDeltaE();
  
  //energy fix <=> BUT NO FIX IS NEEDED since G4ParticleGun setParticleEnergy REQUIRES KINETIC ENERGY...see geant4 source!   
  //G4double restMass = particleGun->GetParticleDefinition()->GetPDGMass();
  //G4double totalEnergy=particleEnergy+restMass;
  //particleGun->SetParticleEnergy(totalEnergy);
  //////////////////////
  //===========================
  G4cout
     << "\n--------------------Start of Run------------------------------\n" <<G4endl;
  G4cout << "### Run " << aRun->GetRunID() << " start. Particle name: "
	  <<particleName<< "; Kinetic energy: "<< G4BestUnit(particleEnergy,"Energy") << G4endl;
  G4String sourceS=Detector->GetSourceType();
  //G4String sourceS2=gun->GetSourceType();
  G4cout <<"Detector -> source summary: " <<sourceS<<G4endl;
  //G4cout <<"Gun -> source summary: " <<sourceS2<<G4endl;
  G4cout <<"Peak window : " <<G4BestUnit(deltaE,"Energy") <<G4endl;
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
    
  //initialize cumulative quantities
  //
  sumEAbs = sum2EAbs =sumEGap = sum2EGap = 0.;
  sumLAbs = sum2LAbs =sumLGap = sum2LGap = 0.; 

  pcum=0.0;ppeak=0.0;pbkg=0.0;
  pcum2=0.0;ppeak2=0.0;pbkg2=0.0;
}

void RunAction::fillPerEvent(G4double EAbs, G4double EGap,
                                  G4double LAbs, G4double LGap)
{
  //accumulate statistic
  //
  sumEAbs += EAbs;  sum2EAbs += EAbs*EAbs;
  sumEGap += EGap;  sum2EGap += EGap*EGap;
  
  sumLAbs += LAbs;  sum2LAbs += LAbs*LAbs;
  sumLGap += LGap;  sum2LGap += LGap*LGap;  
}

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
  
  //compute statistics: mean and rms
  //
  G4double err = basicAnalyze(sumEAbs, sum2EAbs, NbOfEvents);
  sumEAbs /= NbOfEvents; sum2EAbs /= NbOfEvents;
  G4double rmsEAbs =err;
 
  err = basicAnalyze(sumLAbs, sum2LAbs, NbOfEvents);
  sumLAbs /= NbOfEvents; sum2LAbs /= NbOfEvents;
  G4double rmsLAbs =err;
 
  //sumEGap /= NbOfEvents; sum2EGap /= NbOfEvents;
  //G4double rmsEGap = sum2EGap - sumEGap*sumEGap;
  //if (rmsEGap >0.) rmsEGap = std::sqrt(rmsEGap); else rmsEGap = 0.;  
  //sumLGap /= NbOfEvents; sum2LGap /= NbOfEvents;
  //G4double rmsLGap = sum2LGap - sumLGap*sumLGap;
  //if (rmsLGap >0.) rmsLGap = std::sqrt(rmsLGap); else rmsLGap = 0.;
  
  //===================
  G4double mass = Detector->GetDetector_logical()->GetMass();
  G4double dose = sumEAbs/mass;
  G4double rmsDose = rmsEAbs;//%///mass;

  //======================eff========================
  err = basicAnalyze(pcum, pcum2, NbOfEvents);
  G4double efficiency_global = pcum * 100.0 / NbOfEvents;// %
  G4double efficiency_global_error = err;

  err = basicAnalyze(ppeak, ppeak2, NbOfEvents);
  ppeak2 = err;
  err = basicAnalyze(pbkg, pbkg2, NbOfEvents);
  pbkg2 = err;
  // "subtract background from peak"
  ppeak=ppeak-pbkg;
  // "estimate uncertainty on this subtracted value"
  ppeak2 = (ppeak2*ppeak/NbOfEvents/100.)*(ppeak2*ppeak/NbOfEvents/100.)
			+ (pbkg2*pbkg/NbOfEvents/100.)*(pbkg2*pbkg/NbOfEvents/100.)
			+ 2/(NbOfEvents*NbOfEvents)*(ppeak*pbkg)/(NbOfEvents-1.0);
  if (ppeak2 > 0.) {
	ppeak2=sqrt(ppeak2)/(ppeak/NbOfEvents);
  }

  if (ppeak2>0.999){
	ppeak2= 0.999;
  }
  //percent fix
  ppeak2 = 100.*ppeak2;
  //=============
  G4double efficiency = ppeak * 100.0 / NbOfEvents;// %
  G4double efficiency_error = ppeak2;

  //====================
  const G4ParticleGun* particleGun =gun->GetParticleGun();
  G4String particleName = particleGun->GetParticleDefinition()->GetParticleName();
  G4double particleEnergy = particleGun->GetParticleEnergy();
  G4String sourceS=Detector->GetSourceType();//gun->GetSourceType();
  G4double deltaE=gun->GetDeltaE();
  G4double mm = particleGun->GetParticleDefinition()->GetPDGMass();
  //===========================

  G4cout
     << "\n--------------------End of Run------------------------------\n" <<G4endl;
  //G4cout << " Particle name: "  <<particleName<< "; Initial total energy: "<< G4BestUnit(particleEnergy,"Energy") << G4endl;
  G4cout << " Particle name: "  <<particleName<< "; Initial kinetic energy: "<< G4BestUnit(particleEnergy,"Energy") << G4endl;
  G4cout << " Source geometry summary: "  <<sourceS << G4endl;
  G4cout << " Particle mass (rest energy): "  <<G4BestUnit(mm,"Energy") << G4endl;
  G4cout <<" Peak window : " <<G4BestUnit(deltaE,"Energy") <<G4endl;

  G4cout
     //<< "\n--------------------End of Run------------------------------\n"
     << "\n mean Energy in Detector : " << G4BestUnit(sumEAbs,"Energy")
     << " +- "                          <<rmsEAbs<<" %"//<< G4BestUnit(rmsEAbs,"Energy")  
     //<< "\n mean Energy in Gap      : " << G4BestUnit(sumEGap,"Energy")
     //<< " +- "                          << G4BestUnit(rmsEGap,"Energy")
     << G4endl;
     
  G4cout
     << "\n mean trackLength in Detector : " << G4BestUnit(sumLAbs,"Length")
     << " +- "                               <<rmsLAbs<<" %"//<< G4BestUnit(rmsLAbs,"Length")  
     //<< "\n mean trackLength in Gap      : " << G4BestUnit(sumLGap,"Length")
     //<< " +- "                               << G4BestUnit(rmsLGap,"Length")
     << "\n------------------------------------------------------------\n"
     << G4endl;

  G4cout
     << "\n dose in Detector : " << G4BestUnit(dose,"Dose")
     << " +- "                          <<rmsDose<<" %"//<< G4BestUnit(rmsDose,"Dose")  
	 << "\n------------------------------------------------------------\n"
     << G4endl;

  G4cout
     << "\n peak efficiency in Detector (%): " << efficiency
     << " +- "                          <<efficiency_error<<" %"
	 //<< "\n------------------------------------------------------------\n"
     << G4endl;
  G4cout
     << "\n global efficiency in Detector (%): " << efficiency_global
     << " +- "                          <<efficiency_global_error<<" %"
	 << "\n------------------------------------------------------------\n"
     << G4endl;
}
