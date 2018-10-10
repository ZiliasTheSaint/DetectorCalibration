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
// $Id: PrimaryGeneratorAction.cc,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Dan Fulea, JUL 2013, National Institute of Public Health Bucharest, Cluj-Napoca Regional Centre, Romania



#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  Detector = (DetectorConstruction*)
             G4RunManager::GetRunManager()->GetUserDetectorConstruction();  
  
  //create a messenger for this class
  gunMessenger = new PrimaryGeneratorMessenger(this);

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
                    = particleTable->FindParticle(particleName="gamma");//e-");//DEFAULT, change in run1.mac
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));//DEFAULT
  particleGun->SetParticleEnergy(0.662 *MeV);//DEFAULT
  
  rndmFlag = "off";//default
  
  fdeltaE=5.0*keV;//default
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}
//==================================================
void PrimaryGeneratorAction::ProcessFrontalBeam(){

  G4double radius = Detector->GetWorldSizeRadius();
  G4double height = Detector->GetWorldSizeZ();
  G4double beamRadius=Detector->GetFrontalBeamRadius();
  G4double teta =Detector->GetFrontalBeamAngle();
  
  xin=0.0;
  yin=0.0;
  zin=-0.5 * height;
  uin=0.0;
  vin=0.0;
  win=cos(teta);
  
  G4double R2=0.0;
  while(true){
	  xin=G4UniformRand();
	  xin=(2.0 * xin - 1.0)*beamRadius;
	  yin=G4UniformRand();
	  yin=(2.0 * yin - 1.0)*beamRadius;
	  R2=xin*xin+yin*yin;

	  if(R2<=beamRadius*beamRadius && R2>0.0){
		  //shout at random but not good for practical purpose
		  //double cosphi=xin/sqrt(R2);
		  //double sinphi=yin/sqrt(R2);
		  //uin=sin(teta)*cosphi;
		  //vin=sin(teta)*sinphi;

		  //fixed parallel beam with a bias angle, say with respect to x-axis
		  //it is ok due to symmetry of the geometry!! 
		  vin=0.0;
		  uin=sqrt(1.0-win*win-vin*vin);

		  //now the weight
		  fweight=1.0;//here, no var reduction, source is indeed a gun!
		  break;
	  }
  }
}
//====================================================
void PrimaryGeneratorAction::ProcessPointSource(){

  G4double height = Detector->GetWorldSizeZ();  
  G4double diam = Detector->GetWorldSizeRadius();//here, its radius
  G4double distz=Detector->GetPointSourceOrFrontalBeamToDetectorDistance();
  //init-------------
  double pi =3.14159265359;
  xin=0.0;
  yin=0.0;
  zin=0.0;
  uin=0.0;
  vin=0.0;
  win=1.0;
  //-------------POLAR ANGLE-----------
  G4double costmax = distz/sqrt(distz * distz + diam * diam);
  //emmision probability for angles less than theta max:
  G4double dom = (1.0 - costmax) / 2.0;// >0 and <1/2, costmax <1 and >0, thetamax<90!
  //now a random angle less than thetamax
  G4double r=G4UniformRand();
  r = r * dom;//assure costet < costmax
  G4double costet = 1.0 - 2.0 * r;//>0, positive z axis; 2*r-1;//<0 i.e. negativ z axis!!
  //-----------------------------
  G4double sintet = sqrt(1.0 - costet * costet);
  G4double tgtet = sintet / costet;
  G4double teta = abs(atan(tgtet));// -pi/2,pi/2; teta is polar angle and should be >0!
  //----------AZIMUTHAL ANGLE----------  
  r = G4UniformRand();
  G4double phi2 = 2.0 * pi * r;
  //G4cout<<"pi= "<<pi<<"!!!!!!!!@@@@@@@@@@@@@@@@!!!!!!!!"<<G4endl;//YESSS!!!!!!!!! both pi are good!!!!
  uin=sintet * cos(phi2);//sin(teta) * cos(phi2); IN LINUX THIS IS ALWAYS 0!! Must be a c++ compiler bug in linux!
  vin=sintet * sin(phi2);//sin(teta) * sin(phi2); Therefore we must change to SINTET!!!!!!!!!!!
  win = costet;

  xin = 0.0;
  yin = 0.0;
  zin=-0.5 * height;

  //maximum solid angle allowed for emission sampling
  G4double us1 = 2.0 * pi	* (1.0 - costmax);
  fweight = us1 / (4.0 * pi);

  //////////NO VARIANCE REDUCTION//////////
  if (rndmFlag == "on"){
	r = G4UniformRand();
	win = 2.0*r-1.0;//-1,+1
    G4double sinteta=sqrt(1-win*win);
    r = G4UniformRand();
    G4double phi = 2 * pi * r;
    uin = sinteta*cos(phi);
    vin = sinteta*sin(phi);
    fweight=1.0;
  }

}
//====================================================
void PrimaryGeneratorAction::ProcessSarpagan(){

  G4double hdet=Detector->GetDetectorHeight();
  G4double hsource=Detector->GetSourceTotalHeight();
  G4double ethick=Detector->GetSourceEnvelopeThickness();
  G4double asource=0.5*(Detector->GetSourceExternalDiameter());//radius
  G4double adet=0.5*(Detector->GetDetectorDiameter());//radius
  //init-------------
   double pi =3.14159265359;
  //G4cout<<"pi= "<<pi<<"!!!!!!!!@@@@@@@@@@@@@@@@!!!!!!!!"<<G4endl;//YESSS!!!!!!!!! both pi are good!!!!

  xin=0.0;
  yin=0.0;
  zin=0.0;
  uin=0.0;
  vin=0.0;
  win=1.0;
  //===================
  G4double r = 0.0;
  while (true) {
	r = G4UniformRand();
	if (G4UniformRand() < 0.5)
		r = 1.0 - r;
	if (r != 1.0)
		break;
  }
  G4double z1 = (hsource - ethick) * r;// in source:0->hs-ethick
  G4double distz = hsource - z1; 

  r = G4UniformRand();
  G4double ro1 = asource * sqrt(r);// dist to z axis ro1:
  r = G4UniformRand();
  G4double phi1 = 2.0 * pi * r;

  G4double diam = 0.0;  
  if (asource <= adet) {
	diam = 2.0 * adet;
  } else {
	diam = 2.0 * asource;
  }
  //-------------POLAR ANGLE-----------
  G4double costmax = distz/sqrt(distz * distz + diam * diam);
  //emmision probability for angles less than theta max:
  G4double dom = (1.0 - costmax) / 2.0;// >0 and <1/2, costmax <1 and >0, thetamax<90!
  //now a random angle less than thetamax
  r=G4UniformRand();
  r = r * dom;//assure costet < costmax
  G4double costet = 1.0 - 2.0 * r;//>0, positive z axis; 2*r-1;//<0 i.e. negativ z axis!!
  //-----------------------------
  G4double sintet = sqrt(1.0 - costet * costet);
  G4double tgtet = sintet / costet;
  G4double teta = abs(atan(tgtet));// -pi/2,pi/2; teta is polar angle and should be >0!
  //----------AZIMUTHAL ANGLE----------
  r = G4UniformRand();
  G4double phi2 = 2 * pi * r;
  
  G4double height=hdet+hsource;

  //uin=sin(teta) * cos(phi2);
  //vin=sin(teta) * sin(phi2);
  uin=sintet * cos(phi2);//sin(teta) * cos(phi2); IN LINUX THIS IS ALWAYS 0!! Must be a c++ compiler bug in linux!
  vin=sintet * sin(phi2);//sin(teta) * sin(phi2); Therefore we must change to SINTET!!!!!!!!!!!
  win = costet;

  xin = ro1*cos(phi1);
  yin = ro1*sin(phi1);
  zin=-0.5 * height+z1;

  //maximum solid angle allowed for emission sampling
  G4double us1 = 2.0 * pi	* (1.0 - costmax);
  fweight = us1 / (4.0 * pi);

  //////////NO VARIANCE REDUCTION//////////
  if (rndmFlag == "on"){
	r = G4UniformRand();
	win = 2.0*r-1.0;//-1,+1
    G4double sinteta=sqrt(1-win*win);
    r = G4UniformRand();
    G4double phi = 2 * pi * r;
    uin = sinteta*cos(phi);
    vin = sinteta*sin(phi);
    fweight=1.0;
  }

}
//====================================================
void PrimaryGeneratorAction::ProcessMarinelliSup(G4double x0, G4double y0, G4double z0)
{
	G4double hsourceup=Detector->GetSourceUpperHeight();
	G4double asource=0.5*Detector->GetSourceExternalDiameter();
	G4double hsource=Detector->GetSourceTotalHeight();
    G4double hdettot=Detector->GetDetectorHeight();

	G4double height = std::max(hdettot+hsourceup,hsource);

	double pi =3.14159265359;
    //G4cout<<"pi= "<<pi<<"!!!!!!!!@@@@@@@@@@@@@@@@!!!!!!!!"<<G4endl;//YESSS!!!!!!!!! both pi are good!!!!

	G4double distz = hsourceup - z0;
	G4double diam = 2.0 * asource;
	G4double costmax = distz/sqrt(distz * distz + diam * diam);
	G4double us1 = 2.0 * pi	* (1.0 - costmax);

	G4double dom = (1.0 - costmax) / 2.0;
	G4double r=G4UniformRand();
    r = r * dom;//assure costet < costmax
    G4double costet = 1.0 - 2.0 * r;//>0, positive z axis; 2*r-1;//<0 i.e. negativ z axis!!
    //-----------------------------
    G4double sintet = sqrt(1.0 - costet * costet);
    G4double tgtet = sintet / costet;
    G4double teta = abs(atan(tgtet));
	
	r = G4UniformRand();
    G4double phi2 = 2 * pi * r;

	//uin=sin(teta) * cos(phi2);
    //vin=sin(teta) * sin(phi2);
	uin=sintet * cos(phi2);//sin(teta) * cos(phi2); IN LINUX THIS IS ALWAYS 0!! Must be a c++ compiler bug in linux!
    vin=sintet * sin(phi2);//sin(teta) * sin(phi2); Therefore we must change to SINTET!!!!!!!!!!!
    win = costet;

	xin = x0;
    yin = y0;
    zin=-0.5 * height+z0;
	//--------------------------
	fweight = us1 / (4.0 * pi);

	//////////NO VARIANCE REDUCTION//////////
  if (rndmFlag == "on"){
	r = G4UniformRand();
	win = 2.0*r-1.0;//-1,+1
    G4double sinteta=sqrt(1-win*win);
    r = G4UniformRand();
    G4double phi = 2 * pi * r;
    uin = sinteta*cos(phi);
    vin = sinteta*sin(phi);
    fweight=1.0;
  }
}

void PrimaryGeneratorAction::ProcessMarinelliMiddle(G4double x0, G4double y0, G4double z0)
{
	G4double hsourceup=Detector->GetSourceUpperHeight();	
	G4double hsource=Detector->GetSourceTotalHeight();
    G4double hdettot=Detector->GetDetectorHeight();
	G4double ethick=Detector->GetSourceEnvelopeThickness();
	G4double bsource=0.5*(Detector->GetSourceInternalDiameter())+ethick;
	G4double height = std::max(hdettot+hsourceup,hsource);
	G4double adettot=0.5*Detector->GetDetectorDiameter();
	//--------------------------------------------------------------
	double pi =3.14159265359;
    //G4cout<<"pi= "<<pi<<"!!!!!!!!@@@@@@@@@@@@@@@@!!!!!!!!"<<G4endl;//YESSS!!!!!!!!! both pi are good!!!!

	G4double distz =  sqrt(x0*x0+y0*y0)-bsource;
	G4double diam = 2.0 * sqrt(adettot*adettot + (hdettot+ethick) * (hdettot+ethick) / 4.0);
	G4double costmax = distz/sqrt(distz * distz + diam * diam);
	G4double us1 = 2.0 * pi	* (1.0 - costmax);
	
	G4double dom = (1.0 - costmax) / 2.0;
	G4double r=G4UniformRand();
    r = r * dom;//assure costet < costmax
    G4double costet = 1.0 - 2.0 * r;//for instance

	G4double sintet = sqrt(1.0 - costet * costet);
    G4double tgtet = sintet / costet;
    G4double teta = abs(atan(tgtet));//// -pi/2,pi/2	

	r = G4UniformRand();
    G4double phi2 = 2 * pi * r;
	//=================
	G4double distance=sqrt(x0*x0+y0*y0)/costet;//>0
	//Assume we have initial directional cosines:
	G4double u=-x0/distance;
	G4double v=-y0/distance;
	G4double w=0.0;
	//i.e. particle on xy-plane moving against positive x and y axis (towards detector)
	//now the particle suffers a scattering process, having polar angle teta and a random azimuthal angle phi2, relative to its own direction!
	//we have (from theory) the directional cosines relative to lab:

	if (abs(w) > 0.99999)// miuz->normal incident, never happen here!
	{
		//uin = sin(teta) * cos(phi2);// new miux
		//vin = sin(teta) * sin(phi2);// new miuy
		uin=sintet * cos(phi2);//sin(teta) * cos(phi2); IN LINUX THIS IS ALWAYS 0!! Must be a c++ compiler bug in linux!
		vin=sintet * sin(phi2);//sin(teta) * sin(phi2); Therefore we must change to SINTET!!!!!!!!!!!
		if (w < 0.0)
			win = -costet;//win = -cos(teta);
		else
			win = costet;//win = cos(teta);// new miuz
	} else {
		G4double temp = sqrt(1.0 - w * w);
		//uin = sin(teta) * (u * w * cos(phi2) - v * sin(phi2))/ temp + u * cos(teta);
		//vin = sin(teta) * (v * w * cos(phi2) + u * sin(phi2))/ temp + v * cos(teta);
		//win = -sin(teta) * cos(phi2) * temp + w * cos(teta);

		uin = sintet * (u * w * cos(phi2) - v * sin(phi2))/ temp + u * costet;
		vin = sintet * (v * w * cos(phi2) + u * sin(phi2))/ temp + v * costet;
		win = -sintet * cos(phi2) * temp + w * costet;
	}
	
	xin = x0;
    yin = y0;
    zin=-0.5 * height+z0;
	//--------------------------
	fweight = us1 / (4.0 * pi);

	//////////NO VARIANCE REDUCTION//////////
  if (rndmFlag == "on"){
	r = G4UniformRand();
	win = 2.0*r-1.0;//-1,+1
    G4double sinteta=sqrt(1-win*win);
    r = G4UniformRand();
    G4double phi = 2 * pi * r;
    uin = sinteta*cos(phi);
    vin = sinteta*sin(phi);
    fweight=1.0;
  }
}

void PrimaryGeneratorAction::ProcessMarinelliInf(G4double x0, G4double y0, G4double z0)
{
	G4double hsourceup=Detector->GetSourceUpperHeight();
	G4double asource=0.5*Detector->GetSourceExternalDiameter();
	G4double hsource=Detector->GetSourceTotalHeight();
    G4double hdettot=Detector->GetDetectorHeight();

	G4double height = std::max(hdettot+hsourceup,hsource);

	double pi =3.14159265359;
  //G4cout<<"pi= "<<pi<<"!!!!!!!!@@@@@@@@@@@@@@@@!!!!!!!!"<<G4endl;//YESSS!!!!!!!!! both pi are good!!!!

	G4double distz = abs(z0 - (hsourceup + hdettot));
	G4double diam = 2.0 * asource;
	G4double costmax = distz/sqrt(distz * distz + diam * diam);
	G4double us1 = 2.0 * pi	* (1.0 - costmax);

	G4double dom = (1.0 - costmax) / 2.0;
	G4double r=G4UniformRand();
    r = r * dom;//assure costet < costmax
    G4double costet = -1.0 + 2.0 * r;//<0, 2*r-1;//<0 i.e. negativ z axis!!

	G4double sintet = sqrt(1.0 - costet * costet);
    G4double tgtet = sintet / costet;
    G4double teta = abs(atan(tgtet));

	r = G4UniformRand();
    G4double phi2 = 2 * pi * r;

	//uin=sin(teta) * cos(phi2);
    //vin=sin(teta) * sin(phi2);
	uin=sintet * cos(phi2);//sin(teta) * cos(phi2); IN LINUX THIS IS ALWAYS 0!! Must be a c++ compiler bug in linux!
	vin=sintet * sin(phi2);//sin(teta) * sin(phi2); Therefore we must change to SINTET!!!!!!!!!!!
    win = costet;

	xin = x0;
    yin = y0;
    zin=-0.5 * height+z0;
	//--------------------------
	fweight = us1 / (4.0 * pi);

	//////////NO VARIANCE REDUCTION//////////
  if (rndmFlag == "on"){
	r = G4UniformRand();
	win = 2.0*r-1.0;//-1,+1
    G4double sinteta=sqrt(1-win*win);
    r = G4UniformRand();
    G4double phi = 2 * pi * r;
    uin = sinteta*cos(phi);
    vin = sinteta*sin(phi);
    fweight=1.0;
  }
}
////////////===============
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  
  G4int sourceType=Detector->GetSourceTypeCode();
  
  G4double hsource=Detector->GetSourceTotalHeight();
  G4double hsourceup=Detector->GetSourceUpperHeight();
  G4double hdettot=Detector->GetDetectorHeight();
  G4double ethick=Detector->GetSourceEnvelopeThickness();
  G4double asource=0.5*Detector->GetSourceExternalDiameter();
  G4double bsource=0.5*(Detector->GetSourceInternalDiameter())+ethick;

  G4double height = hsource;
  G4int whereIam=-1;

  G4bool infcyl = false;// test if inf cylinder exists
  if (hsource > hsourceup + hdettot)
		infcyl = true;

  if (sourceType==0)
	ProcessFrontalBeam();
  else if (sourceType==1)
    ProcessPointSource();
  else if (sourceType==2)
    ProcessSarpagan();
  else if (sourceType==3){//Marinelli
	  //GES/EGS method for Marinelli sampling=not good since volumes are not equal so
	  //uniform sampling for z1 is not aqurate!
    //pretend we have an outer dummy cylinder in which we sample coordonates then check if we are in Marinelli baker!
    G4double x=0.0;
	G4double y=0.0;
	G4double z=0.0;

	double pi =3.14159265359;
    //G4cout<<"pi= "<<pi<<"!!!!!!!!@@@@@@@@@@@@@@@@!!!!!!!!"<<G4endl;//YESSS!!!!!!!!! both pi are good!!!!

	while(true){
	 
	 G4double r = G4UniformRand();
     G4double ro1 = asource * sqrt(r);
 	 r = G4UniformRand();
     G4double phi1 = 2 * pi * r;
	 x = ro1*cos(phi1);
     y = ro1*sin(phi1);
	 r = G4UniformRand();
     z=height*r;
	 G4double R2=x*x+y*y;
	 //now check
	
	 if (z <= (hsourceup - ethick)) {
		whereIam=0;//ProcessMCylSup
		break;
	 } else if (z > hsourceup) {
		if (infcyl) {
			if (z > hsourceup + hdettot)
			{
				if (bsource<=R2 && R2<=asource){
				   whereIam=2;//ProcessMCylInf
				   break;
				}
			} else// middle region
			{
				if (bsource<=R2 && R2<=asource){
				 whereIam=1;//ProcessMMidInf
				 break;
				}
			}
		} else// middle region
		{
			if (bsource<=R2 && R2<=asource){
		     whereIam=1;//ProcessMMidInf
			 break;
			}
		}
	 }//if (z > hsourceup)
	 else{
		 //top plastic
		 if (bsource<=R2 && R2<=asource){
		     whereIam=1;//ProcessMMidInf
			 break;
		 }
	 }
	}//while

	if (whereIam==0) ProcessMarinelliSup(x,y,z);
	else if (whereIam==1) ProcessMarinelliMiddle(x,y,z);
	else if (whereIam==2) ProcessMarinelliInf(x,y,z);
		
  }//end Marinelli
  else
	ProcessFrontalBeam();//default

  particleGun->SetParticlePosition(G4ThreeVector(xin,yin,zin));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(uin,vin,win));
 
  particleGun->GeneratePrimaryVertex(anEvent);
}
