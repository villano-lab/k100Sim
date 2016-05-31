/*==================k100_PrimaryGeneratorAction.cc================= 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE: Code supporting the generators for the k100 simulation
               all user-defined particle sources should be instantiated
	       in the constructor and executed in GeneratePrimaries(*).
              
              
======================================================================*/

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "G4Neutron.hh"

#include "Randomize.hh"
#include <math.h>

#include "k100_PrimaryGeneratorAction.hh"
#include "k100_DetectorConstruction.hh"
#include "k100_EventInfo.hh"

k100_PrimaryGeneratorAction::k100_PrimaryGeneratorAction()
{
  
  //Get the geometry
  Ndet = (k100_DetectorConstruction*)(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  //set basic particle gun
  particleGun = new G4ParticleGun();

  //set overall rotation to be consistent with geometry
  particleGun->SetParticlePosition(G4ThreeVector(0.0*m,0.0*m,-120.0*cm));
  G4ThreeVector row1 = G4ThreeVector(1,0,0);
  G4ThreeVector row2 = G4ThreeVector(0,0,1.0);
  G4ThreeVector row3 = G4ThreeVector(0,-1.0,0);
  xrot = new G4RotationMatrix(row1,row2,row3);
  xrot->setRows(row1,row2,row3);

}
k100_PrimaryGeneratorAction::~k100_PrimaryGeneratorAction()
{
  delete particleGun;
}
void k100_PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //make an ancillary EventInfo container to store some event info
  //and register it with this event
  k100_EventInfo *anEventInfo = new k100_EventInfo();
  anEvent->SetUserInformation(anEventInfo);

  //use below if you use a custom particle generator which
  //exits based on some condition (other than total number
  //of primaries thrown)
  /*if(!particleGun->GetHasExited())
   particleGun->GeneratePrimaryVertex(anEvent);
  else
  {*/
    // RunManager cannot abort the event from inside
    // UserGeneratePrimaries(), so we do a soft abort
    // to the RunManager, and abort the event ourselves.
    // The result is the same as a hard abort.
  /*     G4cout << "about to abort" << G4endl;
       G4RunManager::GetRunManager()->AbortRun(true); 
       anEvent->SetEventAborted();
  }*/

   particleGun->SetParticleDefinition(G4Neutron::Definition()); 
   std::vector<G4double> adirection = GenerateRandomDirection();
   particleGun->SetParticleMomentumDirection(G4ThreeVector(adirection[0],adirection[1],adirection[2]));
   particleGun->SetParticleEnergy(8.0);
   particleGun->GeneratePrimaryVertex(anEvent);
  
}
std::vector<G4double> k100_PrimaryGeneratorAction::GenerateRandomDirection()
{
  std::vector<G4double> angles;
  G4double rand1,rand2;
  rand1=G4UniformRand();
  rand2=G4UniformRand();
  angles.push_back(rand2*2*pi);            //phi
  G4double costhet = 2*rand1-1.0;
  angles.push_back(costhet);                 //costhet
  angles.push_back(sqrt(1 - costhet*costhet)); //sinthet

  return angles;

}
