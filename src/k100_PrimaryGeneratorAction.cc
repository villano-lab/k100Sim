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
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "G4Neutron.hh"

#include "Randomize.hh"
#include <math.h>
#include <iostream>
#include "k100_PrimaryGeneratorAction.hh"
#include "k100_DetectorConstruction.hh"
#include "k100_EventInfo.hh"


k100_PrimaryGeneratorAction::k100_PrimaryGeneratorAction(G4bool useCapture, G4String InFile):
throwCaptures(useCapture)
{
  
  //Get the geometry
  Ndet = (k100_DetectorConstruction*)(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  //set basic particle gun
  particleGun = new G4ParticleGun();

  //set up a more intricate generator
  particleSource = new G4GeneralParticleSource();

  //set overall rotation to be consistent with geometry
  particleGun->SetParticlePosition(G4ThreeVector(0.0*m,0.0*m,-120.0*cm));
  G4ThreeVector row1 = G4ThreeVector(1,0,0);
  G4ThreeVector row2 = G4ThreeVector(0,0,1.0);
  G4ThreeVector row3 = G4ThreeVector(0,-1.0,0);
  xrot = new G4RotationMatrix(row1,row2,row3);
  xrot->setRows(row1,row2,row3);
  
  if(useCapture) {
    if(InFile == "none") {
      std::cout<<"No input file provided for custom particle gun. Exiting"<<std::endl;
      exit(0);
    }
    // std::cout<<"YOLO ============ Can I do it??  ========="<<std::endl;
    // std::cout<<"YOLO ============ I am doing it  ========="<<std::endl;
    // nrCascadeTree = new k100_nrCascadeTree("/Users/shubhampandey/work/k100sim_analysis/combine_k100sim_nrCascade/python_code/check.root");
    nrCascadeTree = new k100_nrCascadeTree(InFile);
    if(!nrCascadeTree) {
      std::cout<<"nrCascadeTree not loaded. Exiting."<<std::endl;
      exit(0);
    }
    else {
      std::cout<<"nrCascadeTree loaded successfully with total entries = "<<nrCascadeTree->tree->GetEntries()<<std::endl;
    }

  }

}
k100_PrimaryGeneratorAction::~k100_PrimaryGeneratorAction()
{
  delete particleGun;
  delete particleSource;
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

   if(throwCaptures){

     // particleGun->SetParticleDefinition(G4Neutron::Definition()); 
     // //std::vector<G4double> angles = GenerateRandomDirection();
     // //particleGun->SetParticleMomentumDirection(G4ThreeVector(angles[2]*cos(angles[0]),angles[2]*sin(angles[0]),angles[1]));
     // //spandey
     // //particleGun->SetParticleEnergy(8.0);
     // particleGun->SetParticleEnergy(2.5e-8);
     // particleGun->SetParticlePosition(G4ThreeVector(0.0*m,0.0*m,-1.0*cm));
     // particleGun->SetParticleMomentumDirection(G4ThreeVector(-1,0,0));
     // particleGun->GeneratePrimaryVertex(anEvent);


     //std::cout<<"************** CHECK: EventID = "<<anEvent->GetEventID()<<std::endl;
     if(anEvent->GetEventID() > nrCascadeTree->tree->GetEntries()) {
      G4cout<<"Geant4 eventID exceeds input tree eventID. Can not generate more events. Exiting."<<G4endl;
      exit(0);
     }
     nrCascadeTree->tree->GetEntry(anEvent->GetEventID());
     // std::cout<<"Number of particles to be produced in this event = "<<nrCascadeTree->nE_gamma<<std::endl;
     // for(int i = 0; i < nrCascadeTree->nE_gamma; i++){
     //  std::cout<<" \t i : energy :: "<<i<<" : "<<nrCascadeTree->E_gamma[i]<<std::endl;
     // }


     if(nrCascadeTree->nE_gamma == 0) {
      particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("geantino"));
      std::vector<G4double> angles = GenerateRandomDirection();
      particleGun->SetParticlePosition(G4ThreeVector(0.0*m,0.0*m,0.0*m));
      particleGun->SetParticleMomentumDirection(G4ThreeVector(angles[2]*cos(angles[0]),angles[2]*sin(angles[0]),angles[1]));
      particleGun->SetParticleEnergy(1.e-10);
      particleGun->GeneratePrimaryVertex(anEvent);
      //std::cout<<"XXXXXXX Geantino....."<<std::endl;
      
     } 
     else {
      for(int i = 0; i < nrCascadeTree->nE_gamma; i++) {
        particleGun->SetParticleDefinition(G4Gamma::Definition());
        std::vector<G4double> angles = GenerateRandomDirection();
        //particleGun->SetParticlePosition(G4ThreeVector(-3.829808*mm, -46.01114*mm, 22.045949*mm));
        particleGun->SetParticlePosition(G4ThreeVector(nrCascadeTree->x_cap *mm, nrCascadeTree->y_cap *mm, nrCascadeTree->z_cap *mm));
        particleGun->SetParticleMomentumDirection(G4ThreeVector(angles[2]*cos(angles[0]),angles[2]*sin(angles[0]),angles[1]));
        //particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,0.));
        particleGun->SetParticleEnergy(nrCascadeTree->E_gamma[i]);
        //std::cout<<" Event : n : i : Eg : x : y : z :: "<<anEvent->GetEventID()<<" : "<<nrCascadeTree->nE_gamma<<" : "<<i<<" : "<<nrCascadeTree->E_gamma[i]<<" : "<<nrCascadeTree->x_cap<<" : "<<nrCascadeTree->y_cap<<" : "<<nrCascadeTree->z_cap<<std::endl;
        particleGun->GeneratePrimaryVertex(anEvent);

      }
      //std::cout<<"XXXXXXXXXXX  Gamma....."<<std::endl;
     }

     //std::cout<<"XXXXXXXXXXX  HELLOOOOOOOOOOOOOOOOOO....."<<std::endl;

     //particleGun->GeneratePrimaryVertex(anEvent);
     //if(!nrCascadeTree) std::cout<<"File not"<<anEvent->GetEventID()<<std::endl;
   }
   else
   {
     particleSource->GeneratePrimaryVertex(anEvent);
     //found out particleSource returns the energy it JUST threw, not the one it will throw next by looking here:
     //http://www.apc.univ-paris7.fr/~franco/g4doxy/html/classG4GeneralParticleSource.html#cfe1960793b123fe1b50facd722e76d5
     G4double energy = particleSource->GetParticleEnergy();

     //do some particle-gun throwing
     particleGun->SetParticleDefinition(G4Gamma::Definition()); 
     particleGun->SetParticlePosition(particleSource->GetParticlePosition());
     G4int level=0;
     if(energy>9.0){
       level = 0; //ground state
     }
     else if(energy>6.0 && energy<=9.0){
       level = 1;
     }
     else if(energy>1.75 && energy<=6.0){
       level = 2;
     }
     else if(energy<=1.75){
       level = 3;
     }

     //check if the gammas are on/off
     if(!Ndet->GetConstructPuBeSourceAndShield_doPuBeGamma()) level=5;

     //G4cout << "Level: " << level << G4endl; 
     //G4cout << "N Energy: " << energy << G4endl; 
     std::vector<G4double> angles = GenerateRandomDirection();
     switch(level){

       case 0:
         //pure event, do nothing
         //G4cout << "Ground State!" << G4endl;
	 break;
	
       case 1:
         //G4cout << "First State!" << G4endl;
         particleGun->SetParticleEnergy(4.43803);
         angles = GenerateRandomDirection();
         particleGun->SetParticleMomentumDirection(G4ThreeVector(angles[2]*cos(angles[0]),angles[2]*sin(angles[0]),angles[1]));
         particleGun->GeneratePrimaryVertex(anEvent);
	 break;

       case 2:
         //G4cout << "Second State!" << G4endl;
         particleGun->SetParticleEnergy(4.43803);
         angles = GenerateRandomDirection();
         particleGun->SetParticleMomentumDirection(G4ThreeVector(angles[2]*cos(angles[0]),angles[2]*sin(angles[0]),angles[1]));
         particleGun->GeneratePrimaryVertex(anEvent);
         particleGun->SetParticleEnergy(3.21483);
         angles = GenerateRandomDirection();
         particleGun->SetParticleMomentumDirection(G4ThreeVector(angles[2]*cos(angles[0]),angles[2]*sin(angles[0]),angles[1]));
         particleGun->GeneratePrimaryVertex(anEvent);
	 break;

       case 3:
         //G4cout << "Third State!" << G4endl;
         particleGun->SetParticleEnergy(9.637);
         angles = GenerateRandomDirection();
         particleGun->SetParticleMomentumDirection(G4ThreeVector(angles[2]*cos(angles[0]),angles[2]*sin(angles[0]),angles[1]));
         particleGun->GeneratePrimaryVertex(anEvent);
	 break;

       default:
         //G4cout << "No Level Info!" << G4endl;
	 break;

     }


   }
  
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
