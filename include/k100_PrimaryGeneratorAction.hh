/*==================k100_PrimaryGeneratorAction.hh================= 
   
      PROGRAMMER:  Anthony Villano 10/18/12

      UPDATES:      
       

      PURPOSE: Class which defines the generators for the k100
               simulation.
              
======================================================================*/



#ifndef k100_PrimaryGeneratorAction_h
#define k100_PrimaryGeneratorAction_h 1

#include "G4IonTable.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4VUserPrimaryGeneratorAction.hh"

#include <vector>

// ------------------------------------------------

class G4ParticleGun;
class G4GeneralParticleSource;
class G4Event;
class k100_DetectorConstruction;

// ------------------------------------------------
class k100_PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{

public:

  k100_PrimaryGeneratorAction();
  ~k100_PrimaryGeneratorAction();

public:

  void GeneratePrimaries(G4Event* anEvent);

private:

  k100_DetectorConstruction *Ndet; 

  G4ParticleGun* particleGun;
  G4RotationMatrix *xrot;
  //G4GeneralParticleSource* particleGun;

  //helpful functions
  std::vector<G4double> GenerateRandomDirection();

};

#endif
