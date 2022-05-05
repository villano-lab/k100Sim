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
#include "k100_nrCascadeTree.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

//ROOT stuff
#include "TH1D.h"

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

  k100_PrimaryGeneratorAction(G4bool useCapture=false, G4String infile="none");
  ~k100_PrimaryGeneratorAction();

public:

  void GeneratePrimaries(G4Event* anEvent);

private:

  k100_DetectorConstruction *Ndet; 
  k100_nrCascadeTree* nrCascadeTree;
  G4ParticleGun* particleGun;
  G4GeneralParticleSource*      particleSource;	  
  G4RotationMatrix *xrot;
  //G4GeneralParticleSource* particleGun;

  //helpful functions
  std::vector<G4double> GenerateRandomDirection();

  //useful variables
  G4bool                        throwCaptures; // false for GeneralParticleSource
  G4String InFile;

};

#endif
