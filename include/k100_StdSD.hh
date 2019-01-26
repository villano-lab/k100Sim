/*==================k100_StdSD.hh=================================== 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE: Class for specifying a sensitive detector in the k100 
               simulation.
              
======================================================================*/

#ifndef k100_StdSD_h
#define k100_StdSD_h 1

#include "G4VSensitiveDetector.hh"
#include "k100_StdHit.hh"

class G4Step;
class G4HCofThisEvent;
class k100_PrimaryGeneratorAction;

class k100_StdSD : public G4VSensitiveDetector
{
public:
  k100_StdSD(G4String, G4int panelNb);
  ~k100_StdSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  k100_StdHitsCollection* stdCollection;
  G4int genericNb;

  G4bool isNCap(G4Track* track,G4Step *aStep);
  k100_PrimaryGeneratorAction *Ngen;

};

#endif
