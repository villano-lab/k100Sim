// ------------------------------------------------
//
// k100_VetoSD.cc : 2016 
//
// ------------------------------------------------

#ifndef k100_VetoSD_h
#define k100_VetoSD_h 1

#include "G4VSensitiveDetector.hh"
#include "k100_VetoHit.hh"

class G4Step;
class G4HCofThisEvent;

// ------------------------------------------------

class k100_VetoSD : public G4VSensitiveDetector
{
public:
  k100_VetoSD(G4String, G4int panelNb);
  ~k100_VetoSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  k100_VetoHitsCollection* vetoCollection;
  G4int panelNb;

};

// ------------------------------------------------

#endif
