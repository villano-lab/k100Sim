// ------------------------------------------------
//
// k100_ZipSD.cc : Dec 2005
//
// ------------------------------------------------

#ifndef k100_ZipSD_h
#define k100_ZipSD_h 1

#include "G4VSensitiveDetector.hh"
#include "k100_ZipHit.hh"

class G4Step;
class G4HCofThisEvent;

// ------------------------------------------------

class k100_ZipSD : public G4VSensitiveDetector
{
public:
  k100_ZipSD(G4String, G4int towerNb);
  ~k100_ZipSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  k100_ZipHitsCollection* zipCollection;
  G4int towNb;

  G4bool isNCap(G4Track* track,G4Step *aStep);
};

// ------------------------------------------------

#endif
