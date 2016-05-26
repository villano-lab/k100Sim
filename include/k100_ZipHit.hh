// ------------------------------------------------
//
// k100_ZipHit.hh : Dec 2005
//
// ------------------------------------------------

#ifndef k100_ZipHit_h
#define k100_ZipHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

// ------------------------------------------

class k100_ZipHit : public G4VHit
{
public:
  k100_ZipHit();
  k100_ZipHit(const k100_ZipHit&);
	
  ~k100_ZipHit();
	
  const k100_ZipHit& operator=(const k100_ZipHit&);
  int operator==(const k100_ZipHit&) const;
	
  inline void* operator new(size_t);
  inline void  operator delete(void*);
	
  void Draw();
  void Print();
	
public:
  void SetTrackID  (G4int track)       {trackID = track;};
  void SetZipNb    (G4int zip)         {zipNb = zip;};
  void SetEdep     (G4double de)       {edep = de;};
  void SetGlobalTime(G4double gTime) {globalTime = gTime;};
  void SetPos      (G4ThreeVector xyz) {pos = xyz;};
  void SetPName    (G4String pName)    {particleName = pName;};
	
  G4int GetTrackID() {return trackID;};
  G4int GetZipNb() {return zipNb;};
  G4double GetEdep() {return edep;};
  G4double GetGlobalTime() {return globalTime;};
  G4ThreeVector GetPos() {return pos;};
  G4String GetParticle() {return particleName;};
  G4double* GetData() {return dataVector;};

private:
  G4int         trackID;
  G4int         zipNb;
  G4double      edep;
  G4double      globalTime;
  G4ThreeVector pos;
  G4String      particleName;
  G4double*     dataVector;
};

// ------------------------------------------

typedef G4THitsCollection<k100_ZipHit> k100_ZipHitsCollection;

extern G4Allocator<k100_ZipHit> k100_ZipHitAllocator;

// ------------------------------------------

inline void* k100_ZipHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) k100_ZipHitAllocator.MallocSingle();
  return aHit;
}

// ------------------------------------------

inline void k100_ZipHit::operator delete(void *aHit)
{
  k100_ZipHitAllocator.FreeSingle((k100_ZipHit*) aHit);
}

// ------------------------------------------

#endif

