// ------------------------------------------------
//
// k100_ZipHit.cc : Dec 2005
//
// ------------------------------------------------

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "k100_DataStorage.hh" // This is used to get the value of N_DATA
#include "k100_ZipHit.hh"

G4Allocator<k100_ZipHit> k100_ZipHitAllocator;

// ------------------------------------------

k100_ZipHit::k100_ZipHit() 
{
  dataVector = new G4double[N_DATA];
}

// ------------------------------------------

k100_ZipHit::k100_ZipHit(const k100_ZipHit& right) : G4VHit()
{
  trackID    = right.trackID;
  edep       = right.edep;
  globalTime = right.globalTime;
  pos        = right.pos;
  dataVector = right.dataVector;
}

// ------------------------------------------

k100_ZipHit::~k100_ZipHit() 
{
  delete dataVector;
}

// ------------------------------------------

const k100_ZipHit& k100_ZipHit::operator=(const k100_ZipHit& right)
{
  trackID    = right.trackID;
  edep       = right.edep;
  globalTime = right.globalTime;
  pos        = right.pos;
  dataVector = right.dataVector;
  return *this;
}

// ------------------------------------------

int k100_ZipHit::operator==(const k100_ZipHit& /*right*/) const
{
  return 0;
}

// ------------------------------------------

void k100_ZipHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
    {
      G4Circle circle(pos);
      circle.SetScreenSize(10.04);
      circle.SetFillStyle(G4Circle::filled);
      G4Colour colour(0/255.,100/255.,128/255.);
      G4VisAttributes attribs(colour);
      circle.SetVisAttributes(attribs);
      pVVisManager->Draw(circle);
    }
}

// ------------------------------------------

void k100_ZipHit::Print()
{
  G4cout << "  trackID: " << trackID << "  Zip" 
	 << "  energy deposit: " << G4BestUnit(edep,"Energy")
         << "  global time: " << G4BestUnit(globalTime, "Time")
	 << "  position: " << G4BestUnit(pos,"Length") << G4endl;
}

// ------------------------------------------

