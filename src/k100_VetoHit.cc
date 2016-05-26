// ------------------------------------------------
//
// k100_VetoHit.cc : Dec 2005
//
// ------------------------------------------------

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "k100_VetoHit.hh"

G4Allocator<k100_VetoHit> k100_VetoHitAllocator;

// ------------------------------------------

k100_VetoHit::k100_VetoHit() 
{
  dataVector = new G4double[13];
}

// ------------------------------------------

k100_VetoHit::k100_VetoHit(const k100_VetoHit& right) : G4VHit()
{
  trackID    = right.trackID;
  edep       = right.edep;
  globalTime = right.globalTime;
  pos        = right.pos;
  dataVector = right.dataVector;
}

// ------------------------------------------

k100_VetoHit::~k100_VetoHit() 
{
  delete dataVector;
}

// ------------------------------------------

const k100_VetoHit& k100_VetoHit::operator=(const k100_VetoHit& right)
{
  trackID    = right.trackID;
  edep       = right.edep;
  globalTime = right.globalTime;
  pos        = right.pos;
  dataVector = right.dataVector;
  return *this;
}

// ------------------------------------------

int k100_VetoHit::operator==(const k100_VetoHit& /*right*/) const
{
  return 0;
}

// ------------------------------------------

void k100_VetoHit::Draw()
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

void k100_VetoHit::Print()
{
  G4cout << "  trackID: " << trackID << "  Veto" 
	 << "  energy deposit: " << G4BestUnit(edep,"Energy")
         << "  global time: " << G4BestUnit(globalTime,"Time")
	 << "  position: " << G4BestUnit(pos,"Length") << G4endl;
}

// ------------------------------------------
