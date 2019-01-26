/*==================k100_StdHit.cc================================== 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE: Code support for the k100_StdHit object and methods. 
              
======================================================================*/

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "k100_DataStorage.hh" // This is used to get the value of N_DATA
#include "k100_StdHit.hh"

G4Allocator<k100_StdHit> k100_StdHitAllocator;

k100_StdHit::k100_StdHit() 
{
  dataVector = new G4double[N_DATA];
}
k100_StdHit::k100_StdHit(const k100_StdHit& right) : G4VHit()
{
  trackID    = right.trackID;
  edep       = right.edep;
  globalTime = right.globalTime;
  pos        = right.pos;
  dataVector = right.dataVector;
}
k100_StdHit::~k100_StdHit() 
{
  delete dataVector;
}
const k100_StdHit& k100_StdHit::operator=(const k100_StdHit& right)
{
  trackID    = right.trackID;
  edep       = right.edep;
  globalTime = right.globalTime;
  pos        = right.pos;
  dataVector = right.dataVector;
  return *this;
}
int k100_StdHit::operator==(const k100_StdHit& right) const
{
  return 0;
}
void k100_StdHit::Draw()
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
void k100_StdHit::Print()
{
  G4cout << "  trackID: " << trackID << "  Standard Sensitive" 
	 << "  energy deposit: " << G4BestUnit(edep,"Energy")
	 << "  position: " << G4BestUnit(pos,"Length") << G4endl;
}

// ------------------------------------------

