/*==================k100_EventInfo.cc============================ 
   
      PROGRAMMER:  Anthony Villano 04/10/15 

      UPDATES:      
       

      PURPOSE: Code support for the ascii output in the k100 simulation.
              
======================================================================*/

#include "G4RunManager.hh" 
#include "G4UnitsTable.hh"

#include "k100_EventInfo.hh"


k100_EventInfo::k100_EventInfo():
BeLength(-1.0)
{

}
k100_EventInfo::~k100_EventInfo()
{
  
}
void k100_EventInfo::Print() const
{
  G4cout << "The gamma path length through the Be: " << BeLength << G4endl;

  return;
}
// ------------------------------------------------

