/*==================k100_EventInfo.hh========================= 
   
      PROGRAMMER:  Anthony Villano  04/10/15

      UPDATES:      
       

      PURPOSE: Class for a k100 to store some add'l event info
              
======================================================================*/

#ifndef k100_EventInfo_h
#define k100_EventInfo_h 1

// ------------------------------------------------

#include "globals.hh"


// ------------------------------------------------

class k100_EventInfo : public G4VUserEventInformation
{
public :
  k100_EventInfo();
  ~k100_EventInfo();

  //pure virtuals
  void Print() const;

  G4double GetBeLength(){return BeLength;}
  void SetBeLength(G4double length){BeLength = length;}

public :
  
  
private :

  G4double BeLength;

};

#endif
