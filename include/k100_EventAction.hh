/*==================k100_EventAction.hh============================= 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE:  Class to support the track actions carried out
                when an "event" occurs in the simulation.
              
======================================================================*/

#ifndef k100_EventAction_h
#define k100_EventAction_h 1

#include "G4UserEventAction.hh"

class k100_RunAction;
class k100_DataStorage;
class k100_DetectorConstruction;

// ------------------------------------------

class k100_EventAction : public G4UserEventAction
{
public :

  k100_EventAction(G4bool nCapSaveOnly);
  k100_EventAction(k100_RunAction*,G4bool nCapSaveOnly);
  ~k100_EventAction();

public :

  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);

private :


  G4bool saveOnlyNCapture;
  G4bool drawEvent, saveEvent;
  
  //data output class. not used for initial encarnation
  k100_DataStorage* dataOut;
  k100_DetectorConstruction* k100_Detector;
  k100_RunAction* pRunAction;

};

#endif
