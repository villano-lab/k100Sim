/*==================k100_RunActionMessenger.hh===================== 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE: Class to supply an interactive interface to the RunAction
               class.  Examples of tasks include setting the file name
	       prefix for the output. This class will supply commands
	       accessible in standard Geant4 Macros
      
              
======================================================================*/

#ifndef k100_RunActionMessenger_h
#define k100_RunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

// ------------------------------------------------

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;

class k100_RunAction;

// ------------------------------------------------

class k100_RunActionMessenger: public G4UImessenger
{
public:
  k100_RunActionMessenger(k100_RunAction* );
  ~k100_RunActionMessenger();

  void SetNewValue(G4UIcommand*, G4String);

private:
  k100_RunAction*       pRunAction;
  G4UIdirectory*        k100_RunDir;

  // Run Element Activation

  G4UIcmdWithABool*       setDrawEventCmd;
  G4UIcmdWithABool*       setNCapOutputCmd;
  G4UIcmdWithAString*       setRunFileName;
};

#endif
