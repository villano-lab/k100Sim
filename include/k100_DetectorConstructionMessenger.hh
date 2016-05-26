// ------------------------------------------------
//
// k100_DetectorconstructionMessenger.hh : 2016 
//
// ------------------------------------------------

#ifndef k100_DetectorMessenger_h
#define k100_DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;

class k100_DetectorConstruction;

// ------------------------------------------------

class k100_DetectorConstructionMessenger: public G4UImessenger
{
public:
  k100_DetectorConstructionMessenger(k100_DetectorConstruction* );
  ~k100_DetectorConstructionMessenger();

  void SetNewValue(G4UIcommand*, G4String);

private:
  k100_DetectorConstruction* k100_Detector;
  G4UIdirectory*          k100_detDir;

  // Detector Element Activation

  G4UIcmdWithAString*     DetectorActivateCmd;
  G4UIcmdWithAString*     DetectorDeActivateCmd;
  G4UIcmdWithAnInteger*   NbTowersCmd;

  G4UIcmdWithAString*     DrawSolidBox;
  G4UIcmdWithAString*     DrawFrameBox;


};

#endif
