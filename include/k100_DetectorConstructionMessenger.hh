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
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithoutParameter;

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

  G4UIcmdWithoutParameter* UpdateGeometryCmd;
  G4UIcmdWithAString*     DetectorActivateCmd;
  G4UIcmdWithAString*     DetectorDeActivateCmd;
  G4UIcmdWithAnInteger*   NbTowersCmd;
  G4UIcmdWith3VectorAndUnit*   GPSShieldPositionCmd;
  G4UIcmdWith3VectorAndUnit*   GPSShieldSizeCmd;
  G4UIcmdWithAString*     GPSShieldMatCmd;
  G4UIcmdWith3VectorAndUnit*   GeGammaCoinPositionCmd;
  G4UIcmdWith3VectorAndUnit*   GeGammaCoinSizeCmd;
  G4UIcmdWithAString*     SimpleGammaCoinMatCmd;

  G4UIcmdWithAString*     DrawSolidBox;
  G4UIcmdWithAString*     DrawFrameBox;


};

#endif
