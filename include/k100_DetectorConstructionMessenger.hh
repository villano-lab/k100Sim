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
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADouble;

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
  G4UIcmdWith3VectorAndUnit*   DetSizeCmd;
  G4UIcmdWithAString*     GPSShieldMatCmd;
  G4UIcmdWith3VectorAndUnit*   GeGammaCoinPositionCmd;
  G4UIcmdWith3VectorAndUnit*   GeGammaCoinSizeCmd;
  G4UIcmdWithAString*     SimpleGammaCoinMatCmd;
  G4UIcmdWithABool*       ZipConfigureCmd_Mat1;
  G4UIcmdWithABool*       ShieldConfigureCmd_HPGeboron_shield;
  G4UIcmdWithABool*       ShieldConfigureCmd_HPGeboron;
  G4UIcmdWithABool*       ShieldConfigureCmd_SouthNaI;
  G4UIcmdWithABool*       ShieldConfigureCmd_BasePoly;
  G4UIcmdWithABool*       ShieldConfigureCmd_BaseLead;
  G4UIcmdWithAnInteger*       ShieldConfigureCmd_mod;
  G4UIcmdWithABool*       FridgeConfigureCmd_pure3HeBath;
  G4UIcmdWithABool*       PuBeConfigureCmd_westPolySensitive;
  G4UIcmdWithABool*       PuBeConfigureCmd_doPuBeGammas;
  G4UIcmdWithABool*       PuBeConfigureCmd_Barrel;
  G4UIcmdWithABool*       PuBeConfigureCmd_NaIsensitive;
  G4UIcmdWithABool*       PuBeConfigureCmd_R66;
  G4UIcmdWithABool*       PuBeConfigureCmd_R62;
  G4UIcmdWithAnInteger*       PuBeConfigureCmd_mod;
  G4UIcmdWithABool*       PuBeConfigureCmd_Orb;
  G4UIcmdWith3VectorAndUnit*   PuBeConfigureCmd_OrbPos;
  G4UIcmdWithADoubleAndUnit*   PuBeConfigureCmd_OrbRad;

  G4UIcmdWithAString*     DrawSolidBox;
  G4UIcmdWithAString*     DrawFrameBox;

  G4UIcmdWithADouble* SodiumBorateDensityFractionCmd;
  G4UIcmdWithADouble* BoronShieldThicknessCmd;
  G4UIcmdWithAnInteger*   NbBoronVertCmd;
  G4UIcmdWithAnInteger*   NbBoronHoriCmd;


};

#endif
