// ------------------------------------------------
//
// k100_DetectorConstructionMessenger.cc : 2016 
//
// ------------------------------------------------

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"

#include "k100_DetectorConstructionMessenger.hh"
#include "k100_DetectorConstruction.hh"

// ------------------------------------------------

k100_DetectorConstructionMessenger::k100_DetectorConstructionMessenger(k100_DetectorConstruction * k100_Det)
  :k100_Detector(k100_Det)
{ 

  // Directory created under RunActionMessenger <--- NOT YET !!!
  k100_detDir = new G4UIdirectory("/CDMS/");
  k100_detDir->SetGuidance("CDMS specific controls.");
  //
  k100_detDir = new G4UIdirectory("/CDMS/detector/");
  k100_detDir->SetGuidance("CDMS Detector Control.");
  //
  k100_detDir = new G4UIdirectory("/CDMS/rendering/");
  k100_detDir->SetGuidance("CDMS Detector Rendering.");
  //
  k100_detDir = new G4UIdirectory("/CDMS/genericShield/");
  k100_detDir->SetGuidance("CDMS k100 generic shield construction.");
  //
  k100_detDir = new G4UIdirectory("/CDMS/gammaCoin/");
  k100_detDir->SetGuidance("CDMS k100 outside-cryostat gamma coincidence detector.");

  // Turn various components on/off


  UpdateGeometryCmd = new G4UIcmdWithoutParameter("/CDMS/updateGeometry",this);
  UpdateGeometryCmd->SetGuidance("Rebuild the geometry.");

  DetectorActivateCmd = new G4UIcmdWithAString("/CDMS/detector/activate",this);
  DetectorActivateCmd->SetGuidance("Activate CDMS Detector Element.");
  DetectorActivateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  DetectorActivateCmd->SetGuidance("in order for change to take effect.");
  DetectorActivateCmd->SetGuidance("Choices are : Zips/Towers/Veto/Shields/IceBox/NaIArray/BoronShield/PolyBox .");
  DetectorActivateCmd->SetParameterName("choice",false);
  DetectorActivateCmd->AvailableForStates(G4State_Idle);

  DetectorDeActivateCmd = new G4UIcmdWithAString("/CDMS/detector/deactivate",this);
  DetectorDeActivateCmd->SetGuidance("Deactivate CDMS Detector Element.");
  DetectorDeActivateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  DetectorDeActivateCmd->SetGuidance("in order for change to take effect.");
  DetectorDeActivateCmd->SetGuidance("Choices are : Zip/Tower/Veto/Shields/IceBox/NaIArray/BoronShield/polyBox .");
  DetectorDeActivateCmd->SetParameterName("choice",false);
  DetectorDeActivateCmd->AvailableForStates(G4State_Idle);

  ZipConfigureCmd_Mat1 = new G4UIcmdWithABool("/CDMS/Zip1/IsGe",this);
  ZipConfigureCmd_Mat1->SetGuidance("Toggle Ge/Si Material for Zip1.");
  ZipConfigureCmd_Mat1->SetGuidance("This command MUST be applied before \"beamOn\" ");
  ZipConfigureCmd_Mat1->SetGuidance("in order for change to take effect.");
  ZipConfigureCmd_Mat1->SetParameterName("choice",false);
  ZipConfigureCmd_Mat1->AvailableForStates(G4State_Idle);

  ShieldConfigureCmd_HPGeboron_shield = new G4UIcmdWithABool("/CDMS/Shield/HPGeboron_shield",this);
  ShieldConfigureCmd_HPGeboron_shield->SetGuidance("Toggle HPGe boron Shield.");
  ShieldConfigureCmd_HPGeboron_shield->SetGuidance("This command MUST be applied before \"beamOn\" ");
  ShieldConfigureCmd_HPGeboron_shield->SetGuidance("in order for change to take effect.");
  ShieldConfigureCmd_HPGeboron_shield->SetParameterName("choice",false);
  ShieldConfigureCmd_HPGeboron_shield->AvailableForStates(G4State_Idle);

  ShieldConfigureCmd_HPGeboron = new G4UIcmdWithABool("/CDMS/Shield/HPGeboron",this);
  ShieldConfigureCmd_HPGeboron->SetGuidance("Toggle HPGe construction on Shield.");
  ShieldConfigureCmd_HPGeboron->SetGuidance("This command MUST be applied before \"beamOn\" ");
  ShieldConfigureCmd_HPGeboron->SetGuidance("in order for change to take effect.");
  ShieldConfigureCmd_HPGeboron->SetParameterName("choice",false);
  ShieldConfigureCmd_HPGeboron->AvailableForStates(G4State_Idle);

  ShieldConfigureCmd_SouthNaI = new G4UIcmdWithABool("/CDMS/Shield/SouthNaI",this);
  ShieldConfigureCmd_SouthNaI->SetGuidance("Toggle NaI construction on Shield.");
  ShieldConfigureCmd_SouthNaI->SetGuidance("This command MUST be applied before \"beamOn\" ");
  ShieldConfigureCmd_SouthNaI->SetGuidance("in order for change to take effect.");
  ShieldConfigureCmd_SouthNaI->SetParameterName("choice",false);
  ShieldConfigureCmd_SouthNaI->AvailableForStates(G4State_Idle);

  ShieldConfigureCmd_BasePoly = new G4UIcmdWithABool("/CDMS/Shield/BasePoly",this);
  ShieldConfigureCmd_BasePoly->SetGuidance("Toggle base poly panel construction on Shield.");
  ShieldConfigureCmd_BasePoly->SetGuidance("This command MUST be applied before \"beamOn\" ");
  ShieldConfigureCmd_BasePoly->SetGuidance("in order for change to take effect.");
  ShieldConfigureCmd_BasePoly->SetParameterName("choice",false);
  ShieldConfigureCmd_BasePoly->AvailableForStates(G4State_Idle);

  ShieldConfigureCmd_BaseLead = new G4UIcmdWithABool("/CDMS/Shield/BaseLead",this);
  ShieldConfigureCmd_BaseLead->SetGuidance("Toggle base lead layer construction on Shield.");
  ShieldConfigureCmd_BaseLead->SetGuidance("This command MUST be applied before \"beamOn\" ");
  ShieldConfigureCmd_BaseLead->SetGuidance("in order for change to take effect.");
  ShieldConfigureCmd_BaseLead->SetParameterName("choice",false);
  ShieldConfigureCmd_BaseLead->AvailableForStates(G4State_Idle);

  ShieldConfigureCmd_mod = new G4UIcmdWithAnInteger("/CDMS/Shield/Mod",this);
  ShieldConfigureCmd_mod->SetGuidance("Set integer corresponding to shield modification.");
  ShieldConfigureCmd_mod->SetParameterName("Mod",false);
  ShieldConfigureCmd_mod->SetRange("Mod==0 || Mod==1 || Mod==2");
  ShieldConfigureCmd_mod->AvailableForStates(G4State_Idle);

  FridgeConfigureCmd_pure3HeBath = new G4UIcmdWithABool("/CDMS/Fridge/pure3HeBath",this);
  FridgeConfigureCmd_pure3HeBath->SetGuidance("Toggle to 3He in bath.");
  FridgeConfigureCmd_pure3HeBath->SetGuidance("This command MUST be applied before \"beamOn\" ");
  FridgeConfigureCmd_pure3HeBath->SetGuidance("in order for change to take effect.");
  FridgeConfigureCmd_pure3HeBath->SetParameterName("choice",false);
  FridgeConfigureCmd_pure3HeBath->AvailableForStates(G4State_Idle);

  PuBeConfigureCmd_westPolySensitive = new G4UIcmdWithABool("/CDMS/PuBe/westPolySensitive",this);
  PuBeConfigureCmd_westPolySensitive->SetGuidance("Toggle west poly wall sensitivity.");
  PuBeConfigureCmd_westPolySensitive->SetGuidance("This command MUST be applied before \"beamOn\" ");
  PuBeConfigureCmd_westPolySensitive->SetGuidance("in order for change to take effect.");
  PuBeConfigureCmd_westPolySensitive->SetParameterName("choice",false);
  PuBeConfigureCmd_westPolySensitive->AvailableForStates(G4State_Idle);

  PuBeConfigureCmd_doPuBeGammas = new G4UIcmdWithABool("/CDMS/PuBe/doPuBeGammas",this);
  PuBeConfigureCmd_doPuBeGammas->SetGuidance("Do or do not do gammas.");
  PuBeConfigureCmd_doPuBeGammas->SetGuidance("This command MUST be applied before \"beamOn\" ");
  PuBeConfigureCmd_doPuBeGammas->SetGuidance("in order for change to take effect.");
  PuBeConfigureCmd_doPuBeGammas->SetParameterName("choice",false);
  PuBeConfigureCmd_doPuBeGammas->AvailableForStates(G4State_Idle);

  PuBeConfigureCmd_Barrel = new G4UIcmdWithABool("/CDMS/PuBe/Barrel",this);
  PuBeConfigureCmd_Barrel->SetGuidance("Toggle Barrel addition to source.");
  PuBeConfigureCmd_Barrel->SetGuidance("This command MUST be applied before \"beamOn\" ");
  PuBeConfigureCmd_Barrel->SetGuidance("in order for change to take effect.");
  PuBeConfigureCmd_Barrel->SetParameterName("choice",false);
  PuBeConfigureCmd_Barrel->AvailableForStates(G4State_Idle);

  PuBeConfigureCmd_NaIsensitive = new G4UIcmdWithABool("/CDMS/PuBe/NaIsensitive",this);
  PuBeConfigureCmd_NaIsensitive->SetGuidance("Toggle NaI sensitivity.");
  PuBeConfigureCmd_NaIsensitive->SetGuidance("This command MUST be applied before \"beamOn\" ");
  PuBeConfigureCmd_NaIsensitive->SetGuidance("in order for change to take effect.");
  PuBeConfigureCmd_NaIsensitive->SetParameterName("choice",false);
  PuBeConfigureCmd_NaIsensitive->AvailableForStates(G4State_Idle);

  PuBeConfigureCmd_R66 = new G4UIcmdWithABool("/CDMS/PuBe/R66",this);
  PuBeConfigureCmd_R66->SetGuidance("Toggle R66 source construction conditions.");
  PuBeConfigureCmd_R66->SetGuidance("This command MUST be applied before \"beamOn\" ");
  PuBeConfigureCmd_R66->SetGuidance("in order for change to take effect.");
  PuBeConfigureCmd_R66->SetParameterName("choice",false);
  PuBeConfigureCmd_R66->AvailableForStates(G4State_Idle);

  PuBeConfigureCmd_R62 = new G4UIcmdWithABool("/CDMS/PuBe/R62",this);
  PuBeConfigureCmd_R62->SetGuidance("Toggle R62 source construction conditions.");
  PuBeConfigureCmd_R62->SetGuidance("This command MUST be applied before \"beamOn\" ");
  PuBeConfigureCmd_R62->SetGuidance("in order for change to take effect.");
  PuBeConfigureCmd_R62->SetParameterName("choice",false);
  PuBeConfigureCmd_R62->AvailableForStates(G4State_Idle);

  PuBeConfigureCmd_mod = new G4UIcmdWithAnInteger("/CDMS/PuBe/Mod",this);
  PuBeConfigureCmd_mod->SetGuidance("Set integer corresponding to shield modification.");
  PuBeConfigureCmd_mod->SetParameterName("Mod",false);
  PuBeConfigureCmd_mod->SetRange("Mod==0 || Mod==1 || Mod==2");
  PuBeConfigureCmd_mod->AvailableForStates(G4State_Idle);

  PuBeConfigureCmd_Orb = new G4UIcmdWithABool("/CDMS/PuBe/Orb",this);
  PuBeConfigureCmd_Orb->SetGuidance("Do simple lead orb shielding.");
  PuBeConfigureCmd_Orb->SetGuidance("This command MUST be applied before \"beamOn\" ");
  PuBeConfigureCmd_Orb->SetGuidance("in order for change to take effect.");
  PuBeConfigureCmd_Orb->SetParameterName("choice",false);
  PuBeConfigureCmd_Orb->AvailableForStates(G4State_Idle);

  PuBeConfigureCmd_OrbPos = new G4UIcmdWith3VectorAndUnit("/CDMS/PuBe/setOrbPos",this);
  PuBeConfigureCmd_OrbPos->SetGuidance("Set the shielding Orb position.");
  PuBeConfigureCmd_OrbPos->SetParameterName("X","Y","Z",true,true);
  //PuBeConfigureCmd_OrbPos->SetRange("X != 0 && Y != 0 && Z != 0"); //none of these can be zero
  
  PuBeConfigureCmd_OrbRad = new G4UIcmdWithADoubleAndUnit("/CDMS/PuBe/setOrbRadius",this);
  PuBeConfigureCmd_OrbRad->SetGuidance("Set the shielding Orb radius.");
  PuBeConfigureCmd_OrbRad->SetParameterName("R",true,true);
  PuBeConfigureCmd_OrbRad->SetRange("R != 0"); //none of these can be zero
 
  //set parameters for generic shielding for quick sims
  GPSShieldPositionCmd = new G4UIcmdWith3VectorAndUnit("/CDMS/genericShield/setPosition",this);
  GPSShieldPositionCmd->SetGuidance("Set position for generic rectangular shielding.");
  GPSShieldPositionCmd->SetParameterName("Dx","Dy","Dz",true,true);
  GPSShieldPositionCmd->SetRange("Dx != 0 || Dy != 0 || Dz != 0"); //FIXME actually might want to put more constraints? no overlap?

  GPSShieldSizeCmd = new G4UIcmdWith3VectorAndUnit("/CDMS/genericShield/setSize",this);
  GPSShieldSizeCmd->SetGuidance("Set size for generic rectangular shielding.");
  GPSShieldSizeCmd->SetParameterName("Length","Width","Thickness",true,true);
  GPSShieldSizeCmd->SetRange("Length != 0 && Width != 0 && Thickness != 0"); //none of these can be zero

  GPSShieldMatCmd = new G4UIcmdWithAString("/CDMS/genericShield/setMat",this);
  GPSShieldMatCmd->SetGuidance("Set a material, options: 'Lead', 'Poly'");
  // Select between Solid / WireFrame drawing mode
 
  DetSizeCmd = new G4UIcmdWith3VectorAndUnit("/CDMS/Det/setSize",this);
  DetSizeCmd->SetGuidance("Set size for generic ZIP/HV detector.");
  DetSizeCmd->SetParameterName("Diameter","Thickness","unsed",true,true);
  DetSizeCmd->SetRange("Diameter != 0 && Thickness != 0 "); //none of these can be zero

  // Set parameters for gamma-coincidence detector
  GeGammaCoinPositionCmd = new G4UIcmdWith3VectorAndUnit("/CDMS/gammaCoin/setPosition",this);
  GeGammaCoinPositionCmd->SetGuidance("Set position for HPGe coincidence detector.");
  GeGammaCoinPositionCmd->SetParameterName("Dx","Dy","Dz",true,true);
  GeGammaCoinPositionCmd->SetRange("Dx != 0 || Dy != 0 || Dz != 0"); //FIXME actually might want to put more constraints? no overlap?

  GeGammaCoinSizeCmd = new G4UIcmdWith3VectorAndUnit("/CDMS/gammaCoin/setSize",this);
  GeGammaCoinSizeCmd->SetGuidance("Set size for HPGe coincidence detector.");
  GeGammaCoinSizeCmd->SetParameterName("Radius","Thickness","Nothing",true,true);
  GeGammaCoinSizeCmd->SetRange("Radius != 0 && Thickness != 0"); //none of these can be zero

  SimpleGammaCoinMatCmd = new G4UIcmdWithAString("/CDMS/gammaCoin/setMat",this);
  SimpleGammaCoinMatCmd->SetGuidance("Set a material, options: 'HPGe', 'NaI', 'Scint'");

  DrawSolidBox = new G4UIcmdWithAString("/CDMS/rendering/solid",this);
  DrawSolidBox->SetGuidance("Draw solid box for this detector element.");
  DrawSolidBox->SetGuidance("Choices are : Zips/Towers/Veto/Shields/IceBox .");
  DrawSolidBox->SetParameterName("choise", false);
  DrawSolidBox->AvailableForStates(G4State_Idle);


  DrawFrameBox = new G4UIcmdWithAString("/CDMS/rendering/frame",this);
  DrawFrameBox->SetGuidance("Draw frame box for this detector element.");
  DrawFrameBox->SetGuidance("Choices are : Zips/Towers/Veto/Shields/IceBox .");
  DrawFrameBox->SetParameterName("choise", false);
  DrawFrameBox->AvailableForStates(G4State_Idle);

  // Set number of Zips (for future use)

  NbTowersCmd = new G4UIcmdWithAnInteger("/CDMS/detector/setNbOfTowers",this);
  NbTowersCmd->SetGuidance("Set number of Towers detectors.");
  NbTowersCmd->SetParameterName("NbTowers",false);
  NbTowersCmd->SetRange("NbTowers>0 && NbTowers<6");
  NbTowersCmd->AvailableForStates(G4State_Idle);

  SodiumBorateDensityFractionCmd = new G4UIcmdWithADouble("/CDMS/Shield/SodiumBorateDensityFraction",this);
  SodiumBorateDensityFractionCmd->SetGuidance("Set the density fraction of Sodium Tetraborate decahydrate");
  SodiumBorateDensityFractionCmd->SetParameterName("borateDensityFrac",true);
  SodiumBorateDensityFractionCmd->SetRange("borateDensityFrac > 0 && borateDensityFrac <= 1.0");
  SodiumBorateDensityFractionCmd->AvailableForStates(G4State_Idle);

  BoronShieldThicknessCmd = new G4UIcmdWithADouble("/CDMS/Shield/BoronShieldThickness",this);
  BoronShieldThicknessCmd->SetGuidance("Set the thickness of boron shield surrounding NaIArray");
  BoronShieldThicknessCmd->SetParameterName("boronShieldThicness",true);
  BoronShieldThicknessCmd->SetRange("boronShieldThicness > 0");
  BoronShieldThicknessCmd->AvailableForStates(G4State_Idle);


  NbBoronVertCmd = new G4UIcmdWithAnInteger("/CDMS/Shield/NbBoronShieldVertical",this);
  NbBoronVertCmd->SetGuidance("Set number of vertical boron shields.");
  NbBoronVertCmd->SetParameterName("NbBoronVert",false);
  NbBoronVertCmd->SetRange("NbBoronVert==0 || NbBoronVert==1 || NbBoronVert==2 || NbBoronVert==6");
  NbBoronVertCmd->AvailableForStates(G4State_Idle);

  NbBoronHoriCmd = new G4UIcmdWithAnInteger("/CDMS/Shield/NbBoronShieldHorizontal",this);
  NbBoronHoriCmd->SetGuidance("Set number of horizontal boron shields.");
  NbBoronHoriCmd->SetParameterName("NbBoronHori",false);
  NbBoronHoriCmd->SetRange("NbBoronHori==0 || NbBoronHori==1 || NbBoronHori==2 || NbBoronHori==6");
  NbBoronHoriCmd->AvailableForStates(G4State_Idle);

}

// ------------------------------------------------

k100_DetectorConstructionMessenger::~k100_DetectorConstructionMessenger()
{
  delete DetectorActivateCmd;  delete DetectorDeActivateCmd;
  delete DrawSolidBox; delete DrawFrameBox;
  delete NbTowersCmd; 
}

// ------------------------------------------------

void k100_DetectorConstructionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 

  G4String caseZips  = "Zips";
  G4String caseTower = "Towers";
  G4String caseVeto = "Veto";
  G4String caseShields = "Shields";
  G4String caseIceBox = "IceBox";
  G4String caseNaIArray = "NaIArray"; //spandey
  G4String caseBoronShield = "BoronShield"; //spandey
  G4String casePolyBox = "PolyBox"; //spandey
  G4String caseFrame = "Frame";
  G4String caseFloor = "Floor";
  G4String caseWalls = "Walls";
  G4String caseCeiling = "Ceiling";
  G4String caseWestReflector = "WestReflector";
  G4String casePuBeSourceAndShield = "PuBeSourceAndShield";
  G4String caseThermalNeutronBucket = "ShieldBucket";
  G4String caseGPSShielding = "GPSShielding";
  G4String caseHPGeCoincidence = "HPGeCoincidence";

  if( command == UpdateGeometryCmd ) { 
    k100_Detector->UpdateGeometry();
  }

  if( command == DetectorActivateCmd ) { 
    if(newValue == caseZips)          {k100_Detector->SetConstructZipBool(true);}
    else if(newValue == caseTower)    {k100_Detector->SetConstructTowerBool(true);}
    else if(newValue == caseVeto)     {k100_Detector->SetConstructVetoBool(true);}
    else if(newValue == caseShields)  {k100_Detector->SetConstructShieldsBool(true);}
    else if(newValue == caseFloor)  {k100_Detector->SetConstructFloorBool(true);}
    else if(newValue == caseWalls)  {k100_Detector->SetConstructWallsBool(true);}
    else if(newValue == caseCeiling)  {k100_Detector->SetConstructCeilingBool(true);}
    else if(newValue == caseWestReflector)  {k100_Detector->SetConstructWestReflectorBool(true);}
    else if(newValue == caseFrame)  {k100_Detector->SetConstructFrameBool(true);}
    else if(newValue == casePuBeSourceAndShield)  {k100_Detector->SetConstructPuBeSourceAndShieldBool(true);}
    else if(newValue == caseIceBox)   {k100_Detector->SetConstructIceBoxBool(true);}
    else if(newValue == caseThermalNeutronBucket)   {k100_Detector->SetConstructThermalNeutronBoxBool(true);}
    else if(newValue == caseGPSShielding)   {k100_Detector->SetConstructShieldTestEnvironmentBool(true);}
    else if(newValue == caseHPGeCoincidence)   {k100_Detector->SetConstructSimpleGammaCoinBool(true);}
    else if(newValue == caseNaIArray)   {k100_Detector->SetConstructNaIBool(true);} //spandey
    else if(newValue == caseBoronShield)   {k100_Detector->SetConstructBoronShieldBool(true);} //spandey
    else if(newValue == casePolyBox)   {k100_Detector->SetConstructPolyBox(true);} //spandey
  }

  if( command == DetectorDeActivateCmd ) { 
    if(newValue == caseZips)          {k100_Detector->SetConstructZipBool(false);}
    else if(newValue == caseTower)    {k100_Detector->SetConstructTowerBool(false);}
    else if(newValue == caseVeto)     {k100_Detector->SetConstructVetoBool(false);}
    else if(newValue == caseShields)  {k100_Detector->SetConstructShieldsBool(false);}
    else if(newValue == caseFloor)  {k100_Detector->SetConstructFloorBool(false);}
    else if(newValue == caseWalls)  {k100_Detector->SetConstructWallsBool(false);}
    else if(newValue == caseCeiling)  {k100_Detector->SetConstructCeilingBool(false);}
    else if(newValue == caseWestReflector)  {k100_Detector->SetConstructWestReflectorBool(false);}
    else if(newValue == caseFrame)  {k100_Detector->SetConstructFrameBool(false);}
    else if(newValue == casePuBeSourceAndShield)  {k100_Detector->SetConstructPuBeSourceAndShieldBool(false);}
    else if(newValue == caseIceBox)   {k100_Detector->SetConstructIceBoxBool(false);}
    else if(newValue == caseThermalNeutronBucket)   {k100_Detector->SetConstructThermalNeutronBoxBool(false);}
    else if(newValue == caseGPSShielding)   {k100_Detector->SetConstructShieldTestEnvironmentBool(false);}
    else if(newValue == caseHPGeCoincidence)   {k100_Detector->SetConstructSimpleGammaCoinBool(false);}
    else if(newValue == caseNaIArray)   {k100_Detector->SetConstructNaIBool(false);} //spandey
    else if(newValue == caseBoronShield)   {k100_Detector->SetConstructBoronShieldBool(false);} //spandey
    else if(newValue == casePolyBox)   {k100_Detector->SetConstructPolyBox(false);} //spandey

  }

  if( command == ZipConfigureCmd_Mat1 ) { 
    G4bool truth = ZipConfigureCmd_Mat1->GetNewBoolValue(newValue);
    k100_Detector->SetFirstDetGe(truth);
  }

  if( command == ShieldConfigureCmd_HPGeboron_shield ) { 
    G4bool truth = ShieldConfigureCmd_HPGeboron_shield->GetNewBoolValue(newValue);
    k100_Detector->SetConstructShields_HPGeboron_wshield(truth);
  }

  if( command == ShieldConfigureCmd_HPGeboron ) { 
    G4bool truth = ShieldConfigureCmd_HPGeboron->GetNewBoolValue(newValue);
    k100_Detector->SetConstructShields_HPGeboron(truth);
  }

  if( command == ShieldConfigureCmd_SouthNaI ) { 
    G4bool truth = ShieldConfigureCmd_SouthNaI->GetNewBoolValue(newValue);
    k100_Detector->SetConstructShields_addNaISouth(truth);
  }

  if( command == ShieldConfigureCmd_BasePoly ) { 
    G4bool truth = ShieldConfigureCmd_BasePoly->GetNewBoolValue(newValue);
    k100_Detector->SetConstructShields_addBasePoly(truth);
  }

  if( command == ShieldConfigureCmd_BaseLead ) { 
    G4bool truth = ShieldConfigureCmd_BaseLead->GetNewBoolValue(newValue);
    k100_Detector->SetConstructShields_addBaseLead(truth);
  }

  if( command == ShieldConfigureCmd_mod ) { 
    G4int Mod = ShieldConfigureCmd_mod->GetNewIntValue(newValue);
    k100_Detector->SetConstructShields_mod(Mod);
  }

  if( command == FridgeConfigureCmd_pure3HeBath ) { 
    G4bool truth = FridgeConfigureCmd_pure3HeBath->GetNewBoolValue(newValue);
    k100_Detector->SetConstructIceBox_pure3HeBath(truth);
  }

  if( command == PuBeConfigureCmd_westPolySensitive ) { 
    G4bool truth = PuBeConfigureCmd_westPolySensitive->GetNewBoolValue(newValue);
    k100_Detector->SetConstructPuBeSourceAndShield_westPolySensitivity(truth);
  }

  if( command == PuBeConfigureCmd_doPuBeGammas ) { 
    G4bool truth = PuBeConfigureCmd_doPuBeGammas->GetNewBoolValue(newValue);
    k100_Detector->SetConstructPuBeSourceAndShield_doPuBeGamma(truth);
  }

  if( command == PuBeConfigureCmd_Barrel ) { 
    G4bool truth = PuBeConfigureCmd_Barrel->GetNewBoolValue(newValue);
    k100_Detector->SetConstructPuBeSourceAndShield_addBarrel(truth);
  }

  if( command == PuBeConfigureCmd_NaIsensitive ) { 
    G4bool truth = PuBeConfigureCmd_NaIsensitive->GetNewBoolValue(newValue);
    k100_Detector->SetConstructPuBeSourceAndShield_setNaISensitive(truth);
  }

  if( command == PuBeConfigureCmd_R66 ) { 
    G4bool truth = PuBeConfigureCmd_R66->GetNewBoolValue(newValue);
    k100_Detector->SetConstructPuBeSourceAndShield_doR66(truth);
    if(truth){ //set all others to false
      k100_Detector->SetConstructPuBeSourceAndShield_doR62(false);
      k100_Detector->SetConstructPuBeSourceAndShield_doOrb(false);
    }
  }

  if( command == PuBeConfigureCmd_R62 ) { 
    G4bool truth = PuBeConfigureCmd_R62->GetNewBoolValue(newValue);
    k100_Detector->SetConstructPuBeSourceAndShield_doR62(truth);
    if(truth){ //set all others to false
      k100_Detector->SetConstructPuBeSourceAndShield_doR66(false);
      k100_Detector->SetConstructPuBeSourceAndShield_doOrb(false);
    }
  }

  if( command == PuBeConfigureCmd_Orb ) { 
    G4bool truth = PuBeConfigureCmd_Orb->GetNewBoolValue(newValue);
    k100_Detector->SetConstructPuBeSourceAndShield_doOrb(truth);
    if(truth){ //set all others to false
      k100_Detector->SetConstructPuBeSourceAndShield_doR66(false);
      k100_Detector->SetConstructPuBeSourceAndShield_doR62(false);
    }
  }

  if( command == PuBeConfigureCmd_mod ) { 
    G4int Mod = PuBeConfigureCmd_mod->GetNewIntValue(newValue);
    k100_Detector->SetConstructPuBeSourceAndShield_mod(Mod);
  }

  if( command == PuBeConfigureCmd_OrbPos ) { 
    G4ThreeVector pos = PuBeConfigureCmd_OrbPos->GetNew3VectorValue(newValue);
    k100_Detector->SetConstructPuBeSourceAndShield_OrbPos(pos);
  }

  if( command == PuBeConfigureCmd_OrbRad ) { 
    G4double rad = PuBeConfigureCmd_OrbRad->GetNewDoubleValue(newValue);
    k100_Detector->SetConstructPuBeSourceAndShield_OrbRad(rad);
  }

  if( command == GPSShieldPositionCmd ) { 
    G4ThreeVector pos = GPSShieldPositionCmd->GetNew3VectorValue(newValue);
    k100_Detector->SetConstructShieldTestEnvironmentPos(pos.x(),pos.y(),pos.z());
  }

  if( command == GPSShieldSizeCmd ) { 
    G4ThreeVector sz = GPSShieldSizeCmd->GetNew3VectorValue(newValue);
    k100_Detector->SetConstructShieldTestEnvironmentSize(sz.x(),sz.y(),sz.z());
  }

  if( command == GPSShieldMatCmd ) { 
    k100_Detector->SetConstructShieldTestEnvironmentMat(newValue);
  }

  if( command == DetSizeCmd ) { 
    G4ThreeVector sz = GPSShieldSizeCmd->GetNew3VectorValue(newValue);
    k100_Detector->DetSizeMod(sz.x(),sz.y());
  }

  if( command == GeGammaCoinPositionCmd ) { 
    G4ThreeVector pos = GeGammaCoinPositionCmd->GetNew3VectorValue(newValue);
    k100_Detector->SetConstructSimpleGammaCoinPos(pos.x(),pos.y(),pos.z());
  }

  if( command == GeGammaCoinSizeCmd ) { 
    G4ThreeVector sz = GeGammaCoinSizeCmd->GetNew3VectorValue(newValue);
    k100_Detector->SetConstructSimpleGammaCoinSize(sz.x(),sz.y());
  }

  if( command == SimpleGammaCoinMatCmd ) { 
    k100_Detector->SetConstructSimpleGammaCoinMat(newValue);
  }

  if( command == DrawSolidBox) { 
    if(newValue == caseZips)          {k100_Detector->SetDrawSolidZipBool(true);}
    else if(newValue == caseTower)    {k100_Detector->SetDrawSolidTowerBool(true);}
    else if(newValue == caseVeto)     {k100_Detector->SetDrawSolidVetoBool(true);}
    else if(newValue == caseShields)  {k100_Detector->SetDrawSolidShieldsBool(true);}
    else if(newValue == caseIceBox)   {k100_Detector->SetDrawSolidIceBoxBool(true);}
  }

  if( command == DrawFrameBox) { 
    if(newValue == caseZips)          {k100_Detector->SetDrawSolidZipBool(false);}
    else if(newValue == caseTower)    {k100_Detector->SetDrawSolidTowerBool(false);}
    else if(newValue == caseVeto)     {k100_Detector->SetDrawSolidVetoBool(false);}
    else if(newValue == caseShields)  {k100_Detector->SetDrawSolidShieldsBool(false);}
    else if(newValue == caseIceBox)   {k100_Detector->SetDrawSolidIceBoxBool(false);}
  }

  if( command == NbTowersCmd ) { 
    k100_Detector->SetNbOfTowers(NbTowersCmd->GetNewIntValue(newValue));
    k100_Detector->SetNbOfZips((NbTowersCmd->GetNewIntValue(newValue))*6);
  }

  if( command == SodiumBorateDensityFractionCmd) {
    G4double fraction = SodiumBorateDensityFractionCmd->GetNewDoubleValue(newValue);
    k100_Detector->SetSodiumBorateDensityFraction(fraction);
  }

  if( command == NbBoronVertCmd) {
    G4int NbBoronVert = NbBoronVertCmd->GetNewIntValue(newValue);
    //if(NbBoronVert > 2) G4cout<<" Can not place more than two vertical shields"<<G4endl;
    k100_Detector->SetNbBoronShieldVert(NbBoronVert);
  }

  if( command == NbBoronHoriCmd) {
    G4int NbBoronHori = NbBoronHoriCmd->GetNewIntValue(newValue);
    //if(NbBoronHori > 2) G4cout<<" Can not place more than two horizontal shields"<<G4endl;
    k100_Detector->SetNbBoronShieldHori(NbBoronHori);
  }

  if( command == BoronShieldThicknessCmd) {
    G4double thickness = BoronShieldThicknessCmd->GetNewDoubleValue(newValue);
    k100_Detector->SetBoronShieldThickness(thickness);
  }

  // Now call and update the detector
  //k100_Detector->UpdateGeometry();
}

// ------------------------------------------------

