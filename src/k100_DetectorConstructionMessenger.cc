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
#include "G4UIcmdWithoutParameter.hh"

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
  DetectorActivateCmd->SetGuidance("Choices are : Zips/Towers/Veto/Shields/IceBox .");
  DetectorActivateCmd->SetParameterName("choice",false);
  DetectorActivateCmd->AvailableForStates(G4State_Idle);

  DetectorDeActivateCmd = new G4UIcmdWithAString("/CDMS/detector/deactivate",this);
  DetectorDeActivateCmd->SetGuidance("Dectivate CDMS Detector Element.");
  DetectorDeActivateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  DetectorDeActivateCmd->SetGuidance("in order for change to take effect.");
  DetectorDeActivateCmd->SetGuidance("Choices are : Zip/Tower/Veto/Shields/IceBox .");
  DetectorDeActivateCmd->SetParameterName("choice",false);
  DetectorDeActivateCmd->AvailableForStates(G4State_Idle);

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

  // Set parameters for gamma-coincidence detector
  GeGammaCoinPositionCmd = new G4UIcmdWith3VectorAndUnit("/CDMS/gammaCoin/setPosition",this);
  GeGammaCoinPositionCmd->SetGuidance("Set position for HPGe coincidence detector.");
  GeGammaCoinPositionCmd->SetParameterName("Dx","Dy","Dz",true,true);
  GeGammaCoinPositionCmd->SetRange("Dx != 0 || Dy != 0 || Dz != 0"); //FIXME actually might want to put more constraints? no overlap?

  GeGammaCoinSizeCmd = new G4UIcmdWith3VectorAndUnit("/CDMS/gammaCoin/setSize",this);
  GeGammaCoinSizeCmd->SetGuidance("Set size for HPGe coincidence detector.");
  GeGammaCoinSizeCmd->SetParameterName("Radius","Thickness","Nothing",true,true);
  GeGammaCoinSizeCmd->SetRange("Radius != 0 && Thickness != 0"); //none of these can be zero

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
    else if(newValue == caseIceBox)   {k100_Detector->SetConstructIceBoxBool(true);}
    else if(newValue == caseThermalNeutronBucket)   {k100_Detector->SetConstructThermalNeutronBoxBool(true);}
    else if(newValue == caseGPSShielding)   {k100_Detector->SetConstructShieldTestEnvironmentBool(true);}
    else if(newValue == caseHPGeCoincidence)   {k100_Detector->SetConstructSimpleGammaCoinBool(true);}
  }

  if( command == DetectorDeActivateCmd ) { 
    if(newValue == caseZips)          {k100_Detector->SetConstructZipBool(false);}
    else if(newValue == caseTower)    {k100_Detector->SetConstructTowerBool(false);}
    else if(newValue == caseVeto)     {k100_Detector->SetConstructVetoBool(false);}
    else if(newValue == caseShields)  {k100_Detector->SetConstructShieldsBool(false);}
    else if(newValue == caseIceBox)   {k100_Detector->SetConstructIceBoxBool(false);}
    else if(newValue == caseThermalNeutronBucket)   {k100_Detector->SetConstructThermalNeutronBoxBool(false);}
    else if(newValue == caseGPSShielding)   {k100_Detector->SetConstructShieldTestEnvironmentBool(false);}
    else if(newValue == caseHPGeCoincidence)   {k100_Detector->SetConstructSimpleGammaCoinBool(false);}

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

  if( command == GeGammaCoinPositionCmd ) { 
    G4ThreeVector pos = GeGammaCoinPositionCmd->GetNew3VectorValue(newValue);
    k100_Detector->SetConstructSimpleGammaCoinPos(pos.x(),pos.y(),pos.z());
  }

  if( command == GeGammaCoinSizeCmd ) { 
    G4ThreeVector sz = GeGammaCoinSizeCmd->GetNew3VectorValue(newValue);
    k100_Detector->SetConstructSimpleGammaCoinSize(sz.x(),sz.y());
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

  // Now call and update the detector
  k100_Detector->UpdateGeometry();
}

// ------------------------------------------------

