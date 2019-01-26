/*==================k100_RunActionMessenger.cc====================== 
   
      PROGRAMMER:  Anthony Villano 10/18/12

      UPDATES:      
       

      PURPOSE: Code supporting the RunActionMessenger class including
               specific routines which run RunAction methods with various
	       input depending on the values passed via Geant4 macro.
              
======================================================================*/

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"

#include "G4RunManager.hh"

#include "k100_RunActionMessenger.hh"
#include "k100_RunAction.hh"


k100_RunActionMessenger::k100_RunActionMessenger(k100_RunAction* k100_Run):pRunAction(k100_Run)
{ 
  // Create the run directory
  k100_RunDir = new G4UIdirectory("/run/k100/");
  k100_RunDir->SetGuidance("k100 specific run controls.");

  //  run directory already exists

  // Set autoSeed command
  //setAutoSeedCmd = new G4UIcmdWithABool("/run/autoSeed",this);
  //setAutoSeedCmd->SetGuidance("Switch on/off time-based random seeds");
  //setAutoSeedCmd->SetGuidance(" true: run seeds determined by system time");
  //setAutoSeedCmd->SetGuidance("false: use command 'random/resetEngineFrom'");
  //setAutoSeedCmd->SetGuidance("Default = false");
  //setAutoSeedCmd->SetParameterName("autoSeed", false);
  //setAutoSeedCmd->AvailableForStates(G4State_Idle);

  // Set OutputDataToFile
  //setOutputDataToFile = new G4UIcmdWithABool("/run/writedata",this);
  //setOutputDataToFile->SetGuidance("Output results to file.");
  //setOutputDataToFile->SetParameterName("dataout",true,true); 
  // void SetParameterName(const char * theName,G4bool omittable,
  //                      G4bool currentAsDefault=false);
  //setOutputDataToFile->SetDefaultValue(true);
  //setOutputDataToFile->AvailableForStates(G4State_PreInit,G4State_Idle);

  // Set file name
  setRunFileName = new G4UIcmdWithAString("/run/k100/OFPrefix",this);
  setRunFileName->SetGuidance("Set the name of the output files.");
  setRunFileName->SetParameterName("fname",true);
  setRunFileName->SetDefaultValue("Sim");
  setRunFileName->AvailableForStates(G4State_PreInit,G4State_Idle);

  // Set drawEvent
  //  the  G4UIdirectory("/k100/"); directory is defined in DetectorMessenger
  setDrawEventCmd = new G4UIcmdWithABool("/run/k100/Draw",this);
  setDrawEventCmd->SetGuidance("Set drawFlag to Draw an event. (default = true)");
  setDrawEventCmd->SetParameterName("Draw", false);
  setDrawEventCmd->SetDefaultValue(true);

  setNCapOutputCmd = new G4UIcmdWithABool("/run/k100/OnlyNCapOut",this);
  setNCapOutputCmd->SetGuidance("Set true for only NCapture in Zip otuput. (default = false)");
  setNCapOutputCmd->SetParameterName("OnlyNCapOut", false);
  setNCapOutputCmd->SetDefaultValue(false);
}
k100_RunActionMessenger::~k100_RunActionMessenger()
{
   delete setDrawEventCmd;     
}
void k100_RunActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 


  if( command == setDrawEventCmd ) {
    if(setDrawEventCmd->GetNewBoolValue(newValue)) {
      G4cout << "Turning ON event drawing for all events." << G4endl; }
    else {
      G4cout << "Turning OFF event drawing for all events." << G4endl; }
    pRunAction->SetDrawEventCmd(setDrawEventCmd->GetNewBoolValue(newValue));
  }

  if( command == setNCapOutputCmd ) {
    G4bool truth = setNCapOutputCmd->GetNewBoolValue(newValue);
    pRunAction->SetSaveOnlyNCapture(truth);
  }

  if( command == setRunFileName ){
    pRunAction->SetDataFileNamePrefix(newValue);
  }

  
}
