/*==================k100_RunAction.hh================================= 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE: Class to keep track of runs as instantiated from a Geant4
               macro.  Every /run/beamOn n command begins a new run of
	       n thrown events.
              
======================================================================*/

#ifndef k100_RunAction_h
#define k100_RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class k100_DataStorage;
class k100_RunActionMessenger;
class k100_PrimaryGeneratorAction;

class k100_RunAction : public G4UserRunAction
{
public :

  k100_RunAction(k100_PrimaryGeneratorAction*);
  ~k100_RunAction();

public :

  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);

  void SetDataFileNamePrefix(G4String nPrefix) {DataFileNamePrefix = nPrefix;}
  void SetOutputDataToFile(G4bool dataout)     {OutputDataToFile = dataout;}
  void SetDrawEventCmd(G4bool drawBool)        {DrawEventCmd = drawBool;}
  void SetSaltPillOutCmd(G4bool value)         {SaltPillOutCmd = value;}
  void SetAutoSeed (const G4bool value)        {autoSeed    = value;}

  G4bool GetDataFileNamePrefix() {return DataFileNamePrefix;}
  G4bool GetOutputDataToFile()   {return OutputDataToFile;}
  G4bool GetDrawEventCmd()       {return DrawEventCmd;}
  G4bool GetSaltPillOutCmd()     {return SaltPillOutCmd;}
  G4int  GetRunNumber()          {return runN;}

  k100_DataStorage* dataOut;

private :

  k100_RunActionMessenger* runMessenger; 

  G4bool   autoSeed;
  G4String DataFileNamePrefix;
  G4bool   DrawEventCmd;
  G4bool   SaltPillOutCmd;
  G4bool   OutputDataToFile;
  G4int    runN;
  long     randSeed;
  k100_PrimaryGeneratorAction* thisgenerator;

};

// ------------------------------------------


#endif
