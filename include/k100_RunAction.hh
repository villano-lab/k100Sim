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
//#define NON_SD_INFO  //uncomment this line if you want to save all nCap info from the lab
#ifdef NON_SD_INFO
#include "g4root.hh"
#endif

class G4Run;

class k100_DataStorage;
class k100_RunActionMessenger;

class k100_RunAction : public G4UserRunAction
{
public :

  k100_RunAction(G4bool rootOutput);
  ~k100_RunAction();

public :

  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);

  void SetDataFileNamePrefix(G4String nPrefix) {DataFileNamePrefix = nPrefix;}
  void SetOutputDataToFile(G4bool dataout)     {OutputDataToFile = dataout;}
  void SetSaveOnlyNCapture(G4bool value)     {saveOnlyNCapture = value;}
  void SetDrawEventCmd(G4bool drawBool)        {DrawEventCmd = drawBool;}
  void SetSaltPillOutCmd(G4bool value)         {SaltPillOutCmd = value;}
  void SetAutoSeed (const G4bool value)        {autoSeed    = value;}

  G4String GetDataFileNamePrefix() {return DataFileNamePrefix;}
  G4bool GetSaveOnlyNCapture() {return saveOnlyNCapture;}
  G4bool GetOutputDataToFile()   {return OutputDataToFile;}
  G4bool GetDrawEventCmd()       {return DrawEventCmd;}
  G4bool GetSaltPillOutCmd()     {return SaltPillOutCmd;}
  G4int  GetRunNumber()          {return runN;}

  k100_DataStorage* dataOut;

private :

  k100_RunActionMessenger* runMessenger; 

  G4bool OutputRootFlag;
  G4bool   autoSeed;
  G4bool   saveOnlyNCapture; //variable to cull the output list
  G4String DataFileNamePrefix;
  G4bool   DrawEventCmd;
  G4bool   SaltPillOutCmd;
  G4bool   OutputDataToFile;
  G4bool   ResetRun;
  G4int    runN;
  long     randSeed;

};

// ------------------------------------------


#endif
