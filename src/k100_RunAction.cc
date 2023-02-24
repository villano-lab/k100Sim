/*==================k100_RunAction.cc============================== 
   
      PROGRAMMER:  Anthony Villano 10/18/12

      UPDATES:      
       

      PURPOSE: Code supporting the RunAction class including specific
               routines for setting file name prefix and beinginning/end
	       of run actions.
              
======================================================================*/

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"

#include "Randomize.hh"

#include "k100_RunAction.hh"
#include "k100_RunActionMessenger.hh"
#include "k100_DetectorConstruction.hh"
#include "k100_DataStorage.hh"
#include <sys/time.h>
#include <time.h>
char filename[200];


k100_RunAction::k100_RunAction(G4bool rootOutput)
{
 
  //set a default for saveOnlyNCapture
  saveOnlyNCapture = false;
  
  //save output in textFile?
  OutputRootFlag = rootOutput;
  // automatic (time-based) random seeds and filenames for each run
  struct timeval mytime;
  gettimeofday(&mytime, NULL);
  randSeed=mytime.tv_sec-mytime.tv_usec;
  CLHEP::HepRandom::setTheSeed(randSeed);
  CLHEP::HepRandom::showEngineStatus();

  // Set a default file prefix
  sprintf(filename,"k100_%ld",randSeed);
  DataFileNamePrefix = G4String(filename);
  OutputDataToFile = true;
  DrawEventCmd = true;
  SaltPillOutCmd = false;
  ResetRun=false;

  runMessenger = new k100_RunActionMessenger(this);  
}
k100_RunAction::~k100_RunAction()
{

}
void k100_RunAction::BeginOfRunAction(const G4Run* aRun)
{
  


  runN = aRun->GetRunID();
  if ( (runN % 1000 == 0) || (runN<100) ) 
    G4cout << "### Run : " << runN << G4endl;

  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer();
      UI->ApplyCommand("/vis/scene/notifyHandlers");
    }

  // ------ Determine data output file name -----

  // There MUST be a better way of doing this !!!
  char fnum[5];
  sprintf(fnum,"%03d",runN);

  if(OutputDataToFile) {

    // Create the DataOut instance -- with the current data file names
    char file[200];
    sprintf(file,"%ld",randSeed);
    if(runN == 0 || ResetRun){
       DataFileNamePrefix=DataFileNamePrefix+G4String(file);
       ResetRun=false;
     }
    //DataFileNamePrefix=DataFileNamePrefix+G4String(file);
    dataOut = new k100_DataStorage(DataFileNamePrefix,runN,1,OutputRootFlag); 

    #ifdef NON_SD_INFO
    // For Non SD nCap info
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    man->OpenFile("Nsc_"+DataFileNamePrefix);
    man->CreateNtuple("nCapInfo","nCapInfo");
    man->CreateNtupleIColumn("fEvent");
    man->CreateNtupleDColumn("fNeutron_energy");
    man->CreateNtupleIColumn("fPDGID");
    man->CreateNtupleDColumn("fsec_KE");

    man->FinishNtuple(0);

    std::cout<<"################## DataFileNamePrefix = "<<DataFileNamePrefix<<std::endl;
    #endif
  }




}
void k100_RunAction::EndOfRunAction(const G4Run*)
{

  if (G4VVisManager::GetConcreteInstance())
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");

  //Kill the data out class instance
  if(OutputDataToFile) {
    delete dataOut;
  }

  #ifdef NON_SD_INFO
  // For Non SD nCap info
  G4AnalysisManager *man = G4AnalysisManager::Instance(); 
  man->Write();
  man->CloseFile();
  #endif

}
