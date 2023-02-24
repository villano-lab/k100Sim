/*==================k100_EventAction.cc============================ 
   
      PROGRAMMER:  Anthony Villano 10/18/12

      UPDATES:      
       

      PURPOSE:  Code support for the event actions which occur in
                the k100 simulation.
              
======================================================================*/

#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4DigiManager.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

#include "k100_EventAction.hh"
#include "k100_EventInfo.hh"
#include "k100_RunAction.hh"
#include "k100_DetectorConstruction.hh"
#include "k100_StdHit.hh"
#include "k100_DataStorage.hh"

k100_EventAction::k100_EventAction(G4bool nCapSaveOnly):drawEvent(true),saveEvent(true),saveOnlyNCapture(nCapSaveOnly)
{
	k100_Detector = 
    (k100_DetectorConstruction*)(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

}
k100_EventAction::k100_EventAction(k100_RunAction* RunAction,G4bool nCapSaveOnly):drawEvent(true),saveEvent(true),pRunAction(RunAction),saveOnlyNCapture(nCapSaveOnly)
{
	k100_Detector = 
    (k100_DetectorConstruction*)(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

	//make the run action information consistent
	pRunAction->SetSaveOnlyNCapture(saveOnlyNCapture);
}
k100_EventAction::~k100_EventAction()
{

}
void k100_EventAction::BeginOfEventAction(const G4Event* evt)
{
	//be sure the RunAction copy of saveOnlyNCapture is used
	saveOnlyNCapture = pRunAction->GetSaveOnlyNCapture();
	
	drawEvent = pRunAction->GetDrawEventCmd();
	saveEvent = true; // for now this is simple logic
	
        //if you have a sensitive detector, put the collection ID as k100CollID
	G4SDManager * SDman = G4SDManager::GetSDMpointer();

        if(evt->GetEventID()==1){
          std::map<G4String,G4int>::iterator sensitiveItr;
          std::map<G4String,G4int> *namemap = k100_Detector->GetSensitiveList();
          for(sensitiveItr=namemap->begin();sensitiveItr!=namemap->end();++sensitiveItr){
             G4int ID = SDman->GetCollectionID(sensitiveItr->first);
             G4cout << "#### EventAction::BeginOfEventAction : k100CollName for " 
               << sensitiveItr->first << G4endl;
             G4cout << "#### EventAction::BeginOfEventAction : k100CollID for " 
               << ID << G4endl;
          }
        }

        //evt->GetPrimaryVertex()->Print();

	// Periodic printing (should make this a verbosable item)
	if ((evt->GetEventID()%1000) == 1) {
		G4cout << "\n>>>>>> (BeginofEventAction) Event " << evt->GetEventID() << G4endl;
	}
	
}
void k100_EventAction::EndOfEventAction(const G4Event* evt)
{
        //dataOut is a k100_DataStorage object which is a member of
	//the k100_RunAction object.  It handles output data.
	dataOut = pRunAction->dataOut;
	//saveEvent=false;
	
	// Get number trajectories
	//
	G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
	G4int n_trajectories = 0;
	if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

	// Periodic printing (should make this a verbosable item)
	//if (((evt->GetEventID()%1000)==1)) {
	if (trajectoryContainer&&((evt->GetEventID()%1000)==1)) {
		G4cout << "\n>>>>>> (EndOfEventAction) Got trajectory container with " << n_trajectories << " hits." << G4endl;
	}

        //if you have a sensitive detector, put the collection ID as k100CollID
	G4SDManager * SDman = G4SDManager::GetSDMpointer();
        std::map<G4String,G4int>::iterator sensitiveItr;
        std::map<G4String,G4int> *namemap = k100_Detector->GetSensitiveList();
        for(sensitiveItr=namemap->begin();sensitiveItr!=namemap->end();++sensitiveItr){
          G4int ID = SDman->GetCollectionID(sensitiveItr->first);
          G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
          k100_StdHitsCollection* GHC = 0;
          if (HCE) {
             //G4cout << "\n>>> (EndofEventAction) HCE NumColls = " << HCE->GetNumberOfCollections() << G4endl;
             if(k100_Detector->GetConstructGenericSensitiveInt()>0) {
	       GHC = (k100_StdHitsCollection*)(HCE->GetHC(ID));
             } 


           }
           else {
             	std::cout<<"No hit collection for event = "<<evt->GetEventID()<<std::endl;
             } 
 
           // All of the following happens only if there is a hit in the Standard Sensitive-Dets.
           G4int n_hit = -1; // Initialize it with  value of -1

           if(GHC) {
	     //G4cout << "\n>>> (EndofEventAction) GHC entries =  " << GHC->entries() << G4endl;
             if(GHC->entries()) {
	       n_hit = GHC->entries();

	       drawEvent = true;
	       saveEvent = true;
             }
	     else
	       n_hit=0;
           }

           //check the event to be sure at least one hit has an ncapture (parameter 22 of data array unity)
	   G4bool gotNCap = false;
           for (G4int ii=0; ii<n_hit; ii++)
	   {
	     if(!saveOnlyNCapture) break; //skip this if saving all events
	     G4double *quickdata = (*GHC)[ii]->GetData();
	     if(quickdata[22-1]==1){
	       gotNCap = true;
	       break;
	     }
	   }
	   if(!gotNCap && saveOnlyNCapture) saveEvent=false;
    
           if(saveEvent) {


             //check to see if storage of this detector event will overflow current array
             if (dataOut->overflowArray(n_hit))
	     {
	       //if so, then write the array to a file
	       dataOut->writeArray();
	     }
    
             //actual storage loop
             for (G4int ii=0; ii<n_hit; ii++)
	     {
	       ////G4cout << "\n>>> (EndofEventAction) I am HERE 3 " << evt->GetEventID() << G4endl;    

	       //set the current-data structure to contain this hit information
	       dataOut->setData((*GHC)[ii]->GetData());
	  
	       // add aditional data that we can only get here
	       // getData() returns a pointer to a data array so that we can
	       // modify the memory
	       dataOut->getData()[1-1] = evt->GetEventID()+1;
               k100_EventInfo *anEventInfo = (k100_EventInfo*) evt->GetUserInformation();
	       //FIXME this seems like a hack for gamma-Be running!  remove it
	       //dataOut->getData()[18-1] = anEventInfo->GetBeLength();

	       // add this data to the matrix array
	       dataOut->addData();
	     } 
          }
        }
 	
	// Extract the trajectories and draw them
	//
	if(drawEvent){
		if (G4VVisManager::GetConcreteInstance()) {
			for (G4int ii=0; ii<n_trajectories; ii++) {
				G4Trajectory* trj = (G4Trajectory*)
				((*(evt->GetTrajectoryContainer()))[ii]);
				trj->DrawTrajectory();
			}
		}
		
	}    
	
}
