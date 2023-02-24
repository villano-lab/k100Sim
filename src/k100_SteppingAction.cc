#include "k100_SteppingAction.hh"
#include "k100_EventAction.hh"
#include "k100_DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

//#include "G4TrackVector.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                                                                                                                

k100_SteppingAction::k100_SteppingAction(k100_EventAction* eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                                                                                                                

k100_SteppingAction::~k100_SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                                                                                                                

void k100_SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  
  if(aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "nCapture") {
    const G4TrackVector* tkV  = aStep->GetSecondary();
    int nsc = (*tkV).size();
    //std::cout<<"========> nCapture : nsc = "<<nsc<<std::endl;
    // kinetic energy : aStep->GetPreStepPoint()->GetKineticEnergy()
    G4double neutron_energy = aStep->GetPreStepPoint()->GetKineticEnergy();
    for(int i = 0; i < nsc; i++) {
      G4Track *tk = (*tkV)[i];
      if(tk->GetCreatorProcess()->GetProcessName() != "nCapture") continue;
      //std::cout<<"  ========> Daughters: "<< i <<" , PDGID: "<<tk->GetDefinition()->GetPDGEncoding()<<std::endl;
      G4int evnt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

      #ifdef NON_SD_INFO
      // For Non SD nCap info
      G4AnalysisManager *man = G4AnalysisManager::Instance();
      man->FillNtupleIColumn(0,evnt);
      man->FillNtupleDColumn(1,neutron_energy);
      man->FillNtupleIColumn(2,tk->GetDefinition()->GetPDGEncoding()%1000000);
      man->FillNtupleDColumn(3,tk->GetKineticEnergy());
      man->AddNtupleRow(0);
      #endif
    }
    
  } 
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                                                                                                                


