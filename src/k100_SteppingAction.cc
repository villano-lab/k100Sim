#include "k100_SteppingAction.hh"
#include "k100_EventAction.hh"
#include "k100_DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                                                                                                                

k100_SteppingAction::k100_SteppingAction(k100_EventAction* eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                                                                                                                

k100_SteppingAction::~k100_SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                                                                                                                

void k100_SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  
  

  // // if(aStep->GetTrack()->GetDefinition()->GetPDGEncoding() == 22 && aStep->GetTrack()->GetParentID() == 0) {
  // //   if(aStep->GetTrack()->GetNextVolume() == nullptr) {
  // //     std::cout<<"\t \t trackID : KE : x : y : z : currentVol : NextStepVol : Edep :: "<<aStep->GetTrack()->GetTrackID()<<" : "<<aStep->GetPreStepPoint()->GetKineticEnergy()<<" : " << aStep->GetPreStepPoint()->GetPosition().x()<<" : "<< aStep->GetPreStepPoint()->GetPosition().y() << " : " << aStep->GetPreStepPoint()->GetPosition().z()  << " : "<< aStep->GetTrack()->GetVolume()->GetName() << " : NULL : " << aStep->GetTotalEnergyDeposit() <<std::endl;
  // //   } else {
  // //     std::cout<<"\t \t trackID : KE : x : y : z : currentVol : NextStepVol : Edep :: "<<aStep->GetTrack()->GetTrackID()<<" : "<<aStep->GetPreStepPoint()->GetKineticEnergy()<<" : " << aStep->GetPreStepPoint()->GetPosition().x()<<" : "<< aStep->GetPreStepPoint()->GetPosition().y() << " : " << aStep->GetPreStepPoint()->GetPosition().z()  << " : "<< aStep->GetTrack()->GetVolume()->GetName() << " : " <<aStep->GetTrack()->GetNextVolume()->GetName() << " : " << aStep->GetTotalEnergyDeposit() <<std::endl;
  // //   }
    
  // // }

  // // if(((aStep->GetPostStepPoint() == nullptr) || (aStep->GetTrack()->GetNextVolume() == nullptr)) &&
  // //         (aStep->IsLastStepInVolume())) {
  // //   std::cout<<"*************** Can you hear me??? ***************"<<std::endl;
  // //   std::cout<<"pid : KE :: "<<aStep->GetTrack()->GetDefinition()->GetPDGEncoding()
  // //             <<" : " <<aStep->GetPreStepPoint()->GetKineticEnergy()<<std::endl;
  // }

  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                                                                                                                


