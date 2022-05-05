// ------------------------------------------------
//
// k100_ZipSD.cc : Dec 2005
//
// ------------------------------------------------

#include "G4UnitsTable.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "G4Neutron.hh"

#include "k100_ZipSD.hh"
#include <string.h>
// ------------------------------------------------

k100_ZipSD::k100_ZipSD(G4String name, G4int towerNb)
 :G4VSensitiveDetector(name),towNb(towerNb)
{
  G4String HCname;// = "ZipCollection";
  HCname = name;
  collectionName.insert(HCname);
}

// ------------------------------------------

k100_ZipSD::~k100_ZipSD() {}

// ------------------------------------------

void k100_ZipSD::Initialize(G4HCofThisEvent* HCE)
{
  zipCollection = new k100_ZipHitsCollection(SensitiveDetectorName,collectionName[0]);

  G4int HCID = -1;
  if(HCID < 0) {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    //    G4cout << ">>>>>> k100_ZipSD::Initialize HCID = " << HCID << G4endl;
  }
  HCE->AddHitsCollection(HCID, zipCollection);
}

// ------------------------------------------

G4bool k100_ZipSD::ProcessHits(G4Step* aStep, G4TouchableHistory* /*ROhist*/)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double pid;

 // if(edep==0.) return false;

  k100_ZipHit* newHit = new k100_ZipHit();
  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  newHit->SetZipNb(aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber());

  //  G4cout << ">>>>>> ***** (k100_ZipSD::ProcessHits) ZipNb = " 
  //         << aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber() << G4endl;
  newHit->SetEdep(edep);
  newHit->SetGlobalTime(aStep->GetPostStepPoint()->GetGlobalTime());
  newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
  newHit->SetPName(aStep->GetTrack()->GetDefinition()->GetParticleName());

  G4double* dataVector = newHit->GetData();

  if(aStep->GetTrack()->GetDefinition()->GetParticleType() == "nucleus") {
    pid = (aStep->GetTrack()->GetDefinition()->GetPDGCharge()) +
          (1000 * aStep->GetTrack()->GetDefinition()->GetBaryonNumber());
    //  G4cout << "\n>>>> k100_ZipSD::ProcessHits I am a Nucleus! " << G4endl;
    //  G4cout << ">>>> k100_ZipSD::ProcessHits GetPDGEncoding = " << pid << G4endl;
  }
  else {
    pid = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
    //  G4cout << "\n>>>> k100_ZipSD::ProcessHits I am a Particle! " << G4endl;
    //  G4cout << ">>>> k100_ZipSD::ProcessHits GetPDGEncoding = " << pid << G4endl;
  }

  dataVector[1-1]  = 0; // Is not set in this function
  if(!strcmp(aStep->GetPreStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName(),"NaI_tile_LV")) {
    if(abs(aStep->GetPreStepPoint()->GetPosition().x()) < 0.5*(324.0 + 20.0) ) { // For bottom NaI array
      dataVector[2-1] = 1000 + aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber() + 1;
    }
    else { // For vertical NaI array
      dataVector[2-1] = 2000 + aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber() + 1;
    }
  }
  else { // For Si/Ge crystal
    dataVector[2-1]  = (6*(towNb-1)) + 
                     aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber() + 1;
  }
  
  dataVector[3-1]  = 1e5 * aStep->GetTrack()->GetTrackID() +
                           aStep->GetTrack()->GetCurrentStepNumber();
  dataVector[4-1]  = 1e5 * aStep->GetTrack()->GetParentID();
  dataVector[5-1]  = pid; // Particle ID
  dataVector[6-1]  = aStep->GetPreStepPoint()->GetKineticEnergy();
  dataVector[7-1]  = aStep->GetTotalEnergyDeposit();

  dataVector[8-1]  = aStep->GetPostStepPoint()->GetMomentum().x();
  dataVector[9-1]  = aStep->GetPostStepPoint()->GetMomentum().y();
  dataVector[10-1] = aStep->GetPostStepPoint()->GetMomentum().z();
  dataVector[11-1] = aStep->GetPostStepPoint()->GetPosition().x();
  dataVector[12-1] = aStep->GetPostStepPoint()->GetPosition().y();
  dataVector[13-1] = aStep->GetPostStepPoint()->GetPosition().z();
  dataVector[14-1] = aStep->GetPostStepPoint()->GetGlobalTime();

  dataVector[15-1] = aStep->GetPreStepPoint()->GetMomentum().x();
  dataVector[16-1] = aStep->GetPreStepPoint()->GetMomentum().y();
  dataVector[17-1] = aStep->GetPreStepPoint()->GetMomentum().z();
  dataVector[18-1] = aStep->GetPreStepPoint()->GetPosition().x();
  dataVector[19-1] = aStep->GetPreStepPoint()->GetPosition().y();
  dataVector[20-1] = aStep->GetPreStepPoint()->GetPosition().z();
  dataVector[21-1] = aStep->GetPreStepPoint()->GetGlobalTime();

  //get some process information
  G4Track* track = aStep->GetTrack();

  dataVector[22-1] = isNCap(track,aStep);


  // std::cout<<"XXXXXXX"<<std::endl;
  // std::cout<<"X1 : Y1 : Z1 : VolumeName : replicaNumber : DetID :: "<<dataVector[18-1]<<" : "<<dataVector[19-1]<<" : "<<dataVector[20-1]<<" : "<< aStep->GetPreStepPoint()->GetTouchable()->GetVolume()->GetLogicalVolume()->GetName()<<" : "<<aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber()<<" : "<<dataVector[2-1]<<std::endl;

  //for(int i=0;i<21;i++){
  //  G4cout << "vector[" << i << "]: " << dataVector[i] << G4endl;
  //}
  /*
     1. EV (starts with 1)
     2. DT (detector #)
     3. TS (track and step: ttttsssss)
     4. PA (parent track  : tttt00000)
     5. TY (type. -ve =>  : -zzzaaa)
     6. E1 (KE of this track@1)
     7. D3 (energy deposition@3)
     8. PX3 (momentum@3)
     9. PY3
    10. PZ3
    11. X3 (position@3)
    12. Y3
    13. Z3
    14. GT (global time)
  */

  zipCollection->insert(newHit);

  //  newHit->Print();
  newHit->Draw();

  return true;
}

// ------------------------------------------

void k100_ZipSD::EndOfEvent(G4HCofThisEvent* /*HCE*/)
{
  if (verboseLevel>0) {
    G4int NbHits = zipCollection->entries();     
    G4cout << "\n--------> Hits Collection: in this event there are " << NbHits
	   << " hits in the Zip: " << G4endl;

    G4double TotEnergyDep = 0;

    for (G4int ii=0; ii<NbHits; ii++) {
      //(*zipCollection)[ii]->Print();
      TotEnergyDep += (*zipCollection)[ii]->GetEdep();
    }

    G4cout << "--------> Total deposited energy in " 
	   << SensitiveDetectorName << " : " 
	   << std::setw(5) << G4BestUnit(TotEnergyDep,"Energy")
	   << "\n" << G4endl;

  }
}

// ------------------------------------------

G4bool k100_ZipSD::isNCap(G4Track *track, G4Step *aStep){

  static const G4ParticleDefinition* theNeutron = G4Neutron::Definition();

  // Must be a neutron which has been stopped by a nucleus
  if (track->GetDefinition() != theNeutron) return false;
  if (track->GetTrackStatus() != fStopAndKill) return false;
  
  // Only tag if there were no neutron secondaries (identifies a capture)
  const G4TrackVector* secondaries = aStep->GetSecondary();
  size_t nsecs = secondaries->size();
  for(size_t i=0; i < nsecs; ++i) {    //quicker to search in reverse order
    if ((*secondaries)[nsecs-1-i]->GetDefinition() == theNeutron) return false;
  }

  return true;
}
