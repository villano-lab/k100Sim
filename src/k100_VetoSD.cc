// ------------------------------------------------
//
// k100_VetoSD.cc : Dec 2005
//
// ------------------------------------------------

#include "G4UnitsTable.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "k100_VetoSD.hh"

// ------------------------------------------------

k100_VetoSD::k100_VetoSD(G4String name, G4int panelNb):G4VSensitiveDetector(name), panelNb(panelNb)
{
  G4String HCname;// = "VetoCollection";
  HCname = name;
  collectionName.insert(HCname);
}

// ------------------------------------------

k100_VetoSD::~k100_VetoSD() {}

// ------------------------------------------

void k100_VetoSD::Initialize(G4HCofThisEvent* HCE)
{
  //  G4cout << ">>>>>> k100_VetoSD::Initialize SensitiveDetectorName =  " 
  //         << SensitiveDetectorName << G4endl;
  //  G4cout << ">>>>>> k100_VetoSD::Initialize  collectionName = " 
  // << collectionName[0] << G4endl;
  vetoCollection = new k100_VetoHitsCollection(SensitiveDetectorName,collectionName[0]);

  G4int HCID = -1;
  if(HCID < 0) {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    //    G4cout << ">>>>>> k100_VetoSD::Initialize HCID = " << HCID << G4endl;
  }
  HCE->AddHitsCollection(HCID, vetoCollection);
}

// ------------------------------------------

G4bool k100_VetoSD::ProcessHits(G4Step* aStep, G4TouchableHistory* /*ROhist*/)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double pid;

  if(edep==0.) return false;

  k100_VetoHit* newHit = new k100_VetoHit();
  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  newHit->SetVetoNb(aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber());

  //  G4cout << ">>>>>> ***** (k100_VetoSD::ProcessHits) VetoNb = " 
  // << aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber() << G4endl;
  newHit->SetEdep   (edep);
  newHit->SetGlobalTime(aStep->GetPostStepPoint()->GetGlobalTime());
  newHit->SetPos    (aStep->GetPostStepPoint()->GetPosition());
  newHit->SetPName  (aStep->GetTrack()->GetDefinition()->GetParticleName());

  G4double* dataVector = newHit->GetData();

  if(aStep->GetTrack()->GetDefinition()->GetParticleType() == "nucleus") {
    pid = (aStep->GetTrack()->GetDefinition()->GetPDGCharge()) + 
          (1000 * aStep->GetTrack()->GetDefinition()->GetBaryonNumber());
    //  G4cout << "\n>>>> k100_VetoSD::ProcessHits I am a Nucleus! " << G4endl;
    //  G4cout << ">>>> k100_VetoSD::ProcessHits GetPDGEncoding = " << pid << G4endl;
  }
  else {
    pid = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
    //  G4cout << "\n>>>> k100_VetoSD::ProcessHits I am a Particle! " << G4endl;
    //  G4cout << ">>>> k100_VetoSD::ProcessHits GetPDGEncoding = " << pid << G4endl;
  }

  dataVector[1-1]  = 0; // Is not set in this function
  dataVector[2-1]  = 1000 + panelNb; // Veto Panel Number
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

  /*
     1. EV (starts with 1)
     2. DT (detector #)
     3. TS (track and step: ttttsssss)
     4. PA (parent track  : tttt00000)
     5. TY (type. -ve =>  : -zzzaaa)
     6. E1 (KE of this track@1)
     7. D3 (energy deposition@3)
     8. PX3 (position@1)
     9. PY3
     10. PZ3
     11. X3 (position@3)
     12. Y3
     13. Z3
     14  GT (global time)
  */

  vetoCollection->insert(newHit);

  //  newHit->Print();
  newHit->Draw();

  return true;
}

// ------------------------------------------

void k100_VetoSD::EndOfEvent(G4HCofThisEvent* /*HCE*/)
{
  if (verboseLevel>0) {
    G4int NbHits = vetoCollection->entries();     
    G4cout << "\n--------> Hits Collection: in this event there are " << NbHits
	   << " hits in the Veto: " << G4endl;

    G4double TotEnergyDep = 0;

    for (G4int ii=0; ii<NbHits; ii++) {
      //(*vetoCollection)[ii]->Print();
      TotEnergyDep += (*vetoCollection)[ii]->GetEdep();
    }

    G4cout << "--------> Total deposited energy in " 
	   << SensitiveDetectorName << " : " 
	   << std::setw(5) << G4BestUnit(TotEnergyDep,"Energy")
	   << "\n" << G4endl;

  }
}

// ------------------------------------------
