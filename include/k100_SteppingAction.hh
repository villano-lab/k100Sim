#ifndef k100_SteppingAction_h
#define k100_SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "k100_EventAction.hh"

class B1EventAction;

class G4LogicalVolume;

/// Stepping action class                                                                                                                                                                                       
///                                                                                                                                                                                                             

class k100_SteppingAction : public G4UserSteppingAction
{
  public:
    k100_SteppingAction(k100_EventAction* eventAction);
    virtual ~k100_SteppingAction();

    // method from the base class                                                                                                                                                                               
    virtual void UserSteppingAction(const G4Step*);

  private:
    // B1EventAction*  fEventAction;
    // G4LogicalVolume* fScoringVolume;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                                                                                                                

#endif