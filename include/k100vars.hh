// ------------------------------------------------
//
// k100vars.hh : July 2010
//
// ------------------------------------------------

#ifndef k100vars_H
#define k100vars_H 1


#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4Region.hh"


extern G4double frame_x=-10.*2.54*cm, frame_y=-12.0*2.54*cm, frame_z=0.*cm;
extern G4double tower_x=0., tower_y=0., tower_z=-6.5*cm;
//G4double tower_x=0., tower_y=0., tower_z=1.6643*cm; //updated 12/28/17 to get distance from floor right
extern G4double fridge_x=0., fridge_y=0., fridge_z=0.;
//G4double fridge_x=0., fridge_y=0., fridge_z=8.1643*cm; //updated 12/28/17 to get distance from floor right

#endif