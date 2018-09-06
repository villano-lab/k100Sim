// ------------------------------------------------
//
// k100_ZipParametrization.hh
//
// ------------------------------------------------

#ifndef k100_ZipParameterisation_H
#define k100_ZipParameterisation_H 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalConstants.hh"
#include "G4VPVParameterisation.hh"

class G4VPhysicalVolume;
class G4Tubs;
class G4Material;
class G4VisAttributes;

// Dummy declarations to get rid of warnings ... (based on the 4.6.2 N02 example)
class G4Box;
class G4Trd;
class G4Trap;
class G4Cons;
class G4Orb;
class G4Sphere;
class G4Torus;
class G4Para;
class G4Hype;
class G4Polycone;
class G4Polyhedra;

// ------------------------------------------------
// ------------------------------------------------

class k100_ZipParameterisation : public G4VPVParameterisation
{ 

public:

  k100_ZipParameterisation(G4double TotalThk,
		           G4int NbZipsPerTower, 
			   G4double ZSpacing,  
			   G4double ZipRadius,   
			   G4double ZipDepth,
 			   G4int* DetMaterials,
 			   G4Material* zipGeMat, 
 			   G4Material* zipSiMat,
			   G4bool DrawSolidZipBool,
 			   G4int towerNb);

  virtual ~k100_ZipParameterisation();

  void ComputeTransformation (const G4int copyNo, G4VPhysicalVolume* physVol) const;

  void ComputeDimensions (G4Tubs&,const G4int,const G4VPhysicalVolume*) const ;

  //get the coordinates for each copy
  G4ThreeVector GetCoordinates(G4int copyNo);

private:  // Dummy declarations to get rid of warnings ...

  void ComputeDimensions (G4Box & detectorPixel, const G4int copyNo, const G4VPhysicalVolume* physVol) const {};
  void ComputeDimensions (G4Trd&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Trap&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Cons&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Sphere&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Orb&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Torus&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Para&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Hype&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Polycone&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Polyhedra&,const G4int,const G4VPhysicalVolume*) const {}


  virtual G4Material* ComputeMaterial(const G4int copyNo, G4VPhysicalVolume* physVol, const G4VTouchable *parentTouch=0);

private:

  G4int    fNbZipsPerTower;   
  G4double fStartX,fStartY,fStartXY;
  G4double fRadius;    //  The half-width of each pixel
  G4double fHalfDepth;    //  The half-depth of each pixel
  G4double fZSpacing;      //  The space between the pixels' edges
  G4double fVoidThk; //thickness of the entire detector enclosure
  G4double fEdgeSpace; //derived parameter for edge space
  G4int* DetMaterialList;
  G4Material* zipGe;
  G4Material* zipSi;
  G4bool DrawSolidZipBool;
  G4int towerNumber;

  //set the coordinates for a copy
  void SetCoordinates(G4int copyNo, G4ThreeVector vec);

  //map to coordinates of copies
  mutable std::map<G4int,G4ThreeVector> k100ZipParCoords;
};

// ------------------------------------------------
// ------------------------------------------------

#endif


