// ------------------------------------------------
//
// k100_NaIParametrization.hh
//
// ------------------------------------------------

#ifndef k100_NaIParametrization_H
#define k100_NaIParametrization_H 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalConstants.hh"
#include "G4VPVParameterisation.hh"
#include "G4Box.hh"

// class G4VPhysicalVolume;
// class G4Tubs;
// class G4Material;
// class G4VisAttributes;

// Dummy declarations to get rid of warnings ... (based on the 4.6.2 N02 example)
// class G4Box;
// class G4Trd;
// class G4Trap;
// class G4Cons;
// class G4Orb;
// class G4Sphere;
// class G4Torus;
// class G4Para;
// class G4Hype;
// class G4Polycone;
// class G4Polyhedra;

// ------------------------------------------------
// ------------------------------------------------

class k100_NaIParametrization : public G4VPVParameterisation {
public:
	k100_NaIParametrization( 
		G4double TotalHeight,
		G4double TotalLength,
		G4double TotalWidth,
		G4double TileHeight,
		G4double TileLength,
		G4double TileWidth,
		G4int axis,
		G4int NbTiles);

	virtual ~k100_NaIParametrization();
	virtual void ComputeTransformation (const G4int copyNo, G4VPhysicalVolume* physVol) const;
	virtual void ComputeDimensions(G4Box&,const G4int,const G4VPhysicalVolume*) const;
	//get the coordinates for each copy
  	G4ThreeVector GetCoordinates(G4int copyNo);


private:  // Dummy declarations to get rid of warnings ...

	
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

	G4double    fInterTileGapsX;
	G4double    fInterTileGapsY;
	G4double    fInterTileGapsZ;
	G4int       faxis;
	G4double    fTotalHeight; 
	G4double    fTotalLength; 
	G4double    fTotalWidth; 
	G4double    fTileHeight; 
	G4double    fTileLength; 
	G4double    fTileWidth; 
	//G4Material* fTileMat;
	G4int fNbTiles;
	//G4int nStack;

	//set the coordinates for a copy
	void SetCoordinates(G4int copyNo, G4ThreeVector vec);
	//map to coordinates of copies
  	mutable std::map<G4int,G4ThreeVector> k100NaICoords;
};


// ------------------------------------------------
// ------------------------------------------------

#endif