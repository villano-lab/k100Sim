// ------------------------------------------------
//
// k100_NaIParametrization.cc 
//
// ------------------------------------------------

#include "G4UnitsTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4NistManager.hh"

#include "G4ExceptionSeverity.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "k100_NaIParametrization.hh"
//#include "k100vars.hh"

// ------------------------------------------------
// ------------------------------------------------
k100_NaIParametrization :: k100_NaIParametrization( 
		G4double TotalHeight,
		G4double TotalLength,
		G4double TotalWidth,
		G4double TileHeight,
		G4double TileLength,
		G4double TileWidth,
		G4int axis,
		G4int NbTiles) {
	fTotalHeight = TotalHeight;	
	fTotalLength = TotalLength;
	fTotalWidth = TotalWidth;
	fTileHeight = TileHeight;
	fTileLength = TileLength;
	fTileWidth = TileWidth;
	faxis = axis;
	fNbTiles = NbTiles;
	


	fInterTileGapsX = 0;
	fInterTileGapsY = 0;
	fInterTileGapsZ = 0;
	G4int copyNo = 0;
	G4double Zposition = 0;
	G4double Xposition = 0;
	G4double Yposition = 0;
	//nStack = 0;

	//Floor block
	//G4double fridgeHalfHeightToBottomPlate = (12.9045+13.25+0.25)*2.54*cm;
	G4double fridgeHalfHeightToBottomPlate = (12.9045+19.254+0.25)*2.54*cm; //modified 1/1/18 to get floor height right
	//G4double distanceToFloorZ = fridge_z+12.9045*2.54*cm - fridgeHalfHeightToBottomPlate - 21.0*2.54*cm;
	//G4double distanceCenterToFloor = fridgeHalfHeightToBottomPlate + 21.0*2.54*cm;
	G4double distanceCenterToFloor = fridgeHalfHeightToBottomPlate + 21.0*2.54*cm -70.86*mm; //compensate for 70.86mm discrepancy in floor distance 1/1/18
	//G4double floorZ = fridge_z+12.9045*2.54*cm - distanceCenterToFloor;
	G4double floorZ = 0.0+12.9045*2.54*cm - distanceCenterToFloor;

	switch(faxis) {
		case 1: 
			fInterTileGapsX = 9.;//(fTotalWidth - (NbTiles * fTileWidth))/NbTiles;
			Zposition = 0.; 
			break;
		case 3: 
			fInterTileGapsZ = 1.4*mm; // Gap between two tiles in a stack
			//Zposition = -1 * (1255.0/2) + fTileHeight/2; // place first tile of first stack on bottom
			Zposition = (1255.0/2) - fTileHeight/2;      // place first tile of first stack on Top
			break;
		default:
			break;
	}
	if(faxis == 1) {
		if(fTotalWidth*3 < (fTileWidth*NbTiles)) {
			std::cout<<" fTotalWidth : (fTileWidth*NbTiles) :: "<<fTotalWidth<<" : "<<(fTileWidth*NbTiles)<<std::endl;
			G4Exception("Constructor","k100_NaIParametrization",FatalException,"Stack width is more than total");
		}
	}
	else if(faxis == 3) {
		if(fTotalHeight*5 < (fTileHeight*NbTiles)) {
			std::cout<<" fTotalHeight : (fTileHeight*NbTiles) :: "<<fTotalHeight<<" : "<<(fTileHeight*NbTiles)<<std::endl;
			G4Exception("Constructor","k100_NaIParametrization",FatalException,"Stack height is more than total");
		}
	}
	


	G4ThreeVector origin(Xposition,Yposition,Zposition);
	k100NaICoords[copyNo] = origin;
}


k100_NaIParametrization :: ~k100_NaIParametrization() {

}

// ------------------------------------------------
// ------------------------------------------------

void k100_NaIParametrization::ComputeTransformation(const G4int copyNo,
						     G4VPhysicalVolume* physVol) const
{
	G4double Zposition = 0;
	G4double Xposition = 0;
	G4double Yposition = 0;
	G4int nStack;
	//Floor block
	//G4double fridgeHalfHeightToBottomPlate = (12.9045+13.25+0.25)*2.54*cm;
	G4double fridgeHalfHeightToBottomPlate = (12.9045+19.254+0.25)*2.54*cm; //modified 1/1/18 to get floor height right
	//G4double distanceToFloorZ = fridge_z+12.9045*2.54*cm - fridgeHalfHeightToBottomPlate - 21.0*2.54*cm;
	//G4double distanceCenterToFloor = fridgeHalfHeightToBottomPlate + 21.0*2.54*cm;
	G4double distanceCenterToFloor = fridgeHalfHeightToBottomPlate + 21.0*2.54*cm -70.86*mm; //compensate for 70.86mm discrepancy in floor distance 1/1/18
	//G4double floorZ = fridge_z+12.9045*2.54*cm - distanceCenterToFloor;
	G4double floorZ = 0.0+12.9045*2.54*cm - distanceCenterToFloor;
	if(copyNo < 5) nStack = 1;
	else if(copyNo < 10) nStack = 2;
	else if(copyNo < 15) nStack = 3;
	else if(copyNo < 20) nStack = 4;

	if(faxis == 3) {
		
		//Zposition = -1 * (1255.0/2) + fTileHeight/2 + copyNo*(fTileHeight + fInterTileGapsZ) + (nStack-1)*21.75;  // start placing tiles from bottom
		Zposition =  (1255.0/2) - fTileHeight/2 - (copyNo*(fTileHeight + fInterTileGapsZ) + (nStack-1)*21.75);  // start placing tiles from top
		//std::cout<<"copyNo : nStack : Zposition :: "<<copyNo<<" : "<<nStack<<" : "<<Zposition<<std::endl;

	}
	else if(faxis == 1) {
		//Zposition = 0; 
		Xposition = 0. ; // tile placement --> | T1 | T2 | T3 | ---> -Y axis direction
		switch(copyNo) {
			case 0: Xposition = (fTileWidth + fInterTileGapsX); break;
			case 1: Xposition = 0; break;
			case 2: Xposition = -1 * (fTileWidth + fInterTileGapsX); break;
									
			default : G4Exception("ComputeTransformation","k100_NaIParametrization",FatalException,"Can not place more than 3 tiles at the bottom of fridge.");
			break;
		}
		//std::cout<<"copyNo : Xposition :: "<<copyNo<<" : "<<Xposition<<std::endl;
	}
  
  G4ThreeVector origin(Xposition,Yposition,Zposition);
  k100NaICoords[copyNo] = origin;

  
  
  
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);

  //  G4cout << "#### ZipParametrization : ComputeTransformation# Copy Number "  << copyNo << G4endl;
}

// ------------------------------------------------
// ------------------------------------------------

void k100_NaIParametrization::ComputeDimensions(G4Box& NaIDetector, const G4int copyNo,
						 const G4VPhysicalVolume* physVol) const
{
  // All tiles are the same size.

  NaIDetector.SetXHalfLength(fTileWidth/2);
  NaIDetector.SetYHalfLength(fTileLength/2);
  NaIDetector.SetZHalfLength(fTileHeight/2);
  
}

G4Material*  k100_NaIParametrization::ComputeMaterial(const G4int copyNo, G4VPhysicalVolume* physVol, const G4VTouchable* ) 
{
  // Visualization attributes
  G4VisAttributes* VisAttTile = new G4VisAttributes(G4Colour(0/255.,255./255.,0./255.));
  VisAttTile->SetForceSolid(true); 
  G4NistManager* man = G4NistManager::Instance();
  return  man->FindOrBuildMaterial("G4_SODIUM_IODIDE");
 
}

G4ThreeVector k100_NaIParametrization::GetCoordinates(G4int copyNo)
{

   if(k100NaICoords.find(copyNo) != k100NaICoords.end()){
     return k100NaICoords[copyNo];
   }
   else
     return G4ThreeVector(0,0,0);
}
void  k100_NaIParametrization::SetCoordinates(G4int copyNo, G4ThreeVector vec ) 
{
  k100NaICoords[copyNo] = vec; 
}
// ------------------------------------------------
// ------------------------------------------------

