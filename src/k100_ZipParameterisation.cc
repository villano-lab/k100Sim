// ------------------------------------------------
//
// k100_ZipParameterszation.cc 
//
// ------------------------------------------------

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"

#include "k100_ZipParameterisation.hh"

// ------------------------------------------------
// ------------------------------------------------

k100_ZipParameterisation::k100_ZipParameterisation(G4int    NbZipsPerTower,
						   G4double ZSpacing,
						   G4double ZipRadius,
						   G4double ZipDepth,
						   G4int* DetMaterials,
						   G4Material* zipGeMat,
						   G4Material* zipSiMat,
						   G4bool solidZipBool,
						   G4int towerNb)
{
  fNbZipsPerTower = NbZipsPerTower;
  fRadius     = ZipRadius;
  fZSpacing   = ZSpacing;
  fHalfDepth  = ZipDepth*0.5;
  fStartXY    = 0;
  DetMaterialList = DetMaterials;
  zipGe = zipGeMat;
  zipSi = zipSiMat;
  DrawSolidZipBool = solidZipBool;
  towerNumber = towerNb;

  //FIXME kludge to get the shield set up for copyNo zero
  G4int copyNo = 0;
  G4double Zposition = (2.5 - copyNo)*fZSpacing;
  G4double Xposition = 0;
  G4double Yposition = 0;
  G4ThreeVector origin(Xposition,Yposition,Zposition);
  k100ZipParCoords[copyNo] = origin;
  if (fZSpacing < ZipDepth) {
//    G4Exception("k100_ZipParameterisation construction: Thickness>Spacing");
  }

}

// ------------------------------------------------
// ------------------------------------------------

k100_ZipParameterisation::~k100_ZipParameterisation()
{}

// ------------------------------------------------
// ------------------------------------------------

void k100_ZipParameterisation::ComputeTransformation(const G4int copyNo,
						     G4VPhysicalVolume* physVol) const
{
  G4cout << "computing transformation for copyNo: " << copyNo << G4endl;
  G4double Zposition = (2.5 - copyNo)*fZSpacing;
  G4double Xposition = 0;
  G4double Yposition = 0;
  G4ThreeVector origin(Xposition,Yposition,Zposition);
  k100ZipParCoords[copyNo] = origin;
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);

  //  G4cout << "#### ZipParametrization : ComputeTransformation# Copy Number "  << copyNo << G4endl;
}

// ------------------------------------------------
// ------------------------------------------------

void k100_ZipParameterisation::ComputeDimensions(G4Tubs& ZipDetector, const G4int copyNo,
						 const G4VPhysicalVolume* physVol) const
{
  // All Zips are the same size.

  ZipDetector.SetInnerRadius(0);
  ZipDetector.SetOuterRadius(fRadius);
  ZipDetector.SetZHalfLength(fHalfDepth);
  ZipDetector.SetStartPhiAngle(0);
  ZipDetector.SetDeltaPhiAngle(2*pi);
}

// ------------------------------------------------
// ------------------------------------------------

G4Material*  k100_ZipParameterisation::ComputeMaterial(const G4int copyNo, G4VPhysicalVolume* physVol, const G4VTouchable* ) 
{
  // Visualization attributes
  G4VisAttributes* VisAttZipSi = new G4VisAttributes(G4Colour(0/255.,0/255.,128/255.));
  VisAttZipSi->SetForceSolid(DrawSolidZipBool); 

  // Visualization attributes
  G4VisAttributes* VisAttZipGe = new G4VisAttributes(G4Colour(128/255.,0/255.,0/255.));
  VisAttZipGe->SetForceSolid(DrawSolidZipBool); 

  if(DetMaterialList[copyNo + 6*(towerNumber-1)] == 1) {
    physVol->GetLogicalVolume()->SetVisAttributes( VisAttZipGe );
    //    G4cout << "#### ZipParametrization : ComputeMaterial Zip # " << copyNo << " is Ge\n" << G4endl;
    return zipGe;
  } 
  else if(DetMaterialList[copyNo + 6*(towerNumber-1)] == 0) {
    physVol->GetLogicalVolume()->SetVisAttributes( VisAttZipSi );
    //    G4cout << "#### ZipParametrization : ComputeMaterial  Zip # " << copyNo << " is Si\n" << G4endl;
    return zipSi;
  }
  
  return physVol->GetLogicalVolume()->GetMaterial();
}
G4ThreeVector k100_ZipParameterisation::GetCoordinates(G4int copyNo)
{

   G4cout << "got into k100_ZipParameterisation: " << k100ZipParCoords[0].x() << G4endl;
   if(k100ZipParCoords.find(copyNo) != k100ZipParCoords.end()){
     G4cout << "got into FIRST IF k100_ZipParameterisation" << G4endl;
     return k100ZipParCoords[copyNo];
   }
   else
     return G4ThreeVector(0,0,0);
}
void  k100_ZipParameterisation::SetCoordinates(G4int copyNo, G4ThreeVector vec ) 
{
  k100ZipParCoords[copyNo] = vec; 
}

