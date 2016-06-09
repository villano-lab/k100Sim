// ------------------------------------------------
//
// k100_DetectorConstruction.hh : 2016 
//
// ------------------------------------------------

#ifndef k100_DetectorConstruction_H
#define k100_DetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4Region.hh"

class G4VPhysicalVolume;

class k100_DetectorConstructionMessenger;
class k100_ZipSD;
class k100_ZipParameterisation;

  //complicated parameters for complex options
  struct ShieldTest {
    G4double xcntr;
    G4double ycntr;
    G4double zcntr;
    G4double sizel;
    G4double sizew;
    G4double sizethk;
    G4Material* shieldmaterial;
  };

// ------------------------------------------------

class k100_DetectorConstruction : public G4VUserDetectorConstruction
{

public:
  k100_DetectorConstruction();
  ~k100_DetectorConstruction();

public:
  G4VPhysicalVolume* Construct();

  void UpdateGeometry();

  void SetConstructTowerBool(G4bool newVal)      {ConstructTowerBool = newVal;}
  void SetConstructZipBool(G4bool newVal)        {ConstructZipBool = newVal;}
  void SetConstructVetoBool(G4bool newVal)       {ConstructVetoBool = newVal;}
  void SetConstructShieldsBool(G4bool newVal)    {ConstructShieldsBool = newVal;}
  void SetConstructIceBoxBool(G4bool newVal)     {ConstructIceBoxBool = newVal;}
  void SetConstructThermalNeutronBoxBool(G4bool newVal)  {ConstructThermalNeutronBoxBool = newVal&&ConstructZipBool;} //requires construction of Zips
  void SetConstructShieldTestEnvironmentBool(G4bool newVal)      {ConstructShieldTestEnvironmentBool = newVal;}
  void SetConstructShieldTestEnvironmentPos(G4double xcntr,G4double ycntr,G4double zcntr);
  void SetConstructShieldTestEnvironmentSize(G4double sizel,G4double sizew,G4double sizethk);
  void SetConstructShieldTestEnvironmentMat(G4String mat);
  void SetNbOfTowers(G4int newVal)               {NbOfTowers = newVal;}
  void SetNbOfZips(G4int newVal)                 {NbOfZips = newVal;}


  G4bool GetConstructTowerBool()   {return ConstructTowerBool;}
  G4bool GetConstructZipBool()     {return ConstructZipBool;}
  G4bool GetConstructVetoBool()    {return ConstructVetoBool;}
  G4bool GetConstructShieldsBool() {return ConstructShieldsBool;}
  G4bool GetConstructIceBoxBool()  {return ConstructIceBoxBool;}
  G4int GetNbOfTowers()           {return NbOfTowers;}
  G4int GetNbOfZips()             {return NbOfZips;}

  G4bool GetConstructGenericGeometryBool()      {return ConstructGenericGeometryBool;}
  G4bool GetConstructThermalNeutronBoxBool()      {return ConstructThermalNeutronBoxBool;}
  G4bool GetConstructShieldTestEnvironmentBool()      {return ConstructShieldTestEnvironmentBool;}
  struct ShieldTest GetConstructShieldTestEnvironmentParams() {return shieldTestParams;}
  G4String GetConstructShieldTestEnvironmentMat();
  G4int GetConstructGenericTrackerInt() {return ConstructGenericTrackerInt;}
  G4int GetConstructGenericSensitiveInt() {return ConstructGenericSensitiveInt;}
  std::map<G4String,G4int> *GetSensitiveList() { return &k100CollName;}
  G4int    GetNSensitive() {return k100CollName.size();}



  void SetDrawSolidDetBox(G4bool newVal)         {DrawSolidDetBox = newVal;}
  void SetDrawSolidZipBool(G4bool newVal)        {DrawSolidZipBool = newVal;}
  void SetDrawSolidTowerBool(G4bool newVal)      {DrawSolidTowerBool = newVal;}
  void SetDrawSolidVetoBool(G4bool newVal)       {DrawSolidVetoBool = newVal;}
  void SetDrawSolidShieldsBool(G4bool newVal)    {DrawSolidShieldsBool = newVal;}
  void SetDrawSolidIceBoxBool(G4bool newVal)     {DrawSolidIceBoxBool = newVal;}

private:

  G4double world_x, world_y, world_z;
  G4VPhysicalVolume* physicalWorld;

  G4bool DrawSolidDetBox;
  G4bool DrawSolidZipBool, DrawSolidTowerBool, DrawSolidVetoBool, DrawSolidShieldsBool, DrawSolidIceBoxBool;

  G4bool ConstructExperimentBool, ConstructTowerBool, ConstructZipBool;
  G4bool ConstructVetoBool, ConstructShieldsBool, ConstructIceBoxBool;

  G4Material *zipGeMat, *zipSiMat, *towerMat, *scintMat;
  G4Material* polyMat, *mumetalMat;
  G4Material* shieldCuMat, *shieldPbMat;
  G4Material* iceboxCuMat;
  G4Material* defaultMat;
  G4Material* aluminum, *steel, *brass, *helium, *super;

  //SD map
  std::map<G4String,G4int> k100CollName;

  k100_DetectorConstructionMessenger* detectorMessenger;  //pointer to the Messenger
  //  k100_ZipSD*   azipSD;

  G4int NbOfZips, NbOfTowers, NbZipsPerTower;

  G4Region* DetectorRegion;

  G4bool ConstructGenericGeometryBool;
  G4bool ConstructThermalNeutronBoxBool;
  G4bool ConstructShieldTestEnvironmentBool;
  G4int ConstructGenericTrackerInt;
  G4int ConstructGenericSensitiveInt;
  G4int DesignNo;

  //keep the zip parameterization around
  k100_ZipParameterisation *zipParam;

  //structs for grouped parameters
  ShieldTest shieldTestParams;

  void DefineMaterials();
  void ConstructDetector();
  void ConstructZip(G4VPhysicalVolume* physicalDetectorBox, G4ThreeVector positionZip, G4int zipNb);
  void ConstructTower(G4VPhysicalVolume* physicalDetectorBox);
  void ConstructTowerGuts(G4VPhysicalVolume* physicalTower);
  void ConstructEverything(G4LogicalVolume*  logicalWorld);
  void ConstructVeto(G4LogicalVolume*  logicalWorld);
  void ConstructShields(G4LogicalVolume*  logicalWorld);
  void ConstructIceBox(G4LogicalVolume*  logicalWorld);
  void ConstructThermalNeutronBox(G4VPhysicalVolume*  world);
  void ConstructShieldTestEnvironment(G4VPhysicalVolume*  world);
  void FillTheTower(G4VPhysicalVolume* physicalTower, G4int towerNb);

#include "k100_DetectorParameterDef.hh"

public:

  // OK, this realy truly sucks, .... but we have to find a way to generate this string array on the fly
  char** DetCollName;
  char** TowCollName;
  char** VetoCollName;
  G4int* DetMaterials;

};

#endif
