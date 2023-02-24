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
class k100_StdSD;
class k100_ZipParameterisation;

  //complicated parameters for complex options (Fridge)
  struct Fridge {
    G4bool includeMixture; //not yet implemented
    G4bool pure3HeBath; 
  };

  //complicated parameters for complex options (Frame)
  struct Frame {
    G4bool includeSand; //not yet implemented
  };

  //complicated parameters for complex options (Floor)
  struct Floor {
    G4bool addReBar; //not yet implemented
  };

  //complicated parameters for complex options (Shield)
  struct Shield {
    G4bool addLeadSupports; //not yet implemented
    G4bool HPGeboron; //default false -- use if want remove South wall in favor of HPGe 
    G4bool HPGeboron_wshield; //default false -- use if want have boron shield on 
    G4bool addNaISouth; //add 2 NaI detectors to South shielding wall, make sensitive
    G4bool addBasePoly; //add poly panels on the bottom 
    G4bool addBaseLead; //add lead sheets on the bottom 
    G4int mod; //parameter for doing different modifications
  };

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

  //complicated parameters for complex options (gamma coincidence)
  struct GammaCoin {
    G4double xcntr;
    G4double ycntr;
    G4double zcntr;
    G4double sizer;
    G4double sizethk;
    G4Material* coinmaterial;
  };

  //complicated parameters for complex options (PuBe setup)
  struct PuBeNaICoin {
    G4bool westPolySensitivity; //default false
    G4bool doPuBeGamma; //default true
    G4bool addBarrel;
    G4bool NaIsensitivity; //true if dets are to be sensitive
    G4bool doR66;
    G4bool doR62;
    G4int mod; //parameter for doing different modifications:
               // 0 - no mod
	       // 1 - add poly layer 1.6m down from origin, in pit and lead shield on floor
    G4bool doOrb;
    G4ThreeVector OrbPos;
    G4double OrbRad;
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

  void DetSizeMod(G4double R, G4double thk);

  void SetConstructTowerBool(G4bool newVal)      {ConstructTowerBool = newVal;}
  void SetConstructZipBool(G4bool newVal)        {ConstructZipBool = newVal;}
  void SetConstructVetoBool(G4bool newVal)       {ConstructVetoBool = newVal;}
  void SetConstructShieldsBool(G4bool newVal)    {ConstructShieldsBool = newVal;}
  void SetConstructShields_HPGeboron_wshield(G4bool newVal)    {shieldParams.HPGeboron_wshield = newVal;}
  void SetConstructShields_HPGeboron(G4bool newVal)    {
    shieldParams.HPGeboron = newVal; 
    if(newVal)
      shieldParams.addNaISouth = !newVal;
  }
  void SetConstructShields_addNaISouth(G4bool newVal)    {
    shieldParams.addNaISouth = newVal; 
    if(newVal)
      shieldParams.HPGeboron = !newVal;
  }
  void SetConstructShields_addBasePoly(G4bool newVal)    {shieldParams.addBasePoly = newVal;}
  void SetConstructShields_addBaseLead(G4bool newVal)    {shieldParams.addBaseLead = newVal;}
  void SetConstructShields_mod(G4int newVal)    {shieldParams.mod = newVal;}
  void SetConstructIceBoxBool(G4bool newVal)     {ConstructIceBoxBool = newVal;}
  void SetConstructIceBox_pure3HeBath(G4bool newVal)    {fridgeParams.pure3HeBath = newVal;}
  void SetConstructFloorBool(G4bool newVal)     {ConstructFloorBool = newVal;}
  void SetConstructWallsBool(G4bool newVal)     {ConstructWallsBool = newVal;}
  void SetConstructCeilingBool(G4bool newVal)     {ConstructCeilingBool = newVal;}
  void SetConstructWestReflectorBool(G4bool newVal)     {ConstructWestReflectorBool = newVal;}
  void SetConstructFrameBool(G4bool newVal)     {ConstructFrameBool = newVal;}
  void SetConstructPuBeSourceAndShieldBool(G4bool newVal)     {ConstructPuBeSourceAndShieldBool = newVal;}
  void SetConstructPuBeSourceAndShield_westPolySensitivity(G4bool newVal)    {pubeNaIParams.westPolySensitivity = newVal;}
  void SetConstructPuBeSourceAndShield_doPuBeGamma(G4bool newVal)    {pubeNaIParams.doPuBeGamma = newVal;}
  void SetConstructPuBeSourceAndShield_addBarrel(G4bool newVal)    {pubeNaIParams.addBarrel = newVal;}
  void SetConstructPuBeSourceAndShield_setNaISensitive(G4bool newVal)    {pubeNaIParams.NaIsensitivity = newVal;}
  void SetConstructPuBeSourceAndShield_doR66(G4bool newVal)    {pubeNaIParams.doR66 = newVal;}
  void SetConstructPuBeSourceAndShield_doR62(G4bool newVal)    {pubeNaIParams.doR62 = newVal;}
  void SetConstructPuBeSourceAndShield_doOrb(G4bool newVal)    {pubeNaIParams.doOrb = newVal;}
  void SetConstructPuBeSourceAndShield_mod(G4int newVal)    {pubeNaIParams.mod = newVal;}
  void SetConstructPuBeSourceAndShield_OrbPos(G4ThreeVector newVal)    {pubeNaIParams.OrbPos = newVal;}
  void SetConstructPuBeSourceAndShield_OrbRad(G4double newVal)    {pubeNaIParams.OrbRad = newVal;}
  void SetConstructThermalNeutronBoxBool(G4bool newVal)  {ConstructThermalNeutronBoxBool = newVal&&ConstructZipBool;} //requires construction of Zips
  void SetConstructShieldTestEnvironmentBool(G4bool newVal)      {ConstructShieldTestEnvironmentBool = newVal&&ConstructZipBool;} //requires construction of Zips
  void SetConstructShieldTestEnvironmentPos(G4double xcntr,G4double ycntr,G4double zcntr);
  void SetConstructShieldTestEnvironmentSize(G4double sizel,G4double sizew,G4double sizethk);
  void SetConstructShieldTestEnvironmentMat(G4String mat);
  void SetConstructSimpleGammaCoinBool(G4bool newVal)  {ConstructSimpleGammaCoinBool = newVal&&ConstructZipBool;} //requires construction of Zips
  void SetConstructSimpleGammaCoinPos(G4double xcntr,G4double ycntr,G4double zcntr);
  void SetConstructSimpleGammaCoinSize(G4double sizer, G4double sizethk);
  void SetConstructSimpleGammaCoinMat(G4String mat);
  void SetConstructPuBeNaIBool(G4bool newVal)  {ConstructPuBeNaIBool = newVal&&ConstructZipBool;} //requires construction of Zips
  void SetFirstDetGe(G4bool newVal)  {FirstDetGe = newVal;}
  void SetNbOfTowers(G4int newVal)               {NbOfTowers = newVal;}
  void SetNbOfZips(G4int newVal)                 {NbOfZips = newVal;}
  //spandey
  void SetConstructNaIBool(G4bool newVal)     {ConstructNaIArrayBool = newVal;}
  void SetConstructBoronShieldBool(G4bool newVal)     {ConstructBoronShieldBool = newVal;}
  void SetSodiumBorateDensityFraction(G4double newVal) {SodiumBorateDensityFraction = newVal;}
  void SetNbBoronShieldVert(G4int newVal) {NbBoronShieldVert = newVal;}
  void SetNbBoronShieldHori(G4int newVal) {NbBoronShieldHori = newVal;}
  void SetBoronShieldThickness(G4double newVal) { BoronShieldThickness = newVal;}
  void SetConstructPolyBox(G4bool newVal) {ConstructPolyBoxBool = newVal;}


  G4bool GetConstructTowerBool()   {return ConstructTowerBool;}
  G4bool GetConstructZipBool()     {return ConstructZipBool;}
  G4bool GetConstructVetoBool()    {return ConstructVetoBool;}
  G4bool GetConstructShieldsBool() {return ConstructShieldsBool;}
  G4bool GetConstructShields_HPGeboron_wshield()    {return shieldParams.HPGeboron_wshield;}
  G4bool GetConstructShields_HPGeboron()    {return shieldParams.HPGeboron;}
  G4bool GetConstructShields_addNaISouth()    {return shieldParams.addNaISouth;}
  G4bool GetConstructShields_addBasePoly()    {return shieldParams.addBasePoly;}
  G4bool GetConstructShields_addBaseLead()    {return shieldParams.addBaseLead;}
  G4int  GetConstructShields_mod()    {return shieldParams.mod;}
  G4bool GetConstructPuBeSourceAndShield_westPolySensitivity()    {return pubeNaIParams.westPolySensitivity;}
  G4bool GetConstructPuBeSourceAndShield_doPuBeGamma()    {return pubeNaIParams.doPuBeGamma;}
  G4bool GetConstructPuBeSourceAndShield_addBarrel()    {return pubeNaIParams.addBarrel;}
  G4bool GetConstructPuBeSourceAndShield_doR66()    {return pubeNaIParams.doR66;}
  G4bool GetConstructPuBeSourceAndShield_doR62()    {return pubeNaIParams.doR62;}
  G4int  GetConstructPuBeSourceAndShield_mod()    {return pubeNaIParams.mod;}
  G4bool GetConstructIceBoxBool()  {return ConstructIceBoxBool;}
  G4bool GetConstructFloorBool()  {return ConstructFloorBool;}
  G4bool GetConstructWallsBool()  {return ConstructWallsBool;}
  G4bool GetConstructCeilingBool()  {return ConstructCeilingBool;}
  G4bool GetConstructWestReflectorBool()  {return ConstructWestReflectorBool;}
  G4bool GetConstructFrameBool()  {return ConstructFrameBool;}
  G4bool GetConstructPuBeSourceAndShieldBool()  {return ConstructPuBeSourceAndShieldBool;}
  G4int GetNbOfTowers()           {return NbOfTowers;}
  G4int GetNbOfZips()             {return NbOfZips;}

  G4int GetNbBoronShieldVert () {return NbBoronShieldVert;}
  G4int GetNbBoronShieldHori () {return NbBoronShieldHori;}
  G4bool GetConstructPolyBox() {return ConstructPolyBoxBool;}


  G4bool GetConstructGenericGeometryBool()      {return ConstructGenericGeometryBool;}
  G4bool GetConstructThermalNeutronBoxBool()      {return ConstructThermalNeutronBoxBool;}
  G4bool GetConstructShieldTestEnvironmentBool()      {return ConstructShieldTestEnvironmentBool;}
  struct ShieldTest GetConstructShieldTestEnvironmentParams() {return shieldTestParams;}
  G4String GetConstructShieldTestEnvironmentMat();
  G4bool GetConstructSimpleGammaCoinBool()      {return ConstructSimpleGammaCoinBool;}
  struct GammaCoin GetConstructSimpleGammaCoinParams() {return gammaCoinParams;}
  G4bool GetConstructPuBeNaIBool()      {return ConstructPuBeNaIBool;}
  struct PuBeNaICoin GetConstructPuBeNaIParams() {return pubeNaIParams;}
  struct Fridge GetConstructFridgeParams() {return fridgeParams;}
  G4String GetConstructSimpleGammaCoinMat();
  G4int GetConstructGenericTrackerInt() {return ConstructGenericTrackerInt;}
  G4int GetConstructGenericSensitiveInt() {return ConstructGenericSensitiveInt;}
  std::map<G4String,G4int> *GetSensitiveList() { return &k100CollName;}
  G4int    GetNSensitive() {return k100CollName.size();}

  G4double GetSodiumBorateDensityFraction() {return SodiumBorateDensityFraction;}
  G4double GetBoronShieldThickness() {return BoronShieldThickness;}


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
  G4bool ConstructFloorBool;
  G4bool ConstructWallsBool;
  G4bool ConstructCeilingBool;
  G4bool ConstructWestReflectorBool;
  G4bool ConstructFrameBool;
  G4bool ConstructPuBeSourceAndShieldBool;
  //spandey
  G4bool ConstructNaIArrayBool;
  G4bool ConstructBoronShieldBool;
  G4double SodiumBorateDensityFraction;
  G4double BoronShieldThickness;
  G4bool ConstructPolyBoxBool;

  G4Material *zipGeMat, *zipSiMat, *towerMat, *scintMat;
  G4Material* polyMat, *mumetalMat;
  G4Material* d2oMat, *h2oMat, *naiMat;
  G4Material* sio2Mat;
  G4Material* shieldCuMat, *shieldPbMat;
  G4Material* iceboxCuMat;
  G4Material* defaultMat;
  G4Material* aluminum, *steel, *brass, *helium, *super;
  G4Material* stillHe,*MCHe;
  G4Material* blastsand;
  G4Material* carbonsteel;
  G4Material* lightaluminum;
  G4Material* wood;
  G4Material* sodium_borate_anhydrous; //https://www.fishersci.com/shop/products/sodium-tetraborate-anhydrous-99-5-metals-basis-alfa-aesar-2/AA12305A7
  G4Material* G4NISTconcrete,*G4NISTair,*G4NISTNaI,*G4NISTPVC,*G4NISTPE,*G4NISTlucite,*G4NISTparaffin,*G4NISTstainless;
  G4Material* G4NISTAl;
  G4Material* G4NISTGypsum; //drywall
  G4Material* boronShieldMat; // sodium tetraborate decahydrate : https://www.sigmaaldrich.com/US/en/product/sigald/s9640
  //SD map
  std::map<G4String,G4int> k100CollName;
  std::map<G4String,k100_ZipSD*> k100CollPoint;
  std::map<G4String,k100_StdSD*> k100CollPointStd;

  k100_DetectorConstructionMessenger* detectorMessenger;  //pointer to the Messenger
  //  k100_ZipSD*   azipSD;

  G4int NbOfZips, NbOfTowers, NbZipsPerTower;

  G4int NbBoronShieldVert;
  G4int NbBoronShieldHori;

  G4Region* DetectorRegion;


  G4bool ConstructGenericGeometryBool;
  G4bool ConstructThermalNeutronBoxBool;
  G4bool ConstructShieldTestEnvironmentBool;
  G4bool ConstructSimpleGammaCoinBool;
  G4bool ConstructPuBeNaIBool;
  G4bool FirstDetGe;
  G4int ConstructGenericTrackerInt;
  G4int ConstructGenericSensitiveInt;
  G4int DesignNo;

  //keep the zip parameterization around
  k100_ZipParameterisation *zipParam;

  //structs for grouped parameters
  ShieldTest shieldTestParams;
  GammaCoin  gammaCoinParams;
  PuBeNaICoin  pubeNaIParams;
  Fridge	fridgeParams;
  Frame		frameParams;
  Floor		floorParams;
  Shield	shieldParams;

  void DefineMaterials();
  void ConstructDetector();
  void ConstructZip(G4VPhysicalVolume* physicalDetectorBox, G4ThreeVector positionZip, G4int zipNb);
  void ConstructTower(G4VPhysicalVolume* physicalDetectorBox);
  void ConstructTowerGuts(G4VPhysicalVolume* physicalTower);
  void ConstructEverything(G4LogicalVolume*  logicalWorld);
  void ConstructVeto(G4LogicalVolume*  logicalWorld);
  void ConstructShields(G4LogicalVolume*  logicalWorld);
  void ConstructIceBox(G4LogicalVolume*  logicalWorld);
  void ConstructFloor(G4VPhysicalVolume*  world);
  void ConstructWalls(G4VPhysicalVolume*  world);
  void ConstructCeiling(G4VPhysicalVolume*  world);
  void ConstructWestReflector(G4VPhysicalVolume*  world);
  void ConstructFrame(G4VPhysicalVolume*  world);
  void ConstructPuBeSourceAndShield(G4VPhysicalVolume*  world);
  void ConstructThermalNeutronBox(G4VPhysicalVolume*  world);
  void ConstructShieldTestEnvironment(G4VPhysicalVolume*  world);
  void ConstructSimpleGammaCoin(G4VPhysicalVolume*  world);
  void ConstructPuBeNaI(G4VPhysicalVolume*  world);
  void FillTheTower(G4VPhysicalVolume* physicalTower, G4int towerNb);
  //spandey
  void ConstructNaIArray(G4LogicalVolume*  logicalWorld);
#include "k100_DetectorParameterDef.hh"

public:

  // OK, this realy truly sucks, .... but we have to find a way to generate this string array on the fly
  char** DetCollName;
  char** TowCollName;
  char** VetoCollName;
  G4int* DetMaterials;

};

#endif
