// ------------------------------------------------
//
// k100_DetectorConstruction.cc : 2016 
//
// ------------------------------------------------

#include "G4UnitsTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"

#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4ThreeVector.hh"

#include "G4AssemblyVolume.hh"


#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "k100vars.hh"

#include "k100_DetectorConstruction.hh"
#include "k100_DetectorConstructionMessenger.hh"
#include "k100_ZipParameterisation.hh"
#include "k100_ZipSD.hh"
#include "k100_VetoSD.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

// ------------------------------------------------

k100_DetectorConstruction::k100_DetectorConstruction()
{

  // -------- The World ---------
  //world_x = 250.*cm; world_y = 250.*cm; world_z = 250.*cm;
  world_x = 250.*cm; world_y = 250.*cm; world_z = 500.*cm;
  
#include "k100_DetectorParameterDef.icc"
  
  // Create commands for interactive definition of the detector
  detectorMessenger = new k100_DetectorConstructionMessenger(this);
  
  //
  NbOfTowers = 1;
  NbZipsPerTower = 1; // Corrected: Originally 2
  NbOfZips = NbOfTowers * NbZipsPerTower;
  //
  ConstructExperimentBool = true;
  ConstructTowerBool = true;
  ConstructZipBool = true;
  ConstructVetoBool = false;
  ConstructShieldsBool = true;
  ConstructIceBoxBool = true;
  ConstructThermalNeutronBoxBool = false;
  SetConstructShieldTestEnvironmentBool(false); //note, requires construct ZIP bool

  //
  DrawSolidDetBox = true; DrawSolidZipBool = true;
  DrawSolidTowerBool = false; 
  DrawSolidVetoBool = false;
  DrawSolidShieldsBool = false; DrawSolidIceBoxBool = false;


  // ---------Material Definition--------------
  DefineMaterials();

  //more complicated parameters
  shieldTestParams.xcntr = 100.0*cm;
  shieldTestParams.ycntr = 0.0*cm;
  shieldTestParams.zcntr = 100.0*cm;
  shieldTestParams.sizel = 10.0*cm;
  shieldTestParams.sizew = 10.0*cm;
  shieldTestParams.sizethk = 10.0*cm;
  shieldTestParams.shieldmaterial = polyMat; //MUST be after DefineMaterials()


  // ---------Detector Names--------------
  DetCollName = new char*[30];  TowCollName = new char*[5];     DetMaterials = new G4int [30];
  DetCollName[0]  = "zip01";  DetMaterials[0]  = 1; //Ge  ///Ge = 1 Si = 0
  DetCollName[1]  = "zip02";  DetMaterials[1]  = 1; //Ge
  DetCollName[2]  = "zip03";  DetMaterials[2]  = 1; //Ge
  DetCollName[3]  = "zip04";  DetMaterials[3]  = 0; //Si
  DetCollName[4]  = "zip05";  DetMaterials[4]  = 1; //Ge
  DetCollName[5]  = "zip06";  DetMaterials[5]  = 0; //Si
  
  DetCollName[6]  = "zip07";  DetMaterials[6]  = 0; //Si
  DetCollName[7]  = "zip08";  DetMaterials[7]  = 0; //Si
  DetCollName[8]  = "zip09";  DetMaterials[8]  = 1; //Ge
  DetCollName[9]  = "zip10";  DetMaterials[9]  = 0; //Si
  DetCollName[10] = "zip11";  DetMaterials[10] = 1; //Ge
  DetCollName[11] = "zip12";  DetMaterials[11] = 0; //Si 
  
  DetCollName[12] = "zip13";  DetMaterials[12] = 0; //Si
  DetCollName[13] = "zip14";  DetMaterials[13] = 1; //Ge
  DetCollName[14] = "zip15";  DetMaterials[14] = 0; //Si
  DetCollName[15] = "zip16";  DetMaterials[15] = 1; //Ge
  DetCollName[16] = "zip17";  DetMaterials[16] = 1; //Ge
  DetCollName[17] = "zip18";  DetMaterials[17] = 1; //Ge 
  
  DetCollName[18] = "zip19";  DetMaterials[18] = 0; //Si
  DetCollName[19] = "zip20";  DetMaterials[19] = 1; //Ge
  DetCollName[20] = "zip21";  DetMaterials[20] = 0; //Si
  DetCollName[21] = "zip22";  DetMaterials[21] = 1; //Ge
  DetCollName[22] = "zip23";  DetMaterials[22] = 1; //Ge
  DetCollName[23] = "zip24";  DetMaterials[23] = 1; //Ge 

  DetCollName[24] = "zip25";  DetMaterials[24] = 1; //Ge
  DetCollName[25] = "zip26";  DetMaterials[25] = 1; //Ge
  DetCollName[26] = "zip27";  DetMaterials[26] = 0; //Si
  DetCollName[27] = "zip28";  DetMaterials[27] = 1; //Ge
  DetCollName[28] = "zip29";  DetMaterials[28] = 1; //Ge
  DetCollName[29] = "zip30";  DetMaterials[29] = 1; //Ge 

  TowCollName[0]  = "tower1"; 
  TowCollName[1]  = "tower2";
  TowCollName[2]  = "tower3";
  TowCollName[3]  = "tower4";
  TowCollName[4]  = "tower5";

  VetoCollName = new char*[19];
  VetoCollName[1-1]  = "Veto1";   VetoCollName[8-1]   = "Veto8";    VetoCollName[15-1]  = "Veto15"; 
  VetoCollName[2-1]  = "Veto2";   VetoCollName[9-1]   = "Veto9";    VetoCollName[16-1]  = "Veto16"; 
  VetoCollName[3-1]  = "Veto3";   VetoCollName[10-1]  = "Veto10";   VetoCollName[17-1]  = "Veto17"; 
  VetoCollName[4-1]  = "Veto4";   VetoCollName[11-1]  = "Veto11";   VetoCollName[18-1]  = "Veto18"; 
  VetoCollName[5-1]  = "Veto5";   VetoCollName[12-1]  = "Veto12";   VetoCollName[19-1]  = "Veto19"; 
  VetoCollName[6-1]  = "Veto6";   VetoCollName[13-1]  = "Veto13";   
  VetoCollName[7-1]  = "Veto7";   VetoCollName[14-1]  = "Veto14";   


//   //// Maybe some code along these lines 

//   char tmpchar [5]; 
//   char** DetCollName = new char*[12];
  
//   for(G4int ii=0; ii<NumDetCollName; ii++){
//     sprintf(tmpchar , "zip_%.2d" , ii);
//       G4cout << "\n>>> (BeginofRunAction) 2  tmpchar " << tmpchar << G4endl;
//       sprintf(DetCollName[ii] , "%s", tmpchar);
//       G4cout << "\n>>> (BeginofRunAction) 2  DetCollName[ii] " << DetCollName[ii] << G4endl;
//   }
  
//   //// Sadly this does not work

}

// ------------------------------------------------


k100_DetectorConstruction::~k100_DetectorConstruction()
{
  //delete the messenger
  delete detectorMessenger;

  // Clean old geometry, if any
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
}

// ------------------------------------------------

void k100_DetectorConstruction::DefineMaterials()
{

  G4String name, symbol;    
  G4double a, z, density;

  G4int ncomponents, natoms;
  //G4double fractionmass;
  G4double temperature, pressure;


  // ------------------------------------------------
  // Define Elements 
  // ------------------------------------------------
  
  // Define Aluminum
  G4Element* elementAl=new G4Element(name="Aluminum", symbol="Al", z=13., a=26.98*g/mole);

  // Define Hydrogen
  G4Element* elementH=new G4Element(name="Hydrogen", symbol="H", z=1., a=1.01*g/mole);

  // Define Helium
  G4Element* elementHe = new G4Element(name="Helium", symbol="He", z=2., a=4.003*g/mole);

  // Define Carbon
  G4Element* elementC=new G4Element(name="Carbon", symbol="C", z=6., a=12.011*g/mole);

  // Define Silicon
  G4Element* elementSi = new G4Element(name="Silicon", symbol="Si", z=14., a=28.09*g/mole);

  // Define Iron
  G4Element* elementFe = new G4Element(name="Iron", symbol="Fe", z=26., a=55.845*g/mole);

  // Define Copper
  G4Element* elementCu = new G4Element(name="Copper", symbol="Cu", z=29., a=63.5460*g/mole);

  // Define Germanium   
  G4Element* elementGe = new G4Element(name="Germanium", symbol="Ge", z=32., a=72.61*g/mole);

  // Define Lead
  G4Element* elementPb = new G4Element(name="Lead",symbol="Pb", z=82., a=207.2*g/mole);
  
  // Define Phosphorus
  G4Element* elementP = new G4Element(name="Phosphorus",symbol="P", z=15., a=30.97*g/mole);
  
  // Define Sulfur
  G4Element* elementS = new G4Element(name="Sulfur",symbol="S", z=16., a=32.065*g/mole);
  
  // Define Nickel
  G4Element* elementNi = new G4Element(name="Nickel",symbol="Ni", z=28., a=58.6934*g/mole);

  // Define Chromium
  G4Element* elementCr = new G4Element(name="Chromium",symbol="Cr", z=24., a=51.9961*g/mole);

  // Define Manganese
  G4Element* elementMn = new G4Element(name="Manganese",symbol="Mn", z=25., a=54.9380*g/mole);

  // Define Zinc
  G4Element* elementZn = new G4Element(name="Zinc",symbol="Zn", z=30., a=65.38*g/mole);

  // ------------------------------------------------

  // ------------------------------------------------
  // Define Simple Materials
  // ------------------------------------------------

  //   // ------------------------------------------------
  //   // Define a material from elements.   Case 1 : chemical molecule
  //   // ------------------------------------------------

  // Aluminum
  G4Material* Aluminum = new G4Material(name="Aluminum", density = 2.7*g/cm3, ncomponents=1);
  Aluminum->AddElement(elementAl, natoms=1);

  // Silicon 
  G4Material* Silicon = new G4Material(name="Silicon", density = 2.330*g/cm3, ncomponents=1);
  Silicon->AddElement(elementSi, natoms=1);
  
  // Germanium 
  G4Material* Germanium = new G4Material(name="Germanium", density = 5.323*g/cm3, ncomponents=1);
  Germanium->AddElement(elementGe, natoms=1);

  // Copper
  G4Material* Copper = new G4Material(name="Copper", density = 8.920*g/cm3, ncomponents=1);
  Copper->AddElement(elementCu, natoms=1);

  // Lead
  G4Material* Lead = new G4Material(name="Lead_brick", density = 11.35*g/cm3, ncomponents=1);
  Lead->AddElement(elementPb,1);      

  // Sulfur
  G4Material* Sulfur = new G4Material(name="Sulfur", density = 2.07*g/cm3, ncomponents=1);
  Sulfur->AddElement(elementS, natoms=1);

  // Phosphorus
  G4Material* Phosphorus = new G4Material(name="Phosphorus", density = 1.823*g/cm3, ncomponents = 1);
  Phosphorus->AddElement(elementP, natoms=1);

  // Nickel
  G4Material* Nickel = new G4Material(name="Nickel", density = 8.908*g/cm3, ncomponents = 1);
  Nickel->AddElement(elementNi, natoms=1);

  // Chromium
  G4Material* Chromium = new G4Material(name="Chromium", density = 7.19*g/cm3, ncomponents=1);
  Chromium->AddElement(elementCr, natoms=1);

  // Manganese
  G4Material* Manganese = new G4Material(name="Manganese", density=7.21*g/cm3, ncomponents=1);
  Manganese->AddElement(elementMn, natoms=1);

  // Carbon
  G4Material* Carbon = new G4Material(name="Carbon", density=2.267*g/cm3, ncomponents=1);
  Carbon->AddElement(elementC, natoms=1);

  // Iron
  G4Material* Iron = new G4Material(name="Iron", density = 7.874*g/cm3, ncomponents=1);
  Iron->AddElement(elementFe, natoms=1);

  // Zinc
  G4Material* Zinc = new G4Material(name="Zinc", density = 7.14*g/cm3, ncomponents=1);
  Zinc->AddElement(elementZn, natoms=1);

  // Scintillator
  G4Material* Scint =new G4Material(name="Scintillator",density = 1.032*g/cm3, ncomponents=2) ;
  Scint->AddElement(elementH, natoms=11);
  Scint->AddElement(elementC, natoms=10);

  // Poly
  G4Material* poly=new G4Material(name="Poly", density = 0.935*g/cm3, ncomponents=2);
  poly->AddElement(elementH,natoms=2);
  poly->AddElement(elementC,natoms=1);

  // MuMetal
  G4Material* mumetal=new G4Material(name="MuMetal", density=1*g/cm3, ncomponents=2);
  mumetal->AddElement(elementFe, .19);
  mumetal->AddElement(elementNi, .81);

  // Brass (high brass alloy)
  G4Material* Brass=new G4Material(name="Brass",density=8.2*g/cm3,ncomponents=2);
  Brass->AddMaterial(Copper, .65);
  Brass->AddMaterial(Zinc, .35);

  // Liquid Helium
  G4Material* liquidHelium=new G4Material(name="liquidHelium", density=.1412*g/cm3, ncomponents=1);
  liquidHelium->AddElement(elementHe,natoms=2);

  // Stainless Steel (Alloy 304)
  G4Material* Steel = new G4Material(name="Steel", density=8.03*g/cm3, ncomponents=8);
  Steel->AddMaterial(Carbon, .0008);
  Steel->AddMaterial(Manganese, .02);
  Steel->AddMaterial(Chromium, .19);
  Steel->AddMaterial(Nickel, .0925);
  Steel->AddMaterial(Silicon, .01);
  Steel->AddMaterial(Phosphorus, .00045);
  Steel->AddMaterial(Sulfur, .0003);
  Steel->AddMaterial(Iron, .68595);

  // Superinsulation/Vacuum Layer (Modeled as Superinsulation at half normal density.  Normal Density 1.3890*g/cm3)
  G4Material* Super = new G4Material(name="Super", density=.6945*g/cm3, ncomponents=2);
  Super->AddMaterial(Aluminum, .5);
  Super->AddMaterial(poly, .5);

  // Vacuum
  density     = universe_mean_density;   //from PhysicalConstants.h
  pressure    = 1.0e-19*pascal; 
  temperature = 2.73*kelvin;
  G4Material* Vacuum = new G4Material("Vacuum", z=1.0, a=1.01*g/mole, density,
				      kStateGas, temperature, pressure);

  // ------------------------------------------------

  // ------------------------------------------------
  // PROPERTIES OF SCINTILLATOR
  // ------------------------------------------------

  const G4int NUMENTRIES = 32;

  G4double PPCKOV[NUMENTRIES] =
    { 2.034E-9*GeV, 2.068E-9*GeV, 2.103E-9*GeV, 2.139E-9*GeV,
      2.177E-9*GeV, 2.216E-9*GeV, 2.256E-9*GeV, 2.298E-9*GeV,
      2.341E-9*GeV, 2.386E-9*GeV, 2.433E-9*GeV, 2.481E-9*GeV,
      2.532E-9*GeV, 2.585E-9*GeV, 2.640E-9*GeV, 2.697E-9*GeV,
      2.757E-9*GeV, 2.820E-9*GeV, 2.885E-9*GeV, 2.954E-9*GeV,
      3.026E-9*GeV, 3.102E-9*GeV, 3.181E-9*GeV, 3.265E-9*GeV,
      3.353E-9*GeV, 3.446E-9*GeV, 3.545E-9*GeV, 3.649E-9*GeV,
      3.760E-9*GeV, 3.877E-9*GeV, 4.002E-9*GeV, 4.136E-9*GeV };

  G4double RINDEX1[NUMENTRIES] =
    { 
      1.58, 1.58, 1.58, 1.58, 1.58,
      1.58, 1.58, 1.58, 1.58, 1.58,
      1.58, 1.58, 1.58, 1.58, 1.58,
      1.58, 1.58, 1.58, 1.58, 1.58,
      1.58, 1.58, 1.58, 1.58, 1.58,
      1.58, 1.58, 1.58, 1.58, 1.58,
      1.58, 1.58
    } ;

  G4double RINDEX2[NUMENTRIES] =
    { 
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00 
    };

  G4double ABSORPTION1[NUMENTRIES] =
    {
      344.8*cm,  408.2*cm,  632.9*cm,  917.4*cm, 1234.6*cm, 1388.9*cm,
      1515.2*cm, 1724.1*cm, 1886.8*cm, 2000.0*cm, 2631.6*cm, 3571.4*cm,
      4545.5*cm, 4761.9*cm, 5263.2*cm, 5263.2*cm, 5555.6*cm, 5263.2*cm,
      5263.2*cm, 4761.9*cm, 4545.5*cm, 4166.7*cm, 3703.7*cm, 3333.3*cm,
      3000.0*cm, 2850.0*cm, 2700.0*cm, 2450.0*cm, 2200.0*cm, 1950.0*cm,
      1750.0*cm, 1450.0*cm 
    };

  G4double SCINTILLATION[NUMENTRIES] =
    { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00 };

  G4MaterialPropertiesTable* Scint_MPT = 
    new G4MaterialPropertiesTable();

  Scint_MPT->AddProperty("RINDEX", PPCKOV, RINDEX1, NUMENTRIES);
  Scint_MPT->AddProperty("ABSLENGTH",PPCKOV, ABSORPTION1, NUMENTRIES);
  Scint_MPT->AddProperty("SCINTILLATION",PPCKOV, SCINTILLATION, NUMENTRIES);

  // ------------------------------------------------

  // Assign materials to some of the major components

  defaultMat = Vacuum;
  zipGeMat = Germanium;
  zipSiMat = Silicon;
  towerMat = Copper;
  scintMat = Scint;
  polyMat = poly;
  shieldCuMat = Copper;
  shieldPbMat = Lead;
  iceboxCuMat = Copper;
  mumetalMat = mumetal;
  aluminum=Aluminum;
  steel=Steel;
  brass=Brass;
  helium=liquidHelium;
  super=Super;
  // ------------------------------------------------

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

} // ends Material Definitions

// ------------------------------------------------
// ------------------------------------------------

void k100_DetectorConstruction::UpdateGeometry()
{

  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

} // ends UpdateGeometry

// ------------------------------------------------
// ------------------------------------------------

G4VPhysicalVolume* k100_DetectorConstruction::Construct()
{

  // Clean old geometry, if any
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // Prepare to declare sensitive detectors
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  //deactivate the sensitive detectors
  std::map<G4String,G4int>::iterator it;
  for(it=k100CollName.begin();it!=k100CollName.end();++it){
    SDman->Activate(it->first,false);
  }

  // ------------ Construct the Physical world ---------------

  // Construct the World
  G4Box* solidWorld = new G4Box("world_S", world_x, world_y, world_z);
  G4LogicalVolume*  logicalWorld = new G4LogicalVolume(solidWorld,  // The solid
						       defaultMat, // Material
						       "world_L",  // Name
						       0,0,0);
  // Physical volume
  physicalWorld = new G4PVPlacement(0,               // no rotation
				    G4ThreeVector(), // at (0, 0, 0)
				    "world_P",       // name (using the second constructor)
				    logicalWorld,    // the logical volume to use
				    NULL,            // the mother volume
				    false,           // no boolean operation
				    0);              // copy number
  // Visualization attributes
  G4VisAttributes* VisAttWorld = new G4VisAttributes(G4Colour(204/255.,255/255.,255/255.));
  logicalWorld->SetVisAttributes(VisAttWorld);
  // Make Invisible
  logicalWorld->SetVisAttributes(G4VisAttributes::Invisible);

  // --------- End Construct the Physical world --------------

  // ------------ Construct the Whole shebang ---------------
  ConstructEverything(logicalWorld);
  if(ConstructTowerBool) {ConstructTower(physicalWorld);}
  //  if(ConstructExperimentBool)  {ConstructDl();}
  // --------- End Construct the whole shebang --------------
  
  //FIXME construct stuff that can't be constructed before tower
  if(ConstructShieldTestEnvironmentBool)  {ConstructShieldTestEnvironment(physicalWorld);}


  return physicalWorld;
} // ends Construct (for physical world, trigger for others)

// ------------------------------------------------
// ------------------------------------------------
void k100_DetectorConstruction::SetConstructShieldTestEnvironmentPos(G4double xcntr,G4double ycntr,G4double zcntr)
{
  shieldTestParams.xcntr = xcntr;
  shieldTestParams.ycntr = ycntr;
  shieldTestParams.zcntr = zcntr;

}
void k100_DetectorConstruction::SetConstructShieldTestEnvironmentSize(G4double sizel,G4double sizew,G4double sizethk)
{
  shieldTestParams.sizel = sizel;
  shieldTestParams.sizew = sizew;
  shieldTestParams.sizethk = sizethk;

}
G4String k100_DetectorConstruction::GetConstructShieldTestEnvironmentMat()
{
  return shieldTestParams.shieldmaterial->GetName();
}
void k100_DetectorConstruction::SetConstructShieldTestEnvironmentMat(G4String mat)
{
  if(mat=="Lead"){
    shieldTestParams.shieldmaterial = shieldPbMat;
  }
  else if(mat=="Poly"){
    shieldTestParams.shieldmaterial = polyMat;
  }
  else{  //default to poly material
    shieldTestParams.shieldmaterial = polyMat;
  }

}

void k100_DetectorConstruction::ConstructTower(G4VPhysicalVolume* physicalDetectorBox)
{

  //----------------------------------- 
  // Contruct the Tower Logical Volume
  //-----------------------------------

  // Position Vector
  //G4ThreeVector positionTower = G4ThreeVector(xtow[1-1],ytow[1-1],ztow[1-1]);
  G4ThreeVector positionTower = G4ThreeVector(tower_x,tower_y,tower_z);
  G4RotationMatrix r180Rotation;		// flip towers over
  r180Rotation.rotateY(180.*deg);
  G4Transform3D towerflip(r180Rotation, positionTower);
  //r180Rotation->rotateX(M_PI*rad);
  //r180Rotation->rotateZ(M_PI*rad);
  G4Polyhedra* solidTower1 = new G4Polyhedra("Tower1_S",0.*deg,360.*deg,6,Tower_nZcut,Tower_zPcut,Tower_rIcut,Tower_rOcut);

  if(NbOfTowers>=1){
    // Tower 1
    G4LogicalVolume* logicalTower1 = new G4LogicalVolume(solidTower1,   // The solid
							 defaultMat,   // the material
							 "Tower1_L",  // its name
							 0,0,0);
    
    G4VPhysicalVolume* physicalTower1 = new G4PVPlacement(towerflip,                        // flipping with rotation here kills it, even putting no rotation here kills it.  Must use a G4Transform3D, but it works perfectly.
							  /*positionTower,*/           // in towerflip
							  "Tower1_P",             // name
							  logicalTower1,         // the logical volume to use
							  physicalWorld,  // the mother volume
							  false,               // no boolean operation
							  0);                 // copy number
    
    // Visualization attributes
    G4VisAttributes* VisAttTower1 = new G4VisAttributes(G4Colour(215/255.,215/255.,215/255.));
    VisAttTower1->SetForceSolid(false);
    logicalTower1->SetVisAttributes(VisAttTower1);  
    logicalTower1->SetVisAttributes(G4VisAttributes::Invisible);  // Make Invisible
    
    ConstructTowerGuts(physicalTower1);
    if(ConstructZipBool) {FillTheTower(physicalTower1,1);}
  } // ends tower 1

  if(NbOfTowers>=2){
    // Tower 2
    positionTower = G4ThreeVector(xtow[2-1],ytow[2-1],ztow[2-1]);
    G4LogicalVolume* logicalTower2 = new G4LogicalVolume(solidTower1,   // The solid
							 defaultMat,   // the material
							 "Tower2_L",  // its name
							 0,0,0);
    G4VPhysicalVolume* physicalTower2 = new G4PVPlacement(0,                        // no rotation
							  positionTower,           // at (x,y,z)
							  "Tower2_P",             // name
							  logicalTower2,         // the logical volume to use
							  physicalWorld,  // the mother volume
							  false,               // no boolean operation
							  0);                 // copy number
    // Visualization attributes
    G4VisAttributes* VisAttTower2 = new G4VisAttributes(G4Colour(215/255.,215/255.,215/255.));
    VisAttTower2->SetForceSolid(false);
    logicalTower2->SetVisAttributes(VisAttTower2);  
    logicalTower2->SetVisAttributes(G4VisAttributes::Invisible);  // Make Invisible
    
    ConstructTowerGuts(physicalTower2);
    if(ConstructZipBool) {FillTheTower(physicalTower2,2);}
  } // ends if 2 towers

  if(NbOfTowers>=3){
    // Tower 3
    positionTower = G4ThreeVector(xtow[3-1],ytow[3-1],ztow[3-1]);
    G4LogicalVolume* logicalTower3 = new G4LogicalVolume(solidTower1,   // The solid
							 defaultMat,   // the material
							 "Tower3_L",  // its name
							 0,0,0);
    G4VPhysicalVolume* physicalTower3 = new G4PVPlacement(0,                        // no rotation
							  positionTower,           // at (x,y,z)
							  "Tower3_P",             // name
							  logicalTower3,         // the logical volume to use
							  physicalWorld,        // the mother volume
							  false,               // no boolean operation
							  0);                 // copy number
    // Visualization attributes
    G4VisAttributes* VisAttTower3 = new G4VisAttributes(G4Colour(215/255.,215/255.,215/255.));
    VisAttTower3->SetForceSolid(false);
    logicalTower3->SetVisAttributes(VisAttTower3);  
    logicalTower3->SetVisAttributes(G4VisAttributes::Invisible);  // Make Invisible
    
    ConstructTowerGuts(physicalTower3);
    if(ConstructZipBool) {FillTheTower(physicalTower3,3);}
  }  // ends if 3 towers

  if(NbOfTowers>=4){
    // Tower 4
    positionTower = G4ThreeVector(xtow[4-1],ytow[4-1],ztow[4-1]);
    G4LogicalVolume* logicalTower4 = new G4LogicalVolume(solidTower1,   // The solid
							 defaultMat,   // the material
							 "Tower4_L",  // its name
							 0,0,0);
    G4VPhysicalVolume* physicalTower4 = new G4PVPlacement(0,                        // no rotation
							  positionTower,           // at (x,y,z)
							  "Tower4_P",             // name
							  logicalTower4,         // the logical volume to use
							  physicalWorld,        // the mother volume
							  false,               // no boolean operation
							  0);                 // copy number
    // Visualization attributes
    G4VisAttributes* VisAttTower4 = new G4VisAttributes(G4Colour(215/255.,215/255.,215/255.));
    VisAttTower4->SetForceSolid(false);
    logicalTower4->SetVisAttributes(VisAttTower4);  
    logicalTower4->SetVisAttributes(G4VisAttributes::Invisible);  // Make Invisible
    
    ConstructTowerGuts(physicalTower4);
    if(ConstructZipBool) {FillTheTower(physicalTower4,4);}
  } // ends if 4 towers

  if(NbOfTowers>=5){
    // Tower 5
    positionTower = G4ThreeVector(xtow[5-1],ytow[5-1],ztow[5-1]);
    G4LogicalVolume* logicalTower5 = new G4LogicalVolume(solidTower1,   // The solid
							 defaultMat,   // the material
							 "Tower5_L",  // its name
							 0,0,0);
    G4VPhysicalVolume* physicalTower5 = new G4PVPlacement(0,                        // no rotation
							  positionTower,           // at (x,y,z)
							  "Tower5_P",             // name
							  logicalTower5,         // the logical volume to use
							  physicalWorld,        // the mother volume
							  false,               // no boolean operation
							  0);                 // copy number
    // Visualization attributes
    G4VisAttributes* VisAttTower5 = new G4VisAttributes(G4Colour(215/255.,215/255.,215/255.));
    VisAttTower5->SetForceSolid(false);
    logicalTower5->SetVisAttributes(VisAttTower5);  
    logicalTower5->SetVisAttributes(G4VisAttributes::Invisible);  // Make Invisible
    
    ConstructTowerGuts(physicalTower5);
    if(ConstructZipBool) {FillTheTower(physicalTower5,5);}
  } // ends if 5 towers

 //  // Sensor Region (cuts 250 eV) for any particle within the DetBox logical volume
//   // G4Region* DetectorRegion;
//   // if(DetectorRegion) {G4cout << "\n###### DetectorRegion exits." << DetectorRegion << G4endl;}
  
//   DetectorRegion = new G4Region(G4String("Detector"));
//   logicalWorld->SetRegion(DetectorRegion);
//   DetectorRegion->AddRootLogicalVolume(logicalDetectorBox);
 
} // ends Tower construction

// ------------------------------------------------
// ------------------------------------------------

void k100_DetectorConstruction::FillTheTower(G4VPhysicalVolume* physicalTower, G4int towerNb)
{
  G4ThreeVector positionZipArray = G4ThreeVector(0,0,Tower_zPcut[0]-(zPldh[1]-zPldh[0])-zPsdh[1]-Zip_z/2);
 //Corrected from: G4ThreeVector(0,0,Tower_zPcut[0]+(zPldh[1]-zPldh[0]) + zPsdh[1]);
  G4Tubs* solidZipArray = new G4Tubs("ZipArray_S", 0.0, Zip_Rout, 6*Zip_Househeight/2,  0, 2*pi);

  //------------------------------ 
  // Individual Zips
  //------------------------------
  // 
  // An example of Parameterised volumes
  // dummy values for G4Tubs -- modified by parameterised volume

  // simple  ZIP geometry: Plain cylinder
  // G4Tubs* solidZip = new G4Tubs("Zip_S", 0.0, Zip_Rout*0.5, Zip_Househeight*0.5,  0, 2*pi);
  
  //including the real ZIP geometry, by building a Boolean Solid added by T.Bruch 02/07
  G4Tubs* solidZipCircle = new G4Tubs("Zip_Scircle", 0.0, Zip_Rout, Zip_z*0.5,  0, 2*pi);
  G4Box* solidZipBox=new G4Box("Zip_SBox",Zip_Flat_R1,Zip_Flat_R2,Zip_z*0.5);
  G4IntersectionSolid* solidZip=new G4IntersectionSolid("Zip_S",solidZipCircle,solidZipBox);
  
  // Visualization attributes
  G4VisAttributes* VisAttZip = new G4VisAttributes(G4Colour(128/255.,0/255.,0/255.));
  VisAttZip->SetForceSolid(DrawSolidZipBool); 

  // Using the parameterisation defined in the PixelParameteriation class.
  zipParam = 
    new k100_ZipParameterisation(NbZipsPerTower,   // NoZips/Tower
				 Zip_Househeight, // Z spacing of centers
				 Zip_Rout,       // Zip Radius
				 Zip_z,         // Depth of pixel (in Z)
				 DetMaterials, // Zip Materials
				 zipGeMat, zipSiMat,
				 DrawSolidZipBool,
				 towerNb); 
  
  // Prepare to declare sensitive detectors
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  
  if(towerNb==1) {
    G4LogicalVolume* logicalZip1Array = new G4LogicalVolume(solidZipArray, defaultMat, "ZipArray1_L", 0,0,0);
    G4VPhysicalVolume* physicalZip1Array = new G4PVPlacement(0, positionZipArray, "ZipArray1_P", logicalZip1Array,
							     physicalTower, false,  0);
    // Visualization attributes
    G4VisAttributes* VisAttZipArr = new G4VisAttributes(G4Colour(204/255.,255/255.,255/255.));
    VisAttZipArr->SetForceSolid(true);
    logicalZip1Array->SetVisAttributes(VisAttZipArr);  
    logicalZip1Array->SetVisAttributes(G4VisAttributes::Invisible);  // Make Invisible

    G4LogicalVolume* logicalZip1 = new G4LogicalVolume(solidZip, defaultMat, "Zip1_L", 0,0,0);
    logicalZip1->SetVisAttributes(VisAttZip);  

    //G4VPhysicalVolume* physicalZip =  // I gues this part is not needed
    new G4PVParameterised("Zip1_P",          // name
			  logicalZip1,       // logical volume
			  logicalZip1Array,  // mother logical volume
			  kZAxis,            // place along the z axis
			  NbZipsPerTower,    // Number of chambers

			  zipParam,          // Use this parameterisation
			  false);

// Since the ZIPs coordinates are not very clear from above function calls here is some information to
// understand the results:
// 
// zip_center_position = 
//    Tower_coordinate + ZIP_array_position_in_tower + ZIP_position_inZIP_array
//
// all these appear in k100_DetectorConstruction in the following functions:
// physicalTower (positionTower) + physicalZip1Array (positionZipArray) +
//    G4PVParameterised(...,zipParam)
// ==>
// positionTower = xtow ; ztow ; ztow
// positionZipArray = 0 ; 0 ; Tower_zPcut[0]+(zPldh[1]-zPldh[0])+zPsdh[1]
//
// in k100_ZipParametrization.cc 
// zipParam{coordinates} = 0 ; 0 ; (2.5 - copyNo)*fZSpacing
//    where fZSpacing = Zip_Househeight from the function call in k100_DetectorConstruction
//
// in k100_DetectorParameterDef.icc
// xtow[0]=8.32368*cm ; ytow[0]=-4.80568*cm; ztow[0]=Tower_Center_Z=3.2797*cm
// Tower_zPcut[0]+(zPldh[1]-zPldh[0])+zPsdh[1]
//   = (-5.3108-0.2425-4.48691-0.2405/2-0.0685-3.6385-0.07584) + 0.07584*2 + 3.6385 
//   = -13.9433 + 0.1517 + 3.6385 = -10.1531*cm
// Zip_Househeight = 1.20142*cm
//
// ==> coordinated for the center of ZIP1:
// x_ZIP1 = 8.32368*cm 
// y_ZIP1 = -4.80568*cm 
// z_ZIP1 = 3.2797 - 10.1531 + 2.5*1.20142 = -3.8699*cm

	
	
    //------------------------------------------------ 
    // Sensitive detectors
    //------------------------------------------------ 
    
    G4String detectorZipSDname = "tower1";
    G4int collID = -1; collID = SDman->GetCollectionID(detectorZipSDname);
    k100_ZipSD* azipSD1;
    ConstructGenericSensitiveInt=1; 
    azipSD1 = new k100_ZipSD(detectorZipSDname, towerNb);
    k100CollName[detectorZipSDname] = towerNb;
    SDman->AddNewDetector(azipSD1);
    //    G4cout << "#### DetCon : zipCollID[ii]  " << SDman->GetCollectionID(detectorZipSDname) << G4endl;
    logicalZip1->SetSensitiveDetector(azipSD1);
  }


  if(towerNb==2) {
    G4LogicalVolume* logicalZip2Array = new G4LogicalVolume(solidZipArray, defaultMat, "ZipArray2_L", 0,0,0);
    G4VPhysicalVolume* physicalZip2Array = new G4PVPlacement(0, positionZipArray, "ZipArray2_P", logicalZip2Array,
							     physicalTower, false,  0);
    // Visualization attributes
    G4VisAttributes* VisAttZipArr = new G4VisAttributes(G4Colour(204/255.,255/255.,255/255.));
    VisAttZipArr->SetForceSolid(true);
    logicalZip2Array->SetVisAttributes(VisAttZipArr);  
    logicalZip2Array->SetVisAttributes(G4VisAttributes::Invisible);  // Make Invisible

    G4LogicalVolume* logicalZip2 = new G4LogicalVolume(solidZip, defaultMat, "Zip2_L", 0,0,0);
    logicalZip2->SetVisAttributes(VisAttZip);  

    //G4VPhysicalVolume* physicalZip =  // I guess this part is not needed
    new G4PVParameterised("Zip2_P", logicalZip2, physicalZip2Array, kZAxis, NbZipsPerTower, zipParam);
    
    //------------------------------------------------ 
    // Sensitive detectors
    //------------------------------------------------ 
        
    G4String detectorZipSDname = "tower2";
    G4int collID = -1; collID = SDman->GetCollectionID(detectorZipSDname);
    k100_ZipSD* azipSD2;
    if(collID<0){ 
      azipSD2 = new k100_ZipSD(detectorZipSDname, towerNb);
      SDman->AddNewDetector(azipSD2);
    }
    //    G4cout << "#### DetCon : zipCollID[ii]  " << SDman->GetCollectionID(detectorZipSDname) << G4endl;
    logicalZip2->SetSensitiveDetector(azipSD2);
  }


  if(towerNb==3) {
    G4LogicalVolume* logicalZip3Array = new G4LogicalVolume(solidZipArray, defaultMat, "ZipArray3_L", 0,0,0);
    G4VPhysicalVolume* physicalZip3Array = new G4PVPlacement(0, positionZipArray, "ZipArray3_P", logicalZip3Array,
							     physicalTower, false,  0);
    // Visualization attributes
    G4VisAttributes* VisAttZipArr = new G4VisAttributes(G4Colour(204/255.,255/255.,255/255.));
    VisAttZipArr->SetForceSolid(true);
    logicalZip3Array->SetVisAttributes(VisAttZipArr);  
    logicalZip3Array->SetVisAttributes(G4VisAttributes::Invisible);  // Make Invisible

    G4LogicalVolume* logicalZip3 = new G4LogicalVolume(solidZip, defaultMat, "Zip3_L", 0,0,0);
    logicalZip3->SetVisAttributes(VisAttZip);  

    //G4VPhysicalVolume* physicalZip =  // I guess this part is not needed
    new G4PVParameterised("Zip3_P", logicalZip3, physicalZip3Array, kZAxis, NbZipsPerTower, zipParam);
    
    //------------------------------------------------ 
    // Sensitive detectors
    //------------------------------------------------ 
        
    G4String detectorZipSDname = "tower3";
    G4int collID = -1; collID = SDman->GetCollectionID(detectorZipSDname);
    k100_ZipSD* azipSD3;
    if(collID<0){ 
      azipSD3 = new k100_ZipSD(detectorZipSDname, towerNb);
      SDman->AddNewDetector(azipSD3);
    }
    //    G4cout << "#### DetCon : zipCollID[ii]  " << SDman->GetCollectionID(detectorZipSDname) << G4endl;
    logicalZip3->SetSensitiveDetector(azipSD3);
  }

  if(towerNb==4) {
    G4LogicalVolume* logicalZip4Array = new G4LogicalVolume(solidZipArray, defaultMat, "ZipArray4_L", 0,0,0);
    G4VPhysicalVolume* physicalZip4Array = new G4PVPlacement(0, positionZipArray, "ZipArray4_P", logicalZip4Array,
							     physicalTower, false,  0);
    // Visualization attributes
    G4VisAttributes* VisAttZipArr = new G4VisAttributes(G4Colour(204/255.,255/255.,255/255.));
    VisAttZipArr->SetForceSolid(true);
    logicalZip4Array->SetVisAttributes(VisAttZipArr);  
    logicalZip4Array->SetVisAttributes(G4VisAttributes::Invisible);  // Make Invisible

    G4LogicalVolume* logicalZip4 = new G4LogicalVolume(solidZip, defaultMat, "Zip4_L", 0,0,0);
    logicalZip4->SetVisAttributes(VisAttZip);  

    //G4VPhysicalVolume* physicalZip =  // I guess this part is not needed
    new G4PVParameterised("Zip4_P", logicalZip4, physicalZip4Array, kZAxis, NbZipsPerTower, zipParam);
    
    //------------------------------------------------ 
    // Sensitive detectors
    //------------------------------------------------ 
        
    G4String detectorZipSDname = "tower4";
    G4int collID = -1; collID = SDman->GetCollectionID(detectorZipSDname);
    k100_ZipSD* azipSD4;
    if(collID<0){ 
      azipSD4 = new k100_ZipSD(detectorZipSDname, towerNb);
      SDman->AddNewDetector(azipSD4);
    }
    //    G4cout << "#### DetCon : zipCollID[ii]  " << SDman->GetCollectionID(detectorZipSDname) << G4endl;
    logicalZip4->SetSensitiveDetector(azipSD4);
  }

  if(towerNb==5) {
    G4LogicalVolume* logicalZip5Array = new G4LogicalVolume(solidZipArray, defaultMat, "ZipArray5_L", 0,0,0);
    G4VPhysicalVolume* physicalZip5Array = new G4PVPlacement(0, positionZipArray, "ZipArray5_P", logicalZip5Array,
							     physicalTower, false,  0);
    // Visualization attributes
    G4VisAttributes* VisAttZipArr = new G4VisAttributes(G4Colour(204/255.,255/255.,255/255.));
    VisAttZipArr->SetForceSolid(true);
    logicalZip5Array->SetVisAttributes(VisAttZipArr);  
    logicalZip5Array->SetVisAttributes(G4VisAttributes::Invisible);  // Make Invisible

    G4LogicalVolume* logicalZip5 = new G4LogicalVolume(solidZip, defaultMat, "Zip5_L", 0,0,0);
    logicalZip5->SetVisAttributes(VisAttZip);  

    //G4VPhysicalVolume* physicalZip =  // I guess this part is not needed
    new G4PVParameterised("Zip5_P", logicalZip5, physicalZip5Array, kZAxis, NbZipsPerTower, zipParam);
    
    //------------------------------------------------ 
    // Sensitive detectors
    //------------------------------------------------ 
        
    G4String detectorZipSDname = "tower5";
    G4int collID = -1; collID = SDman->GetCollectionID(detectorZipSDname);
    k100_ZipSD* azipSD5;
    if(collID<0){ 
      azipSD5 = new k100_ZipSD(detectorZipSDname, towerNb);
      SDman->AddNewDetector(azipSD5);
    }
    //    G4cout << "#### DetCon : zipCollID[ii]  " << SDman->GetCollectionID(detectorZipSDname) << G4endl;
    logicalZip5->SetSensitiveDetector(azipSD5);
  }

} // ends FillTheTower

// ------------------------------------------------
// ------------------------------------------------

void k100_DetectorConstruction::ConstructTowerGuts(G4VPhysicalVolume* physicalTower)
{

  // Visualization attributes
  G4VisAttributes* VisAttCu1 = new G4VisAttributes(G4Colour(180/255.,140/255.,60/255.)); VisAttCu1->SetForceSolid(DrawSolidTowerBool); 
  G4VisAttributes* VisAttCu2 = new G4VisAttributes(G4Colour(220/255.,165/255.,55/255.)); VisAttCu2->SetForceSolid(DrawSolidTowerBool); 
  G4VisAttributes* VisAttCu3 = new G4VisAttributes(G4Colour(210/255.,165/255.,70/255.)); VisAttCu3->SetForceSolid(DrawSolidTowerBool); 
  G4VisAttributes* VisAttCu4 = new G4VisAttributes(G4Colour(128/255.,64/255.,0/255.)); VisAttCu4->SetForceSolid(DrawSolidTowerBool); 

  //-------------------------------------------------------------------------
  //Working from the top down: this is the `copper upper tower', mass 2.0 kg
  //-------------------------------------------------------------------------
  G4ThreeVector position_top = G4ThreeVector(0,0,+Tower_zPcut[1] - zPcut[1]);
  G4Polyhedra* cu_tops=new G4Polyhedra("cu_top",0.*deg,360.*deg,6,nZcut,zPcut,rIcut,rOcut);
  G4LogicalVolume* cu_topl1 = new G4LogicalVolume(cu_tops,towerMat,"cutl1");
  G4PVPlacement* cu_topp1 = new G4PVPlacement(0,position_top,"cutp1",cu_topl1,physicalTower,false,0);
  cu_topl1->SetVisAttributes(VisAttCu1);  
  
  //----------------------------------------------------------------------
  //Next, the lower cap of the `copper upper tower', mass 0.233 kg/tower
  //----------------------------------------------------------------------
  G4ThreeVector position_clc = G4ThreeVector(0,0,+Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - zPclc[1]);
  G4Polyhedra* cu_clcs=new G4Polyhedra("cu_clc",0.*deg,360.*deg,6,nZclc,zPclc,rIclc,rOclc);
  G4LogicalVolume* cu_clcl1 = new G4LogicalVolume(cu_clcs,towerMat,"clcl1");
  G4PVPlacement* cu_clcp1 = new G4PVPlacement(0,position_clc,"clcp1",cu_clcl1,physicalTower,false,0);
  cu_clcl1->SetVisAttributes(VisAttCu2); 
     
  //------------------------------------------------------------------
  //Next, connector tube, mass 0.140 kg/tower (next ring excluded...)
  //------------------------------------------------------------------
  G4ThreeVector position_ctu = G4ThreeVector(0,0,+Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - zHctu);
  G4Tubs* cu_ctus=new G4Tubs("cu_ctu",rIctu,rOctu,zHctu,0.*deg,360.*deg);
  G4LogicalVolume* cu_ctul1=new G4LogicalVolume(cu_ctus,towerMat,"ctul1");
  G4PVPlacement* cu_ctup1=new G4PVPlacement(0,position_ctu,"ctup1",cu_ctul1,physicalTower,false,0);
  cu_ctul1->SetVisAttributes(VisAttCu3); 


  //--------------------------------------------------------------
  //Next, ring at base of connector tube, mass 0.051 kg/tower
  //--------------------------------------------------------------
  G4ThreeVector position_rb = G4ThreeVector(0,0, Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu); // QUESTION : Is this ring centered on the bottom of the tower ??
  G4Tubs* cu_rbs=new G4Tubs("cu_rb",rIrb,rOrb,zHrb,0.*deg,360.*deg);
  G4LogicalVolume* cu_rbl1=new G4LogicalVolume(cu_rbs,towerMat,"crbl1");
  G4PVPlacement* cu_rbp1=new G4PVPlacement(0,position_rb,"crbp1",cu_rbl1,physicalTower,false,0);
  cu_rbl1->SetVisAttributes(VisAttCu1); 

  //------------------------------------------------------------------
  //Next, the upper cap of the detector housing, mass 0.071 kg/tower
  // (about 0.044 kg/tower actually included in ring at base above)
  //------------------------------------------------------------------
  G4ThreeVector position_udh = G4ThreeVector(0,0, Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu - zHrb - zPudh[1]); 
  G4Polyhedra* cu_udhs=new G4Polyhedra("cu_udh",0.*deg,360.*deg,6,nZudh,zPudh,rIudh,rOudh);
  G4LogicalVolume* cu_udhl1 = new G4LogicalVolume(cu_udhs,towerMat,"udhl1");
  G4PVPlacement* cu_udhp1 = new G4PVPlacement(0, position_udh,"udhp1",cu_udhl1,physicalTower,false,0);
  cu_udhl1->SetVisAttributes(VisAttCu2);
  

  //--------------------------------------------------------------
  //Next, the side detector housing, mass 0.384 kg/tower
  //--------------------------------------------------------------
  G4ThreeVector position_sdh = G4ThreeVector(0,0, Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu - zHrb - (zPudh[1]-zPudh[0]) - zPsdh[1]); 
  G4Polyhedra* cu_sdhs=new G4Polyhedra("cu_sdh",0.*deg,360.*deg,6,nZsdh,zPsdh,rIsdh,rOsdh);
  G4LogicalVolume* cu_sdhl1 = new G4LogicalVolume(cu_sdhs,towerMat,"sdhl1");
  G4PVPlacement* cu_sdhp1 = new G4PVPlacement(0,position_sdh,"sdhp1",cu_sdhl1,physicalTower,false,0);
  cu_sdhl1->SetVisAttributes(VisAttCu3);
  //cu_sdhl1->SetVisAttributes(G4VisAttributes::Invisible);  // Make Invisible


  //------------------------------------------------------------------
  //Next, the lower cap of the detector housing, mass 0.078 kg/tower
  //------------------------------------------------------------------
  G4ThreeVector position_ldh = G4ThreeVector(0,0, Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu - zHrb - (zPudh[1]-zPudh[0]) - (zPsdh[1]-zPsdh[0]) - zPldh[1] ); 
  G4Polyhedra* cu_ldhs=new G4Polyhedra("cu_ldh",0.*deg,360.*deg,6,nZldh,zPldh,rIldh,rOldh);
  G4LogicalVolume* cu_ldhl1 = new G4LogicalVolume(cu_ldhs,towerMat,"ldhl1");
  G4PVPlacement* cu_ldhp1 = new G4PVPlacement(0,position_ldh,"ldhp1",cu_ldhl1,physicalTower,false,0);
  cu_ldhl1->SetVisAttributes(VisAttCu2);



  //------------------------------------------------------------------
  //Next, the Side Coaxes
  //------------------------------------------------------------------

  //rotation matrices for side coaxes                                                     
  G4ThreeVector newz(0.0,0.0,1.0);  //z-axis doesn't change for any coax                  
  
  G4ThreeVector newx1(-cos(60*deg),sin(60*deg),0.0);
  G4ThreeVector newy1(-sin(60*deg),-cos(60*deg),0.0);
  G4RotationMatrix* coaxrotation1=new G4RotationMatrix(newx1,newy1,newz);
  
  G4ThreeVector newx2(-cos(60*deg),-sin(60*deg),0.0);
  G4ThreeVector newy2(sin(60*deg),-cos(60*deg),0.0);
  G4RotationMatrix* coaxrotation2=new G4RotationMatrix(newx2,newy2,newz);

  G4Box* solidcoax1=new G4Box("scoax1", sdcx_widthH, sdcx_thicknessH, sdcx_lenH[1-1]);
  G4Box* solidcoax2=new G4Box("scoax2", sdcx_widthH, sdcx_thicknessH, sdcx_lenH[2-1]);
  G4Box* solidcoax3=new G4Box("scoax3", sdcx_widthH, sdcx_thicknessH, sdcx_lenH[3-1]);
  G4Box* solidcoax4=new G4Box("scoax4", sdcx_widthH, sdcx_thicknessH, sdcx_lenH[4-1]);
  G4Box* solidcoax5=new G4Box("scoax5", sdcx_widthH, sdcx_thicknessH, sdcx_lenH[5-1]);
  G4Box* solidcoax6=new G4Box("scoax6", sdcx_widthH, sdcx_thicknessH, sdcx_lenH[6-1]);
  G4LogicalVolume* logiccoax1=new G4LogicalVolume(solidcoax1,towerMat,"LogCoax1");
  G4LogicalVolume* logiccoax2=new G4LogicalVolume(solidcoax2,towerMat,"LogCoax2");
  G4LogicalVolume* logiccoax3=new G4LogicalVolume(solidcoax3,towerMat,"LogCoax3");
  G4LogicalVolume* logiccoax4=new G4LogicalVolume(solidcoax4,towerMat,"LogCoax4");
  G4LogicalVolume* logiccoax5=new G4LogicalVolume(solidcoax5,towerMat,"LogCoax5");
  G4LogicalVolume* logiccoax6=new G4LogicalVolume(solidcoax6,towerMat,"LogCoax6");
  G4LogicalVolume* logiccoax7=new G4LogicalVolume(solidcoax1,towerMat,"LogCoax7");
  G4LogicalVolume* logiccoax8=new G4LogicalVolume(solidcoax2,towerMat,"LogCoax8");
  G4LogicalVolume* logiccoax9=new G4LogicalVolume(solidcoax3,towerMat,"LogCoax9");
  G4LogicalVolume* logiccoax10=new G4LogicalVolume(solidcoax4,towerMat,"LogCoax10");
  G4LogicalVolume* logiccoax11=new G4LogicalVolume(solidcoax5,towerMat,"LogCoax11");
  G4LogicalVolume* logiccoax12=new G4LogicalVolume(solidcoax6,towerMat,"LogCoax12");
  
  // Visualization line
  logiccoax1->SetVisAttributes(VisAttCu4);
  logiccoax2->SetVisAttributes(VisAttCu4);
  logiccoax3->SetVisAttributes(VisAttCu4);
  logiccoax4->SetVisAttributes(VisAttCu4);
  logiccoax5->SetVisAttributes(VisAttCu4);
  logiccoax6->SetVisAttributes(VisAttCu4);
    
  //work down the tower, placing housing, detectors, coaxes as we go

  //   G4double ztower=-4.2927*cm;
  //   G4double househeight=1.20142*cm;
  //   G4double coaxpos=ztower+4.86*cm;

  G4double coaxpos;
  coaxpos = (Tower_zPcut[0]+zPldh[1] + (6-0.5-5)*Zip_Househeight + sdcx_lenH[0]);
  G4PVPlacement* physcoax1 = new G4PVPlacement(0,G4ThreeVector(0.*cm,-(sdcx_thicknessH+rOsdh[0]),coaxpos),"coax1",logiccoax1,physicalTower,false,0);

  //G4PVPlacement* physcoax1 = new G4PVPlacement(0,G4ThreeVector(0.*cm,-4.18*cm,coaxpos),"coax1",logiccoax1,physicalTower,false,0);

  coaxpos = (Tower_zPcut[0]+zPldh[1] + (6-0.5-1)*Zip_Househeight + sdcx_lenH[1]);
 // G4PVPlacement* physcoax2 = new G4PVPlacement(coaxrotation1,G4ThreeVector(0+3.62*cm,0-2.09*cm,coaxpos),"coax2",logiccoax2,physicalTower,false,0);

  coaxpos = (Tower_zPcut[0]+zPldh[1] + (6-0.5-2)*Zip_Househeight + sdcx_lenH[2]);
  //G4PVPlacement* physcoax3 = new G4PVPlacement(coaxrotation2,G4ThreeVector(0+3.62*cm,0+2.09*cm,coaxpos),"coax3",logiccoax3,physicalTower,false,0);

  coaxpos = (Tower_zPcut[0]+zPldh[1] + (6-0.5-5)*Zip_Househeight + sdcx_lenH[3]);
  G4PVPlacement* physcoax4 = new G4PVPlacement(0,G4ThreeVector(0,sdcx_thicknessH+rOsdh[0],coaxpos),"coax4",logiccoax4,physicalTower,false,0);  

  //G4PVPlacement* physcoax4 = new G4PVPlacement(0,G4ThreeVector(0,0+4.18*cm,coaxpos),"coax4",logiccoax4,physicalTower,false,0);

  coaxpos = (Tower_zPcut[0]+zPldh[1] + (6-0.5-4)*Zip_Househeight + sdcx_lenH[4]);
  //G4PVPlacement* physcoax5 = new G4PVPlacement(coaxrotation1,G4ThreeVector(0-3.62*cm,0+2.09*cm,coaxpos),"coax5",logiccoax5,physicalTower,false,0);

  coaxpos = (Tower_zPcut[0]+zPldh[1] + (6-0.5-5)*Zip_Househeight + sdcx_lenH[5]);
  G4PVPlacement* physcoax6 = new G4PVPlacement(coaxrotation2,G4ThreeVector(-(sdcx_thicknessH+rOsdh[0])*sin(60*deg),-(sdcx_thicknessH+rOsdh[0])*cos(60*deg),coaxpos),"coax6",logiccoax6,physicalTower,false,0);

  //G4PVPlacement* physcoax6 = new G4PVPlacement(coaxrotation2,G4ThreeVector(0-3.62*cm,0-2.09*cm,coaxpos),"coax6",logiccoax6,physicalTower,false,0);

} // ends Tower Guts construction

// ------------------------------------------------
// ------------------------------------------------

void k100_DetectorConstruction::ConstructEverything(G4LogicalVolume*  logicalWorld)
{

  //------------------------------ 
  // Here we construct the enitre
  // experimental apparatus
  //------------------------------ 

  if(ConstructVetoBool)    {ConstructVeto(logicalWorld);}
  if(ConstructShieldsBool) {ConstructShields(logicalWorld);}
  if(ConstructIceBoxBool)  {ConstructIceBox(logicalWorld);}
  if(ConstructThermalNeutronBoxBool)  {ConstructThermalNeutronBox(physicalWorld);}

} // ends ConstructEverything

// ------------------------------------------------
// ------------------------------------------------

void k100_DetectorConstruction::ConstructVeto(G4LogicalVolume*  logicalWorld)
{
} // ends Veto panel construction

// ------------------------------------------------
// ------------------------------------------------


void k100_DetectorConstruction::ConstructShields(G4LogicalVolume*  logicalWorld)
{
  // BJORN WILL USE THIS AREA FOR NEUTRON SHIELD PANELS
  
  // ------------------- Aluminum frames ------------------
  // setup and initial solids
  G4ThreeVector framePosition = G4ThreeVector(frame_x,frame_y,frame_z); // TODO:  Set proper position in include/k100vars.hh
  G4Box* vertBar = new G4Box("vertBar",2.54*cm,2.54*cm,27.5*2.54*cm); //Corrected:Originally 1x1x27 in
  G4Box* horizBar = new G4Box("horizBar",9*2.54*cm,2.54*cm,2.54*cm); //Corrected: Originally 11x1x1 in
  G4Box* base = new G4Box("base",27.0*2.54*cm,27*2.54*cm,2.0*2.54*cm); //New:Base of frame 
  G4Box*  baseHole1 = new G4Box("baseHole1",9*2.45*cm,11.5*2.54*cm,2.5*2.54*cm); //New: Holes in Base
  G4Box*  baseHole2 = new G4Box("baseHole2",7*2.45*cm,11.5*2.54*cm,2.5*2.54*cm); //New: Holes in Base
  

  G4RotationMatrix noRot;
  G4ThreeVector barPos(20*2.54*cm,0,0);
  G4Transform3D barTrans(noRot, barPos);
  G4UnionSolid* oldFrame = new G4UnionSolid("oldFrame",vertBar,vertBar,barTrans);
  G4UnionSolid* frame;

  // Loop to add other two vertical bars
  for(int i=0; i<2; i++)
    {
      barPos = G4ThreeVector(20*2.54*i*cm,24.0*2.54*cm,0); //Corrected: Originally: 20ix24x0 in 
      barTrans=G4Transform3D(noRot, barPos); // do we need this line? yes
      frame=new G4UnionSolid("frame", oldFrame, vertBar, barTrans);
      oldFrame=frame;
    }
  
  // Loop to add four horizontal bars
  for(int i=0; i<2; i++)
    {
      for(int j=0; j<2; j++)
	{
	  barPos = G4ThreeVector(10*2.54*cm,i*24*2.54*cm,(j*33.0-12.5)*2.54*cm); //Corrected: Originally 10x24ix40j-20 in
	  barTrans=G4Transform3D(noRot, barPos); // do we need this line? yes
	  frame=new G4UnionSolid("frame", oldFrame, horizBar, barTrans);
	  oldFrame=frame;
	}
    }

  //New: Add holes to base
  G4ThreeVector baseHolePos(0*2.54*cm,13.5*2.54*cm,0*2.54*cm);
  G4Transform3D baseOffset(noRot,baseHolePos);
  G4SubtractionSolid* oldBase= new G4SubtractionSolid("oldBase",base,baseHole1,baseOffset);
  G4SubtractionSolid* newBase;

   baseHolePos = G4ThreeVector(0*2.54*cm,-13.5*2.54*cm,0*cm);
   baseOffset = G4Transform3D(noRot,baseHolePos);
   newBase = new G4SubtractionSolid("newBase",oldBase,baseHole1,baseOffset);
   oldBase=newBase;

  for(int i=0;i<2;i++)
    {
      for(int j=0;j<2;j++)
        {
          baseHolePos = G4ThreeVector((2*i-1)*18*2.54*cm,(2*j-1)*13.5*2.54*cm,0*cm);
           baseOffset = G4Transform3D(noRot,baseHolePos);
           newBase = new G4SubtractionSolid("newBase",oldBase,baseHole2,baseOffset);
           oldBase=newBase;
        }
    }

  //New: Add base
  G4ThreeVector basePos(10*2.54*cm,12*2.54*cm,-29.5*2.54*cm);
  G4Transform3D baseTrans(noRot,basePos);
  frame = new G4UnionSolid("frame",oldFrame,oldBase,baseTrans);
  oldFrame=frame;

  //Add rods
  G4Tubs* rod = new G4Tubs("rod", 0, .5*2.54*cm, 5.1*2.54*cm, 0, 2*pi);
  G4RotationMatrix rodRot;
  rodRot.rotateY(pi/2*rad);
  
  // Loop to add rods (end rods)
  for(int i=0; i<2; i++)
    {
      for(int j=0; j<2; j++)
	{
	  for(int k=0; k<2; k++)
	    {
	      G4ThreeVector rodPos((i*29.8-4.9)*2.54*cm,j*24*2.54*cm,(k*25-8.5)*2.54*cm); //Corrected: Originally 30i-5x24jx25k-8.5 in
	      G4Transform3D rodTrans(rodRot,rodPos);
	      frame=new G4UnionSolid("frame", oldFrame, rod, rodTrans);
	      oldFrame=frame;
	    }
	}
    }

  rodRot.rotateZ(pi/2*rad);

  // Loop to add rods (side rods)
  for(int i=0; i<2; i++)
    {
      for(int j=0; j<2; j++)
	{
	  for(int k=0; k<2; k++)
	    {
	      G4ThreeVector rodPos(i*20*2.54*cm,(j*33.8-4.9)*2.54*cm,(k*25-6.5)*2.54*cm); //Corrected: Originally 20ix34j-5x25k-6.5 in
	      G4Transform3D rodTrans(rodRot,rodPos);
	      frame=new G4UnionSolid("frame", oldFrame, rod, rodTrans);
	      oldFrame=frame;
	    }
	}
    }

  // place the aluminum frames
  G4LogicalVolume* logicFrame = new G4LogicalVolume(frame,aluminum,"logicFrame",0,0,0);
  new G4PVPlacement(0,framePosition,"physicFrame", logicFrame, physicalWorld, false, 0);
  
  // frame visuals
  G4VisAttributes* frameVis = new G4VisAttributes(G4Colour(128/255.,128/255.,128/255.));
  frameVis->SetForceSolid(true);
  logicFrame->SetVisAttributes(frameVis);

  // --------------------- Poly Panels --------------------------
  // poly visuals
  G4VisAttributes* polyVis = new G4VisAttributes(G4Colour(255/255.,255/255.,255/255.));
  polyVis->SetForceSolid(false);
  
  // base panel setup
  G4ThreeVector panelPosition = G4ThreeVector(frame_x+10*2.54*cm,frame_y+12*2.54*cm,frame_z-23.5*2.54*cm);
  G4Box* baseShield = new G4Box("baseShield", 11*2.54*cm, 11*2.54*cm, 4*2.54*cm);
  G4LogicalVolume* logicBase = new G4LogicalVolume(baseShield, polyMat, "logicBase",0,0,0);
  new G4PVPlacement(0,panelPosition,"physicBase",logicBase,physicalWorld,false,0);
  logicBase->SetVisAttributes(polyVis);

  // square side panels
  // initial setup
  G4Box* squareBase = new G4Box("squareBase",11*2.54*cm,4*2.54*cm,27.5*2.54*cm);
  G4Tubs* hole = new G4Tubs("hole", 0, .55*2.54*cm, 2.1*2.54*cm, 0, 2*pi);
  G4ThreeVector off(1*m,1*m,1*m);
  G4RotationMatrix holeRot;
  holeRot.rotateX(pi/2*rad);
  G4Transform3D offset(holeRot,off);
  G4SubtractionSolid* oldSquare = new G4SubtractionSolid("oldSquare", squareBase, hole, offset);
  G4SubtractionSolid* square;

  // Loop to put in holes
  for(int i=-1; i<2; i++)
    {
      for(int j=-1; j<2; j++)
	{
	  if((i!=0)&&(j!=0))
	    {
	      G4ThreeVector holePos(i*10*2.54*cm,0,(j*12.5+6)*2.54*cm);
	      G4Transform3D holeTrans(holeRot,holePos);
	      square = new G4SubtractionSolid("square", oldSquare, hole, holeTrans);
	      oldSquare = square;
	    }
	}
    }

  // place squares
  G4LogicalVolume* logicSquare = new G4LogicalVolume(square,polyMat,"logicSquare",0,0,0);
  panelPosition=G4ThreeVector(frame_x+10*2.54*cm,frame_y-5*2.54*cm,frame_z);
  new G4PVPlacement(0,panelPosition,"physicSquare",logicSquare,physicalWorld,false,0);
  panelPosition=G4ThreeVector(frame_x+10*2.54*cm,frame_y+29*2.54*cm,frame_z);
  new G4PVPlacement(0,panelPosition,"physicSquare1",logicSquare,physicalWorld,false,1);
  logicSquare->SetVisAttributes(polyVis);

  // rect side panels
  // initial setup
  G4Box* rectBase = new G4Box("rectBase",4*2.54*cm,21*2.54*cm,27.5*2.54*cm);
  holeRot.rotateZ(pi/2*rad);
  G4Transform3D offs(holeRot,off);
  G4SubtractionSolid* oldRect = new G4SubtractionSolid("oldRect", rectBase, hole, offs);
  G4SubtractionSolid* rect;

  // Loop to put in holes
  for(int i=-1; i<2; i++)
    {
      for(int j=-1; j<2; j++)
	{
	  if((i!=0)&&(j!=0))
	    {
	      G4ThreeVector holePos(0,i*12*2.54*cm,(j*12.5+4)*2.54*cm);
	      G4Transform3D holeTrans(holeRot,holePos);
	      rect = new G4SubtractionSolid("rect", oldRect, hole, holeTrans);
	      oldRect = rect;
	    }
	}
    }

  // place rects
  G4LogicalVolume* logicRect = new G4LogicalVolume(rect,polyMat,"logicRect",0,0,0);
  panelPosition=G4ThreeVector(frame_x-5*2.54*cm,frame_y+12*2.54*cm,frame_z);
  new G4PVPlacement(0,panelPosition,"physicRect",logicRect,physicalWorld,false,0);
  panelPosition=G4ThreeVector(frame_x+25*2.54*cm,frame_y+12*2.54*cm,frame_z);
  new G4PVPlacement(0,panelPosition,"physicRect1",logicRect,physicalWorld,false,1);
  logicRect->SetVisAttributes(polyVis);


 // --------------------- Lead Frame Panels --------------------------
  // This section contains the aluminum fram that surrounds the lead shield.

  // lead frame visuals
  G4VisAttributes* lframeVis = new G4VisAttributes(G4Colour(128/255.,128/255.,128/255.));
  lframeVis->SetForceSolid(false);
 
 // --------------------- Lead Panels --------------------------

  // lead visuals
  G4VisAttributes* leadVis = new G4VisAttributes(G4Colour(102/255.,0/255.,187/255.));
  leadVis->SetForceSolid(false);

  // base lead panel setup
  G4ThreeVector leadPosition = G4ThreeVector(frame_x+10*2.54*cm,frame_y+12*2.54*cm,frame_z-27.75*2.54*cm);
  G4Box* leadBaseShield = new G4Box("leadBaseShield", 21*2.54*cm, 23*2.54*cm, 0.25*2.54*cm);
  G4Box* leadBaseHole = new G4Box("leadBaseHole", 1.1*2.54*cm,1.1*2.54*cm,.3*2.54*cm);
  
  G4ThreeVector baseHole(-10*2.54*cm,-12*2.54*cm,0*2.54*cm);
  G4Transform3D frameOffSet(noRot,baseHole);
  G4SubtractionSolid* oldLeadBase = new G4SubtractionSolid("oldLeadBase",leadBaseShield,leadBaseHole,frameOffSet);
  G4SubtractionSolid* leadBase;
  for(int i=0;i<2;i++)
    {
      for(int j=0;j<2;j++)
        {
        G4ThreeVector holePos((2*i-1)*10*2.54*cm,(2*j-1)*12*2.54*cm,0*2.54*cm);
        G4Transform3D frameOffSet(noRot,holePos);
        leadBase = new G4SubtractionSolid("leadBase",oldLeadBase,leadBaseHole,frameOffSet);
        oldLeadBase = leadBase;
        }
    } 

  G4LogicalVolume* logicLeadBase = new G4LogicalVolume(leadBase, shieldPbMat, "logicLeadBase",0,0,0);
  new G4PVPlacement(0,leadPosition,"physicLeadBase",logicLeadBase,physicalWorld,false,0);
  logicLeadBase->SetVisAttributes(leadVis);

  // square side panels
  // initial setup
  
  G4Box* squareLeadBase = new G4Box("squareLeadBase",21*2.54*cm,0.25*2.54*cm,27.5*2.54*cm);
  
  // place squares
  G4LogicalVolume* logicLeadSquare = new G4LogicalVolume(squareLeadBase,shieldPbMat,"logicLeadSquare",0,0,0);
  leadPosition=G4ThreeVector(frame_x+10*2.54*cm,frame_y-(11.25)*2.54*cm,frame_z);
  new G4PVPlacement(0,leadPosition,"physicLeadSquare",logicLeadSquare,physicalWorld,false,0);
  leadPosition=G4ThreeVector(frame_x+10*2.54*cm,frame_y+(35.25)*2.54*cm,frame_z);
  new G4PVPlacement(0,leadPosition,"physicLeadSquare1",logicLeadSquare,physicalWorld,false,1);
  logicLeadSquare->SetVisAttributes(leadVis);


  // rect side panels
  // initial setup

  G4Box* rectLeadBase = new G4Box("rectLeadBase",0.25*2.54*cm,23*2.54*cm,27.5*2.54*cm);
  
  //place squares
  G4LogicalVolume* logicLeadRect = new G4LogicalVolume(rectLeadBase,shieldPbMat,"logicLeadRect",0,0,0);
  leadPosition=G4ThreeVector(frame_x-(11.25)*2.54*cm,frame_y+12*2.54*cm,frame_z);
  new G4PVPlacement(0,leadPosition,"physicLeadRect",logicLeadRect,physicalWorld,false,0);
  leadPosition=G4ThreeVector(frame_x+(31.25)*2.54*cm,frame_y+12*2.54*cm,frame_z);
  new G4PVPlacement(0,leadPosition,"physicLeadRect1",logicLeadRect,physicalWorld,false,1);
  logicLeadRect->SetVisAttributes(leadVis);

                          
} // ends Shield Construction


// ------------------------------------------------
// ------------------------------------------------

void k100_DetectorConstruction::ConstructIceBox(G4LogicalVolume*  logicalWorld)
{
  // THIS AREA FOR THE K100 FRIDGE ITSELF
  // just to test...some zero crosshairs:
  /*G4ThreeVector crosspos(0,0,0);
  G4Box* xcross = new G4Box("xcross", 10.*cm, .01*cm, .01*cm);
  G4LogicalVolume* logicross = new G4LogicalVolume(xcross,steel,"logicross",0,0,0);

  new G4PVPlacement(0,crosspos,"physicross",logicross,physicalWorld,false,0);

  G4RotationMatrix crot;
  crot.rotateZ(90.*deg);
  G4Transform3D transcross(crot, crosspos);
  new G4PVPlacement(transcross, "physicross2",logicross,physicalWorld,false,1);

  crot.rotateX(90.*deg);
  transcross = G4Transform3D(crot, crosspos);
  new G4PVPlacement(transcross, "physicross3",logicross,physicalWorld,false,2);*/
  
  // --------------------------- visuals ------------------------------
  G4VisAttributes* copperVis = new G4VisAttributes(G4Colour(184/255.,115/255.,51/255.));
  copperVis->SetForceSolid(false);
  G4VisAttributes* brassVis = new G4VisAttributes(G4Colour(255/255.,215/255.,0/255.));
  brassVis->SetForceSolid(false);
  G4VisAttributes* steelVis = new G4VisAttributes(G4Colour(105/255.,105/255.,105/255.));
  steelVis->SetForceSolid(false);
  G4VisAttributes* heliumVis = new G4VisAttributes(G4Colour(50/255.,205/255.,50/255.));
  heliumVis->SetForceSolid(false);
  G4VisAttributes* superVis = new G4VisAttributes(G4Colour(0/255.,0/255.,255/255.));
  superVis->SetForceSolid(false);

  
  // ------------------------- 30 mK shield ---------------------------
  G4ThreeVector mK30Pos(fridge_x,fridge_y,fridge_z);
  G4Tubs* mK30Solid = new G4Tubs("mK30Solid", .5*4.438*2.54*cm,.5*4.5*2.54*cm, .5*9.75*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicmK30 = new G4LogicalVolume(mK30Solid,towerMat,"logicmK30",0,0,0);
  new G4PVPlacement(0,mK30Pos,"physicmK30",logicmK30,physicalWorld,false,0);
  logicmK30->SetVisAttributes(copperVis);

  // upper endcap
  G4ThreeVector mK30topPos(fridge_x,fridge_y,fridge_z+5.*2.54*cm);
  G4Tubs* mK30top = new G4Tubs("mK30top", 0, .5*4.5*2.54*cm, .5*.25*2.54*cm,0,2*pi);
  G4LogicalVolume* logicmK30top = new G4LogicalVolume(mK30top,towerMat,"logicmK30top",0,0,0);
  new G4PVPlacement(0,mK30topPos,"physicmK30top",logicmK30top,physicalWorld,false,0);
  logicmK30top->SetVisAttributes(copperVis);

  // lower endcap
  G4ThreeVector mK30lowPos(fridge_x,fridge_y,fridge_z-5.*2.54*cm);
  G4Tubs* mK30lowbase = new G4Tubs("mK30lowbase",0,.5*4.5*2.54*cm,.5*.25*2.54*cm,0,2*pi);
  G4Polyhedra* towerHole = new G4Polyhedra("towerHole",0.*deg,360.*deg,6,Tower_nZcut,Tower_zPcut,Tower_rIcut,Tower_rOcut);
  G4SubtractionSolid* mK30low = new G4SubtractionSolid("mK30low",mK30lowbase,towerHole);
  G4LogicalVolume* logicmK30low=new G4LogicalVolume(mK30low,towerMat,"logicmK30low",0,0,0);
  new G4PVPlacement(0,mK30lowPos,"physicmK30low",logicmK30low,physicalWorld,false,0);
  logicmK30low->SetVisAttributes(copperVis);
  

  // ----------------------- 100 mK shield -----------------------
  G4ThreeVector mK100pos(fridge_x,fridge_y,fridge_z-.2935*2.54*cm);
  G4Tubs* mK100 = new G4Tubs("mK100", .5*8.256*2.54*cm, .5*8.320*2.54*cm, .5*11.263*2.54*cm,0,2*pi);
  G4LogicalVolume* logicmK100 = new G4LogicalVolume(mK100,towerMat,"logicmK100",0,0,0);
  new G4PVPlacement(0,mK100pos,"physicmK100",logicmK100,physicalWorld,false,0);
  logicmK100->SetVisAttributes(copperVis);

  // upper endcap
  G4ThreeVector mK100topPos(fridge_x,fridge_y,fridge_z+5.463*2.54*cm);
  G4Tubs* mK100top = new G4Tubs("mK100top", 0, .5*8.320*2.54*cm, .5*.25*2.54*cm,0,2*pi);
  G4LogicalVolume* logicmK100top = new G4LogicalVolume(mK100top,towerMat,"logicmK100top",0,0,0);
  new G4PVPlacement(0,mK100topPos,"physicmK100top",logicmK100top,physicalWorld,false,0);
  logicmK100top->SetVisAttributes(copperVis);

  // lower endcap
  G4ThreeVector mK100lowPos(fridge_x,fridge_y,fridge_z-6.175*2.54*cm);
  G4Tubs* mK100lowbase = new G4Tubs("mK100lowbase",0,.5*8.320*2.54*cm, .5*.5*2.54*cm,0,2*pi);
  G4SubtractionSolid* mK100low = new G4SubtractionSolid("mK100low",mK100lowbase,towerHole);
  G4LogicalVolume* logicmK100low = new G4LogicalVolume(mK100low,brass,"logicmK100low",0,0,0);
  new G4PVPlacement(0,mK100lowPos,"physicmK100low",logicmK100low,physicalWorld,false,0);
  logicmK100low->SetVisAttributes(brassVis);
  
  // -------------------- upper copper and brass pieces ---------------------
  // lower endcap
  G4ThreeVector uBrLowPos(fridge_x,fridge_y,fridge_z+5.838*2.54*cm);
  G4Tubs* uBrLow = new G4Tubs("uBrLow", 0, .5*3.5*2.54*cm, .5*.5*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicuBrLow = new G4LogicalVolume(uBrLow,brass,"logicuBrLow",0,0,0);
  new G4PVPlacement(0,uBrLowPos,"physicuBrLow",logicuBrLow,physicalWorld,false,0);
  logicuBrLow->SetVisAttributes(brassVis); 

  // copper cylinder
  G4ThreeVector uCuPos(fridge_x,fridge_y,fridge_x+9.06*2.54*cm);
  G4Tubs* uCu = new G4Tubs("uCu",.5*3.375*2.54*cm,.5*3.5*2.54*cm,.5*5.944*2.54*cm,0,2*pi);
  G4LogicalVolume* logicuCu = new G4LogicalVolume(uCu,towerMat,"logicuCu",0,0,0);
  new G4PVPlacement(0,uCuPos,"physicuCu",logicuCu,physicalWorld,false,0);
  logicuCu->SetVisAttributes(copperVis);

  // upper endcap
  G4ThreeVector uBrTopPos(fridge_x,fridge_y,fridge_z+12.407*2.54*cm);
  G4Tubs* uBrTop = new G4Tubs("uBrTop", 0, .5*3.5*2.54*cm, .5*.75*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicuBrTop = new G4LogicalVolume(uBrTop,brass,"logicuBrTop",0,0,0);
  new G4PVPlacement(0,uBrTopPos,"physicuBrTop",logicuBrTop,physicalWorld,false,0);
  logicuBrTop->SetVisAttributes(brassVis);

  // -------------------- inner vacuum can ---------------------------
  G4ThreeVector ivcPos(fridge_x,fridge_y,fridge_z-1.3035*2.54*cm);
  G4Tubs* ivc = new G4Tubs("ivc", .5*8.87*2.54*cm, .5*9.*2.54*cm, .5*13.483*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicivc = new G4LogicalVolume(ivc,steel,"logicivc",0,0,0);
  new G4PVPlacement(0,ivcPos,"physicivc",logicivc,physicalWorld,false,0);
  logicivc->SetVisAttributes(steelVis); 
  
  // upper endcap
  G4ThreeVector ivcTopPos(fridge_x,fridge_y,fridge_z+5.888*2.54*cm);
  G4Tubs* ivcTop = new G4Tubs("ivcTop", .5*4.*2.54*cm, .5*9.*2.54*cm, .5*.9*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicivcTop = new G4LogicalVolume(ivcTop,steel,"logicivcTop",0,0,0);
  new G4PVPlacement(0,ivcTopPos,"physicivcTop",logicivcTop,physicalWorld,false,0);
  logicivcTop->SetVisAttributes(steelVis); 

  // upper pipe
  G4ThreeVector uppPos(fridge_x,fridge_y,fridge_z+14.191*2.54*cm);
  G4Tubs* upp = new G4Tubs("upp", .5*3.834*2.54*cm, .5*4.*2.54*cm, .5*17.506*2.54*cm,0,2*pi);
  G4LogicalVolume* logicupp = new G4LogicalVolume(upp,steel,"logicupp",0,0,0);
  new G4PVPlacement(0,uppPos,"physicupp",logicupp,physicalWorld,false,0);
  logicupp->SetVisAttributes(steelVis);

  // upper pipe cap
  G4ThreeVector upcPos(fridge_x,fridge_y,fridge_z+23.144*2.54*cm);
  G4Tubs* upc = new G4Tubs("upc", 0,.5*4.*2.54*cm,.5*.4*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicupc = new G4LogicalVolume(upc,steel,"logicupc",0,0,0);
  new G4PVPlacement(0,upcPos,"physicupc",logicupc,physicalWorld,false,0);
  logicupc->SetVisAttributes(steelVis);

  // lower endcap
  G4ThreeVector ivcLowPos(fridge_x,fridge_y,fridge_z-8.495*2.54*cm);
  G4Tubs* ivcLowBase = new G4Tubs("ivcLowBase", 0, .5*9.*2.54*cm, .5*.9*2.54*cm,0,2*pi);
  G4SubtractionSolid* ivcLow = new G4SubtractionSolid("ivcLow",ivcLowBase,towerHole);
  G4LogicalVolume* logicivcLow = new G4LogicalVolume(ivcLow,steel,"logicivcLow",0,0,0);
  new G4PVPlacement(0,ivcLowPos,"physicivcLow",logicivcLow,physicalWorld,false,0);
  logicivcLow->SetVisAttributes(steelVis); 

  // electronics box
  G4ThreeVector eboxPos(fridge_x,fridge_y,fridge_z-10.785*2.54*cm);
  G4Tubs* ebox = new G4Tubs("ebox", .5*7.87*2.54*cm, .5*8.*2.54*cm, .5*3.68*2.54*cm,0,2*pi);
  G4LogicalVolume* logicebox = new G4LogicalVolume(ebox,steel,"logicebox",0,0,0);
  new G4PVPlacement(0,eboxPos,"physicebox",logicebox,physicalWorld,false,0);
  logicebox->SetVisAttributes(steelVis);

  // lower cap for electronics box
  G4ThreeVector ecapPos(fridge_x,fridge_y,fridge_z-12.75*2.54*cm);
  G4Tubs* ecap = new G4Tubs("ecap", 0, .5*8.*2.54*cm, .5*.25*2.54*cm,0,2*pi);
  G4LogicalVolume* logicecap=new G4LogicalVolume(ecap,steel,"logicecap",0,0,0);
  new G4PVPlacement(0,ecapPos,"physicebox",logicecap,physicalWorld,false,0);
  logicecap->SetVisAttributes(steelVis);

  // ---------------------------- Liquid Helium --------------------------
  // top piece
  G4ThreeVector topHePos(fridge_x,fridge_y,fridge_z+35.9075*2.54*cm);
  G4Tubs* topHe = new G4Tubs("topHe", 0, .5*10.*2.54*cm, .5*18.563*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logictopHe = new G4LogicalVolume(topHe,helium,"logictopHe",0,0,0);
  new G4PVPlacement(0,topHePos,"physictopHe",logictopHe,physicalWorld,false,0);
  logictopHe->SetVisAttributes(heliumVis);

  // top center piece
  G4ThreeVector tcHePos(fridge_x,fridge_y,fridge_z+24.985*2.54*cm);
  G4Tubs* tcHe = new G4Tubs("tcHe", 0, .5*14.*2.54*cm, .5*3.282*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logictcHe = new G4LogicalVolume(tcHe,helium,"logictcHe",0,0,0);
  new G4PVPlacement(0,tcHePos,"physictcHe",logictcHe,physicalWorld,false,0);
  logictcHe->SetVisAttributes(heliumVis);

  // center piece
  G4ThreeVector ceHePos(fridge_x,fridge_y,fridge_z+14.841*2.54*cm);
  G4Tubs* ceHe = new G4Tubs("ceHe", .5*4.*2.54*cm, .5*14.*2.54*cm, .5*17.006*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logiceHe = new G4LogicalVolume(ceHe,helium,"logiceHe",0,0,0);
  new G4PVPlacement(0,ceHePos,"physiceHe",logiceHe,physicalWorld,false,0);
  logiceHe->SetVisAttributes(heliumVis);
  
  // lower center piece
  G4ThreeVector lcHePos(fridge_x,fridge_y,fridge_z+5.8255*2.54*cm);
  G4Tubs* lcHe = new G4Tubs("lcHe", .5*9.*2.54*cm, .5*14.*2.54*cm, .5*1.025*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logiclcHe = new G4LogicalVolume(lcHe,helium,"logiclcHe",0,0,0);
  new G4PVPlacement(0,lcHePos,"physiclcHe",logiclcHe,physicalWorld,0,0,0);
  logiclcHe->SetVisAttributes(heliumVis);
  
  // upper lower piece
  G4ThreeVector ulHePos(fridge_x,fridge_y,fridge_z-1.816*2.54*cm);
  G4Tubs* ulHe = new G4Tubs("ulHe", .5*9.*2.54*cm, .5*10.*2.54*cm, .5*14.258*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logiculHe = new G4LogicalVolume(ulHe,helium,"logiculHe",0,0,0);
  new G4PVPlacement(0,ulHePos,"physiculHe",logiculHe,physicalWorld,0,0,0);
  logiculHe->SetVisAttributes(heliumVis);

  // center lower piece
  G4ThreeVector clHePos(fridge_x,fridge_y,fridge_z-10.91*2.54*cm);
  G4Tubs* clHe = new G4Tubs("clHe", .5*8.*2.54*cm, .5*10.*2.54*cm, .5*3.93*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logiclHe = new G4LogicalVolume(clHe,helium,"logiclHe",0,0,0);
  new G4PVPlacement(0,clHePos,"physiclHe",logiclHe,physicalWorld,0,0,0);
  logiclHe->SetVisAttributes(heliumVis);

  // bottom piece
  G4ThreeVector botHePos(fridge_x,fridge_y,fridge_z-13.*2.54*cm);
  G4Tubs* botHe = new G4Tubs("botHe", 0, .5*10.*2.54*cm, .5*.25*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicbotHe = new G4LogicalVolume(botHe,helium,"logicbotHe",0,0,0);
  new G4PVPlacement(0,botHePos,"physicbotHe",logicbotHe,physicalWorld,0,0,0);
  logicbotHe->SetVisAttributes(heliumVis);

  // ------------------------- Shell ----------------------
  // Outer Edge
  G4ThreeVector oeShellPos(fridge_x,fridge_y,fridge_z+12.9045*2.54*cm);
  G4Tubs* oeShell = new G4Tubs("oeShell", .5*21.5*2.54*cm, .5*22.*2.54*cm, .5*64.067*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicoeShell = new G4LogicalVolume(oeShell,steel,"logicoeShell",0,0,0);
  new G4PVPlacement(0,oeShellPos,"physicoeShell",logicoeShell,physicalWorld,false,0);
  logicoeShell->SetVisAttributes(steelVis);

  // Top Edge
  G4ThreeVector teShellPos(fridge_x,fridge_y,fridge_z+45.063*2.54*cm);
  G4Tubs* teShell = new G4Tubs("teShell", .5*10.*2.54*cm, .5*22.*2.54*cm, .5*.25*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicteShell = new G4LogicalVolume(teShell,steel,"logicteShell",0,0,0);
  new G4PVPlacement(0,teShellPos,"physicteShell",logicteShell,physicalWorld,false,0);
  logicteShell->SetVisAttributes(steelVis);

  // Bottom Edge
  G4ThreeVector btShellPos(fridge_x,fridge_y,fridge_z-19.254*2.54*cm);
  G4Tubs* btShell = new G4Tubs("btShell", 0, .5*22.*2.54*cm, .5*.25*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicbtShell = new G4LogicalVolume(btShell,steel,"logicbtShell",0,0,0);
  new G4PVPlacement(0,btShellPos,"physicbtShell",logicbtShell,physicalWorld,false,0);
  logicbtShell->SetVisAttributes(steelVis);

  // Upper Inner (core)
  G4ThreeVector ucShellPos(fridge_x,fridge_y,fridge_z+35.9075*2.54*cm);
  G4Tubs* ucShell = new G4Tubs("ucShell", .5*10.*2.54*cm, .5*10.5*2.54*cm, .5*18.063*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicucShell = new G4LogicalVolume(ucShell,steel,"logicucShell",0,0,0);
  new G4PVPlacement(0,ucShellPos,"physicucShell",logicucShell,physicalWorld,false,0);
  logicucShell->SetVisAttributes(steelVis);

  // Upper Interface
  G4ThreeVector uiShellPos(fridge_x,fridge_y,fridge_z+26.751*2.54*cm);
  G4Tubs* uiShell = new G4Tubs("uiShell", .5*10.*2.54*cm, .5*14.5*2.54*cm, .5*.25*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicuiShell = new G4LogicalVolume(uiShell,steel,"logicuiShell",0,0,0);
  new G4PVPlacement(0,uiShellPos,"physicuiShell",logicuiShell,physicalWorld,false,0);
  logicuiShell->SetVisAttributes(steelVis);

  // Middle Inner (core)
  G4ThreeVector mcShellPos(fridge_x,fridge_y,fridge_z+15.9695*2.54*cm);
  G4Tubs* mcShell = new G4Tubs("mcShell", .5*14.*2.54*cm, .5*14.5*2.54*cm, .5*21.313*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicmcShell = new G4LogicalVolume(mcShell,steel,"logicmcShell",0,0,0);
  new G4PVPlacement(0,mcShellPos,"physicmcShell",logicmcShell,physicalWorld,false,0);
  logicmcShell->SetVisAttributes(steelVis);

  // Lower Interface
  G4ThreeVector liShellPos(fridge_x,fridge_y,fridge_z+5.188*2.54*cm);
  G4Tubs* liShell = new G4Tubs("liShell", .5*10.*2.54*cm, .5*14.5*2.54*cm, .5*.25*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicliShell = new G4LogicalVolume(liShell,steel,"logicliShell",0,0,0);
  new G4PVPlacement(0,liShellPos,"physicliShell",logicliShell,physicalWorld,false,0);
  logicliShell->SetVisAttributes(steelVis);

  // Lower Inner (core)
  G4ThreeVector lcShellPos(fridge_x,fridge_y,fridge_z-4.031*2.54*cm);
  G4Tubs* lcShell = new G4Tubs("lcShell", .5*10.*2.54*cm, .5*10.5*2.54*cm, .5*18.188*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logiclcShell = new G4LogicalVolume(lcShell,steel,"logiclcShell",0,0,0);
  new G4PVPlacement(0,lcShellPos,"physiclcShell",logiclcShell,physicalWorld,false,0);
  logiclcShell->SetVisAttributes(steelVis);

  // Lower Cap
  G4ThreeVector lpShellPos(fridge_x,fridge_y,fridge_z-13.25*2.54*cm);
  G4Tubs* lpShell = new G4Tubs("lpShell", 0, .5*10.5*2.54*cm, .5*.25*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logiclpShell = new G4LogicalVolume(lpShell,steel,"logiclpShell",0,0,0);
  new G4PVPlacement(0,lpShellPos,"physiclpShell",logiclpShell,physicalWorld,false,0);
  logiclpShell->SetVisAttributes(steelVis);

  // ------------------------ Insulation --------------------
  // Upper piece
  G4ThreeVector upInPos(fridge_x,fridge_y,fridge_z+35.9075*2.54*cm);
  G4Tubs* upIn = new G4Tubs("upIn", .5*10.5*2.54*cm, .5*21.5*2.54*cm, .5*18.063*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicupIn = new G4LogicalVolume(upIn,super,"logicupIn",0,0,0);
  new G4PVPlacement(0,upInPos,"physicupIn",logicupIn,physicalWorld,false,0);
  logicupIn->SetVisAttributes(superVis);

  // Middle piece
  G4ThreeVector miInPos(fridge_x,fridge_y,fridge_z+15.9695*2.54*cm);
  G4Tubs* miIn = new G4Tubs("miIn", .5*14.5*2.54*cm, .5*21.5*2.54*cm, .5*21.813*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicmiIn = new G4LogicalVolume(miIn,super,"logicmiIn",0,0,0);
  new G4PVPlacement(0,miInPos,"physicmiIn",logicmiIn,physicalWorld,false,0);
  logicmiIn->SetVisAttributes(superVis);

  // Lower piece
  G4ThreeVector loInPos(fridge_x,fridge_y,fridge_z-4.156*2.54*cm);
  G4Tubs* loIn = new G4Tubs("loIn", .5*10.5*2.54*cm, .5*21.5*2.54*cm, .5*18.438*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicloIn = new G4LogicalVolume(loIn,super,"logicloIn",0,0,0);
  new G4PVPlacement(0,loInPos,"physicloIn",logicloIn,physicalWorld,false,0);
  logicloIn->SetVisAttributes(superVis);

  // Bottom piece
  G4ThreeVector btInPos(fridge_x,fridge_y,fridge_z-16.252*2.54*cm);
  G4Tubs* btIn = new G4Tubs("btIn", 0, .5*21.5*2.54*cm, .5*5.754*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicbtIn = new G4LogicalVolume(btIn,super,"logicbtIn",0,0,0);
  new G4PVPlacement(0,btInPos,"physicbtIn",logicbtIn,physicalWorld,false,0);
  logicbtIn->SetVisAttributes(superVis);
} // ends IceBox Construction

// ---------------- end of class ---------------
void k100_DetectorConstruction::ConstructThermalNeutronBox(G4VPhysicalVolume *world)
{
        //displacement in world
	G4ThreeVector disp = G4ThreeVector(0.0,0.0,-1200); //2m down

	//define dimensions
	G4double leadThick = 50.8;
	G4double polyD = 500.0;
	G4double polyThick = 500.0;

	//create a lead cylinder
	G4Tubs* leadCylinder = new G4Tubs("leadCyl_S",0.0,polyD/2+leadThick,polyThick/2.0+leadThick/2.0,0.0,2*pi);
	G4LogicalVolume* logicalLeadCylinder;
        logicalLeadCylinder = new G4LogicalVolume(leadCylinder,shieldPbMat,"polyCyl_L",0,0,0);
	G4VPhysicalVolume* cylinderLeadWorld = new G4PVPlacement(0, 
								disp,
								"leadCyl_P",
								logicalLeadCylinder,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttLeadCyl = new G4VisAttributes(G4Colour(3.,3.,3.));
	VisAttLeadCyl->SetForceWireframe(false);  //I want a Wireframe of the me
	logicalLeadCylinder->SetVisAttributes(VisAttLeadCyl);  
	// Make Invisible
	//logicalLeadCylinder->SetVisAttributes(G4VisAttributes::Invisible);
	
        //create a poly cylinder
	G4ThreeVector polydisp = G4ThreeVector(0.0,0.0,leadThick/2.0);
	G4Tubs* polyCylinder = new G4Tubs("polyCyl_S",0.0,polyD/2,polyThick/2.0,0.0,2*pi);
	G4LogicalVolume* logicalPolyCylinder;
        logicalPolyCylinder = new G4LogicalVolume(polyCylinder,polyMat,"polyCyl_L",0,0,0);
	G4VPhysicalVolume* cylinderPolyWorld = new G4PVPlacement(0, 
								polydisp,
								"polyCyl_P",
								logicalPolyCylinder,
								cylinderLeadWorld,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttPolyCyl = new G4VisAttributes(G4Colour(7.,3.,0.));
	VisAttPolyCyl->SetForceWireframe(false);  //I want a Wireframe of the me
	logicalPolyCylinder->SetVisAttributes(VisAttPolyCyl);  
	// Make Invisible
	//logicalPolyCylinder->SetVisAttributes(G4VisAttributes::Invisible);

        return;
}
void k100_DetectorConstruction::ConstructShieldTestEnvironment(G4VPhysicalVolume *world)
{
        //Get the coordinates of copies and source point and things
	G4ThreeVector stackrel_detorigin = zipParam->GetCoordinates(0);  //is the first one is 0
	G4double towerassy_shift = -6.5*cm;
	G4double zipstack_shift = -(Tower_zPcut[0] - (zPldh[1] - zPldh[0]) - zPsdh[1] -Zip_z/2.0); //negative because tower is flipped
	G4double instack_shift = -stackrel_detorigin.z(); //negative because tower is flipped
	G4ThreeVector detorigin = G4ThreeVector(0,0,towerassy_shift + zipstack_shift +instack_shift);
	G4ThreeVector point(shieldTestParams.xcntr,shieldTestParams.ycntr,shieldTestParams.zcntr);
	G4ThreeVector relative = point - detorigin;

	G4cout << "xdet: " << detorigin.x() << " ydet: " << detorigin.y() << " zdet: " << detorigin.z() << G4endl;
	G4cout << "towerassy_shift: " << towerassy_shift << " zipstack_shift: " << zipstack_shift << " instack_shift: " << instack_shift << G4endl;

	//Do the appropriate rotations
	G4RotationMatrix *shieldrot = new G4RotationMatrix;
	shieldrot->rotateZ(-relative.phi());
	shieldrot->rotateY(-relative.theta());

	//create the simple rectangular shield
	G4Box* shieldBox = new G4Box("shieldBox_S",shieldTestParams.sizel/2.0,shieldTestParams.sizew/2.0,shieldTestParams.sizethk/2.0);
	G4LogicalVolume* logicalShieldBox;
        logicalShieldBox = new G4LogicalVolume(shieldBox,shieldTestParams.shieldmaterial,"shieldBox_L",0,0,0);
	G4VPhysicalVolume* shieldBoxWorld = new G4PVPlacement(shieldrot, 
								point - (relative.unit()*shieldTestParams.sizethk/2.0), //put source at edge of shielding always
								"shieldBox_P",
								logicalShieldBox,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttShieldBox = new G4VisAttributes(G4Colour(7.,3.,0.));
	VisAttShieldBox->SetForceWireframe(false);  //I want a Wireframe of the me
	logicalShieldBox->SetVisAttributes(VisAttShieldBox);  
}
