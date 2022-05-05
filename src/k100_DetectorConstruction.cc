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
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
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
#include "k100_StdSD.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

//spandey
#include "k100_NaIParametrization.hh"

#include <vector>
#include <ostream>

// ------------------------------------------------

k100_DetectorConstruction::k100_DetectorConstruction()
{

  // -------- The World ---------
  //world_x = 250.*cm; world_y = 250.*cm; world_z = 250.*cm;
  //contrive the overall dimensions to fit the walls and stuff
  //world_x = 250.*cm; world_y = 350.*cm; world_z = 700.*cm;
  world_y = 160.*cm + 5.5*2.54*cm; world_x = 188.*cm + 5.5*2.54*cm; world_z = 1100.*cm;
  world_y*=2; //above are x and y half widths, but this should be full width
  world_x*=2; 

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
  ConstructShieldsBool = false;
  ConstructIceBoxBool = false;
  ConstructFloorBool = false;
  ConstructWallsBool = false;
  ConstructCeilingBool = false;
  ConstructWestReflectorBool = false;
  ConstructFrameBool = false;
  ConstructPuBeSourceAndShieldBool = false;
  ConstructNaIArrayBool = false; //spandey
  SetConstructThermalNeutronBoxBool(false); //note, requires construct ZIP bool
  SetConstructShieldTestEnvironmentBool(false); //note, requires construct ZIP bool
  SetConstructSimpleGammaCoinBool(false); //note, requires construct ZIP bool
  SetConstructPuBeNaIBool(false); //note, requires construct ZIP bool
  SetFirstDetGe(true); //we only use the first detector in the array for now, should it be Ge?


  //
  DrawSolidDetBox = true; DrawSolidZipBool = true;
  DrawSolidTowerBool = true; 
  DrawSolidVetoBool = false;
  DrawSolidShieldsBool = false; DrawSolidIceBoxBool = false;


  // ---------Material Definition--------------
  DefineMaterials();

  //more complicated parameters
  shieldTestParams.xcntr = 100.0*cm;
  shieldTestParams.ycntr = 0.0*cm;
  shieldTestParams.zcntr = 0.0*cm;
  shieldTestParams.sizel = 10.0*cm;
  shieldTestParams.sizew = 10.0*cm;
  shieldTestParams.sizethk = 10.0*cm;
  shieldTestParams.shieldmaterial = polyMat; //MUST be after DefineMaterials()

  gammaCoinParams.xcntr = -100.0*cm; //start off opposed to shield test
  gammaCoinParams.ycntr = 0.0*cm;
  gammaCoinParams.zcntr = 0.0*cm;
  gammaCoinParams.sizer = 2.54*3*cm;
  gammaCoinParams.sizethk = 5.0*cm;
  gammaCoinParams.coinmaterial = zipGeMat; //MUST be after DefineMaterials()

  frameParams.includeSand = false;

  fridgeParams.includeMixture = false;
  fridgeParams.pure3HeBath = false;

  shieldParams.addNaISouth = false;
  shieldParams.HPGeboron = false;
  shieldParams.HPGeboron_wshield = false;
  shieldParams.addBasePoly = false;
  shieldParams.addBaseLead = false;
  shieldParams.mod = 0; 
 
  pubeNaIParams.westPolySensitivity = false; //default is not to make this sensitive 
  pubeNaIParams.doPuBeGamma = true; //default to do PuBe gammas
  pubeNaIParams.addBarrel = true; //default to use barrel
  pubeNaIParams.NaIsensitivity = true; //default to have NaI sensitive 
  pubeNaIParams.doR66 = true; //default to R66 shield 
  pubeNaIParams.doR62 = false; //FIXME not yet implemented
  pubeNaIParams.mod = 0; 
  pubeNaIParams.doOrb = false; 
  pubeNaIParams.OrbPos = G4ThreeVector(0,0,0);
  pubeNaIParams.OrbRad = 10*cm;

  // ---------Detector Names--------------
  DetCollName = new char*[30];  TowCollName = new char*[5];     DetMaterials = new G4int [30];
  DetCollName[0]  = "zip01";  
  if(FirstDetGe)
    DetMaterials[0]  = 1; //Ge  ///Ge = 1 Si = 0
  else
    DetMaterials[0]  = 0; //Si  ///Ge = 1 Si = 0
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

  // Define Deuterium 
  G4Element* elementD=new G4Element(name="Deuterium", symbol="D", z=1., a=2.014*g/mole);

  // Define Helium
  G4Element* elementHe = new G4Element(name="Helium", symbol="He", z=2., a=4.003*g/mole);
  G4IsotopeVector *ivec;
  G4double *relabvec;
  ivec = elementHe->GetIsotopeVector();
  relabvec = elementHe->GetRelativeAbundanceVector();
  for(int i=0;i<elementHe->GetNumberOfIsotopes();i++){
    G4cout << "He isotope " << i+1 << " is A = " << (*ivec)[i]->GetN() <<
	    ", relative abundance: " << relabvec[i] << G4endl;
  }
  G4cout << "Is the He natural abundance distribution? " << elementHe->GetNaturalAbundanceFlag() << G4endl;

  // Define specialized Helium
  G4Isotope *isoHe4 = new G4Isotope("4He",2,4); //when it asks for N does it really mean (integer) A?
  G4Isotope *isoHe3 = new G4Isotope("3He",2,3); 
  G4Element *stillLiquidHe = new G4Element("Liquid3He","L3He",1);
  stillLiquidHe->AddIsotope(isoHe3,1.0);
  G4Element *MCLiquidHe = new G4Element("MCLiquidHe","LHeMix",2);
  MCLiquidHe->AddIsotope(isoHe3,0.12);
  MCLiquidHe->AddIsotope(isoHe4,0.88);

  // Define Boron 
  G4Element* elementB=new G4Element(name="Boron", symbol="B", z=5., a=10.811*g/mole);

  // Define Carbon
  G4Element* elementC=new G4Element(name="Carbon", symbol="C", z=6., a=12.011*g/mole);

  // Define Oxygen 
  G4Element* elementO=new G4Element(name="Oxygen", symbol="O", z=8., a=15.9994*g/mole);

  // Define Sodium
  G4Element* elementNa = new G4Element(name="Sodium", symbol="Na", z=11., a=22.9898*g/mole);

  // Define Silicon
  G4Element* elementSi = new G4Element(name="Silicon", symbol="Si", z=14., a=28.09*g/mole);

  // Define Iron
  G4Element* elementFe = new G4Element(name="Iron", symbol="Fe", z=26., a=55.845*g/mole);

  // Define Copper
  G4Element* elementCu = new G4Element(name="Copper", symbol="Cu", z=29., a=63.5460*g/mole);

  // Define Germanium   
  G4Element* elementGe = new G4Element(name="Germanium", symbol="Ge", z=32., a=72.61*g/mole);

  // Define Iodine   
  G4Element* elementI = new G4Element(name="Iodine", symbol="I", z=53., a=126.90447*g/mole);

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

  // LightAluminum -- i.e. a plate of aluminum with a 2" pitch array of 1/4-20 holes
  G4Material* LightAluminum = new G4Material(name="Aluminum", density = 0.989*2.7*g/cm3, ncomponents=1);
  LightAluminum->AddElement(elementAl, natoms=1);

  //wood
  G4Material* WOOD = new G4Material(name="wood", density=0.9*g/cm3, ncomponents=3);
  WOOD->AddElement(elementH , 4);
  WOOD->AddElement(elementO , 1);
  WOOD->AddElement(elementC , 2);

  //sodium borate anhydrous
  G4Material* sba = new G4Material(name="sodium_borate_anhydrous", density=2.367*g/cm3, ncomponents=3);
  sba->AddElement(elementB , 4);
  sba->AddElement(elementO , 7);
  sba->AddElement(elementNa , 2);

  // Silicon 
  G4Material* Silicon = new G4Material(name="Silicon", density = 2.330*g/cm3, ncomponents=1);
  Silicon->AddElement(elementSi, natoms=1);
  
  // Germanium 
  G4Material* Germanium = new G4Material(name="Germanium", density = 5.323*g/cm3, ncomponents=1);
  Germanium->AddElement(elementGe, natoms=1);

  // Sodium Iodide 
  G4Material* NaI = new G4Material(name="NaI", density = 3.67*g/cm3, ncomponents=2);
  NaI->AddElement(elementNa, natoms=1);
  NaI->AddElement(elementI, natoms=1);

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

  // H2O 
  G4Material* h2o=new G4Material(name="H2O", density = 1.0*g/cm3, ncomponents=2);
  h2o->AddElement(elementH,natoms=2);
  h2o->AddElement(elementO,natoms=1);

  // D2O 
  G4Material* d2o=new G4Material(name="D2O", density = 1.11*g/cm3, ncomponents=2);
  d2o->AddElement(elementD,natoms=2);
  d2o->AddElement(elementO,natoms=1);

  // SiO2 amorphous
  G4Material* sio2=new G4Material(name="SiO2", density = 2.196*g/cm3, ncomponents=2);
  sio2->AddElement(elementO,natoms=2);
  sio2->AddElement(elementSi,natoms=1);

  // SiO2 amorphous (blasting sand with 64% packing)
  G4Material* sio2_blastSand=new G4Material(name="SiO2BlastSand", density = 0.64*2.196*g/cm3, ncomponents=2);
  sio2_blastSand->AddElement(elementO,natoms=2);
  sio2_blastSand->AddElement(elementSi,natoms=1);

  // Poly
  //G4Material* poly=new G4Material(name="Poly", density = 0.935*g/cm3, ncomponents=2);
  //use average of density measures for R66 sim N-MISC-17-003 pg 31
  G4Material* poly=new G4Material(name="Poly", density = 0.973*g/cm3, ncomponents=2); 
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

  // More complex Helium Liquid
  G4Material* stillLiquid = new G4Material(name="stillLiquid",density=0.081*g/cm3,ncomponents=1);
  stillLiquid->AddElement(stillLiquidHe,natoms=2);
  G4Material* MCLiquid = new G4Material(name="MCLiquid",density=0.1412*g/cm3,ncomponents=1); //assumed to be same as LHe
  MCLiquid->AddElement(MCLiquidHe,natoms=2);

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

  // Steel (Standard 1% carbon)
  G4Material* StandardSteel = new G4Material(name="StandardSteel", density=7.85*g/cm3, ncomponents=2);
  StandardSteel->AddMaterial(Carbon, .0100);
  StandardSteel->AddMaterial(Iron, .99);

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
  // Define Database Materials
  // ------------------------------------------------
  G4NistManager* man = G4NistManager::Instance(); 
 
  G4NISTconcrete  = man->FindOrBuildMaterial("G4_CONCRETE");
  G4NISTair  = man->FindOrBuildMaterial("G4_AIR");
  G4NISTNaI  = man->FindOrBuildMaterial("G4_SODIUM_IODIDE");
  G4NISTPVC  = man->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
  G4NISTPE  = man->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4NISTAl  = man->FindOrBuildMaterial("G4_Al");
  G4NISTlucite  = man->FindOrBuildMaterial("G4_PLEXIGLASS");
  G4NISTparaffin  = man->FindOrBuildMaterial("G4_PARAFFIN");
  G4NISTGypsum  = man->FindOrBuildMaterial("G4_GYPSUM");
  G4NISTstainless  = man->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  // ------------------------------------------------
  // end define database materials

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
  d2oMat = d2o;
  sio2Mat = sio2;
  h2oMat = h2o;
  naiMat = NaI;
  blastsand = sio2_blastSand;
  shieldCuMat = Copper;
  shieldPbMat = Lead;
  iceboxCuMat = Copper;
  mumetalMat = mumetal;
  aluminum=Aluminum;
  lightaluminum=LightAluminum;
  sodium_borate_anhydrous=sba;
  wood=WOOD;
  steel=Steel;
  carbonsteel=StandardSteel;
  brass=Brass;
  helium=liquidHelium;
  stillHe=stillLiquid;
  MCHe=MCLiquid;
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
  //FIXME: clean sensitive detectors?
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
    //delete k100CollPoint[it->first]; //trouble! never
  }
  k100CollName.clear();
  k100CollPoint.clear();

  // ------------ Construct the Physical world ---------------

  // Construct the World
  G4Box* solidWorld = new G4Box("world_S", 0.5*world_x, 0.5*world_y, 0.5*world_z);
  G4cout << "WORLD MOTHER VOLUME: " << solidWorld->GetCubicVolume()/1e9 << G4endl;
  G4LogicalVolume*  logicalWorld = new G4LogicalVolume(solidWorld,  // The solid
						       G4NISTair, // Material
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
  if(ConstructSimpleGammaCoinBool)  {ConstructSimpleGammaCoin(physicalWorld);}


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
void k100_DetectorConstruction::SetConstructShieldTestEnvironmentMat(G4String mat)
{
  if(mat=="Lead"){
    shieldTestParams.shieldmaterial = shieldPbMat;
  }
  else if(mat=="Poly"){
    shieldTestParams.shieldmaterial = polyMat;
  }
  else if(mat=="D2O"){
    shieldTestParams.shieldmaterial = d2oMat;
  }
  else if(mat=="H2O"){
    shieldTestParams.shieldmaterial = h2oMat;
  }
  else{  //default to poly material
    shieldTestParams.shieldmaterial = polyMat;
  }

}
void k100_DetectorConstruction::SetConstructSimpleGammaCoinPos(G4double xcntr,G4double ycntr,G4double zcntr)
{
  gammaCoinParams.xcntr = xcntr;
  gammaCoinParams.ycntr = ycntr;
  gammaCoinParams.zcntr = zcntr;

}
void k100_DetectorConstruction::SetConstructSimpleGammaCoinSize(G4double sizer,G4double sizethk)
{
  gammaCoinParams.sizer = sizer;
  gammaCoinParams.sizethk = sizethk;

}
void k100_DetectorConstruction::SetConstructSimpleGammaCoinMat(G4String mat)
{
  if(mat=="NaI"){
    gammaCoinParams.coinmaterial = naiMat;
  }
  else if(mat=="HPGe"){
    gammaCoinParams.coinmaterial = zipGeMat;
  }
  else if(mat=="Scint"){
    gammaCoinParams.coinmaterial = scintMat;
  }
  else{  //default to HPGe 
    gammaCoinParams.coinmaterial = zipGeMat;
  }

}
G4String k100_DetectorConstruction::GetConstructShieldTestEnvironmentMat()
{
  return shieldTestParams.shieldmaterial->GetName();
}
G4String k100_DetectorConstruction::GetConstructSimpleGammaCoinMat()
{
  return gammaCoinParams.coinmaterial->GetName();
}

void k100_DetectorConstruction::ConstructTower(G4VPhysicalVolume* physicalDetectorBox)
{

  //----------------------------------- 
  // Contruct the Tower Logical Volume
  //-----------------------------------
  //some measurements of the floor and fridge
  //G4double fridgeHalfHeightToBottomPlate = (12.9045+13.25+0.25)*2.54*cm;
  G4double fridgeHalfHeightToBottomPlate = (12.9045+19.254+0.25)*2.54*cm; //modified 1/1/18 to get floor height right
  //G4double distanceCenterToFloor = fridgeHalfHeightToBottomPlate + 21.0*2.54*cm;
  G4double distanceCenterToFloor = fridgeHalfHeightToBottomPlate + 21.0*2.54*cm -70.86*mm; //compensate for 70.86mm discrepancy in floor distance 1/1/18
  G4double floorZ = fridge_z+12.9045*2.54*cm - distanceCenterToFloor;

  // Position Vector
  //G4ThreeVector positionTower = G4ThreeVector(xtow[1-1],ytow[1-1],ztow[1-1]);
  G4ThreeVector positionTower = G4ThreeVector(tower_x,tower_y,tower_z);

  G4cout << "Position of Tower 1 in World: (" << tower_x << "," << tower_y << "," << tower_z << ")" << G4endl;
  G4cout << "Tower 1 Height Above Floor: " << tower_z-floorZ  << G4endl;
  G4RotationMatrix r180Rotation;		// flip towers over
  r180Rotation.rotateY(180.*deg);
  G4Transform3D towerflip(r180Rotation, positionTower);
  //r180Rotation->rotateX(M_PI*rad);
  //r180Rotation->rotateZ(M_PI*rad);
  G4cout << "Tower_nZcut: " << Tower_nZcut << G4endl;
  G4cout << "Tower_zPcut: " << Tower_zPcut[0] << "\t" << Tower_zPcut[1] << G4endl;
  G4cout << "Tower_rIcut: " << Tower_rIcut[0] << "\t" << Tower_rIcut[1] << G4endl;
  G4cout << "Tower_rOcut: " << Tower_rOcut[0] << "\t" << Tower_rOcut[1] << G4endl;
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
    //logicalTower1->SetVisAttributes(G4VisAttributes::Invisible);  // Make Invisible
    
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

  //check the null case
  if(NbZipsPerTower<1) return;

  //some measurements of the floor and fridge
  //G4double fridgeHalfHeightToBottomPlate = (12.9045+13.25+0.25)*2.54*cm;
  G4double fridgeHalfHeightToBottomPlate = (12.9045+19.254+0.25)*2.54*cm; //modified 1/1/18 to get floor height right
  //G4double distanceCenterToFloor = fridgeHalfHeightToBottomPlate + 21.0*2.54*cm;
  G4double distanceCenterToFloor = fridgeHalfHeightToBottomPlate + 21.0*2.54*cm -70.86*mm; //compensate for 70.86mm discrepancy in floor distance 1/1/18
  G4double floorZ = fridge_z+12.9045*2.54*cm - distanceCenterToFloor;

  G4double zHeightAboveFloor = tower_z-floorZ;

  //set up some variables
  //G4ThreeVector position_ldh = G4ThreeVector(0,0, Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu - zHrb - (zPudh[1]-zPudh[0]) - (zPsdh[1]-zPsdh[0]) - zPldh[1] ); 
  G4double voidThk = NbZipsPerTower*Zip_z + (NbZipsPerTower-1)*Zip_Space + (zPsdh[1]-zPsdh[0]) - Zip_z; //total height plus tiny extra space 
  G4double z0 = Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu - zHrb - (zPudh[1]-zPudh[0]); //z position after top housing of position 1 

  //G4ThreeVector positionZipArray = G4ThreeVector(0,0,Tower_zPcut[0]-(zPldh[1]-zPldh[0])-zPsdh[1]-Zip_z/2);
  //G4ThreeVector positionZipArray = G4ThreeVector(0,0, Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu - zHrb - (zPudh[1]-zPudh[0]) - zPsdh[1]); 
  G4ThreeVector positionZipArray = G4ThreeVector(0,0, z0-voidThk/2.0); 
  //G4cout << "Zip Array Position In Tower: " << Tower_zPcut[0]-(zPldh[1]-zPldh[0])-zPsdh[1]-Zip_z/2 << G4endl;
  //G4cout << "Zip Array Position In Tower: " << Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu - zHrb - (zPudh[1]-zPudh[0]) - zPsdh[1] << G4endl; 
  G4cout << "Zip Array Position In Tower: " << z0-voidThk/2.0 << G4endl;
  G4cout << "Zip Height Above Floor: " << tower_z-floorZ-z0+Zip_z/2.0+(zPsdh[1]-zPsdh[0]-Zip_z)/2.0  << G4endl;
 //Corrected from: G4ThreeVector(0,0,Tower_zPcut[0]+(zPldh[1]-zPldh[0]) + zPsdh[1]);
  //G4Tubs* solidZipArray = new G4Tubs("ZipArray_S", 0.0, Zip_Rout, 6*Zip_Househeight/2,  0, 2*pi);
  G4Tubs* solidZipArray = new G4Tubs("ZipArray_S", 0.0, Zip_Rout, voidThk/2,  0, 2*pi);

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
    new k100_ZipParameterisation(voidThk, //total thickness
		                 NbZipsPerTower,   // NoZips/Tower
				 Zip_Space, // Z spacing between 
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

    zHeightAboveFloor-=positionZipArray.z();            //minus sign because tower rotated
    zHeightAboveFloor-=zipParam->GetCoordinates(0).z(); //minus sign because tower rotated
    G4cout << "Zip Height Above Floor (2nd calc): " << zHeightAboveFloor  << G4endl;
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


    G4ThreeVector Zip0Coord = zipParam->GetCoordinates(0);
    G4cout << "Coordinates of Zip w/ copy no = 0: " << Zip0Coord.x() << "," << Zip0Coord.y() << "," << Zip0Coord.z() << G4endl;
    G4ThreeVector Zip1Coord = zipParam->GetCoordinates(1);
    G4cout << "Coordinates of Zip w/ copy no = 1: " << Zip1Coord.x() << "," << Zip1Coord.y() << "," << Zip1Coord.z() << G4endl;

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

    //if(collID==-1){
    if(true){
      azipSD1 = new k100_ZipSD(detectorZipSDname, towerNb);
      G4cout << "Tower 1 is detector " << towerNb << G4endl;
      k100CollName[detectorZipSDname] = towerNb;
      k100CollPoint[detectorZipSDname] = azipSD1;
      SDman->AddNewDetector(azipSD1);
      //    G4cout << "#### DetCon : zipCollID[ii]  " << SDman->GetCollectionID(detectorZipSDname) << G4endl;
    }
    else{
      azipSD1 = k100CollPoint[detectorZipSDname];
    }
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
      //SDman->AddNewDetector(azipSD2);
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
      //SDman->AddNewDetector(azipSD3);
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
      //SDman->AddNewDetector(azipSD4);
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
      //SDman->AddNewDetector(azipSD5);
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
  G4cout << "z shift for Tower Guts: " << Tower_zPcut[1] - zPcut[1] << G4endl;
  G4cout << "Height of copper upper tower: " << zPcut[1] - zPcut[0] << G4endl;
  G4Polyhedra* cu_tops=new G4Polyhedra("cu_top",0.*deg,360.*deg,6,nZcut,zPcut,rIcut,rOcut);
  G4LogicalVolume* cu_topl1 = new G4LogicalVolume(cu_tops,towerMat,"cutl1");
  G4PVPlacement* cu_topp1 = new G4PVPlacement(0,position_top,"cutp1",cu_topl1,physicalTower,false,0);
  cu_topl1->SetVisAttributes(VisAttCu1);  
  
  //----------------------------------------------------------------------
  //Next, the lower cap of the `copper upper tower', mass 0.233 kg/tower
  //----------------------------------------------------------------------
  G4ThreeVector position_clc = G4ThreeVector(0,0,+Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - zPclc[1]);
  G4cout << "z shift for Tower lower cap: " << Tower_zPcut[1] - (zPcut[1]-zPcut[0]) << G4endl;
  G4cout << "Height of copper lower cap: " << zPclc[1] - zPclc[0] << G4endl;
  G4Polyhedra* cu_clcs=new G4Polyhedra("cu_clc",0.*deg,360.*deg,6,nZclc,zPclc,rIclc,rOclc);
  G4LogicalVolume* cu_clcl1 = new G4LogicalVolume(cu_clcs,towerMat,"clcl1");
  G4PVPlacement* cu_clcp1 = new G4PVPlacement(0,position_clc,"clcp1",cu_clcl1,physicalTower,false,0);
  cu_clcl1->SetVisAttributes(VisAttCu2); 
     
  //------------------------------------------------------------------
  //Next, connector tube, mass 0.140 kg/tower (next ring excluded...)
  //------------------------------------------------------------------
  G4ThreeVector position_ctu = G4ThreeVector(0,0,+Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - zHctu);
  G4cout << "z shift for spool: " << Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - zHctu << G4endl;
  G4cout << "Height of spool: " << 2*zHctu << G4endl;
  G4Tubs* cu_ctus=new G4Tubs("cu_ctu",rIctu,rOctu,zHctu,0.*deg,360.*deg);
  G4LogicalVolume* cu_ctul1=new G4LogicalVolume(cu_ctus,towerMat,"ctul1");
  G4PVPlacement* cu_ctup1=new G4PVPlacement(0,position_ctu,"ctup1",cu_ctul1,physicalTower,false,0);
  cu_ctul1->SetVisAttributes(VisAttCu3); 


  //--------------------------------------------------------------
  //Next, ring at base of connector tube, mass 0.051 kg/tower
  //--------------------------------------------------------------
  G4ThreeVector position_rb = G4ThreeVector(0,0, Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu); // QUESTION : Is this ring centered on the bottom of the tower ??
  G4cout << "z shift for base ring: " << Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu << G4endl;
  G4cout << "Height of base ring: " << 2*zHrb << G4endl;
  G4Tubs* cu_rbs=new G4Tubs("cu_rb",rIrb,rOrb,zHrb,0.*deg,360.*deg);
  G4LogicalVolume* cu_rbl1=new G4LogicalVolume(cu_rbs,towerMat,"crbl1");
  G4PVPlacement* cu_rbp1=new G4PVPlacement(0,position_rb,"crbp1",cu_rbl1,physicalTower,false,0);
  cu_rbl1->SetVisAttributes(VisAttCu1); 

  //------------------------------------------------------------------
  //Next, the upper cap of the detector housing, mass 0.071 kg/tower
  // (about 0.044 kg/tower actually included in ring at base above)
  //------------------------------------------------------------------
  G4ThreeVector position_udh = G4ThreeVector(0,0, Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu - zHrb - zPudh[1]); 
  G4cout << "z shift for upper det housing cap: " << Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu - zHrb - zPudh[1] << G4endl;
  G4cout << "Height of upper det housing cap: " << zPudh[1] - zPudh[0] << G4endl;
  G4Polyhedra* cu_udhs=new G4Polyhedra("cu_udh",0.*deg,360.*deg,6,nZudh,zPudh,rIudh,rOudh);
  G4LogicalVolume* cu_udhl1 = new G4LogicalVolume(cu_udhs,towerMat,"udhl1");
  G4PVPlacement* cu_udhp1 = new G4PVPlacement(0, position_udh,"udhp1",cu_udhl1,physicalTower,false,0);
  cu_udhl1->SetVisAttributes(VisAttCu2);
  

  //--------------------------------------------------------------
  //Next, the side detector housing, mass 0.384 kg/tower
  //--------------------------------------------------------------
  G4ThreeVector position_sdh = G4ThreeVector(0,0, Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu - zHrb - (zPudh[1]-zPudh[0]) - zPsdh[1]); 
  G4cout << "z shift for side det housing: " << Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu - zHrb - (zPudh[1]-zPudh[0]) - zPsdh[1] << G4endl;
  G4cout << "Height of side det housing: " << zPsdh[1] - zPsdh[0] << G4endl;
  G4Polyhedra* cu_sdhs=new G4Polyhedra("cu_sdh",0.*deg,360.*deg,6,nZsdh,zPsdh,rIsdh,rOsdh);
  G4LogicalVolume* cu_sdhl1 = new G4LogicalVolume(cu_sdhs,towerMat,"sdhl1");
  G4PVPlacement* cu_sdhp1 = new G4PVPlacement(0,position_sdh,"sdhp1",cu_sdhl1,physicalTower,false,0);
  cu_sdhl1->SetVisAttributes(VisAttCu3);
  //cu_sdhl1->SetVisAttributes(G4VisAttributes::Invisible);  // Make Invisible


  //------------------------------------------------------------------
  //Next, the lower cap of the detector housing, mass 0.078 kg/tower
  //------------------------------------------------------------------
  G4ThreeVector position_ldh = G4ThreeVector(0,0, Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu - zHrb - (zPudh[1]-zPudh[0]) - (zPsdh[1]-zPsdh[0]) - zPldh[1] ); 
  G4cout << "z shift for lower det housing cap: " << Tower_zPcut[1] - (zPcut[1]-zPcut[0]) - (zPclc[1]-zPclc[0]) - 2*zHctu - zHrb - (zPudh[1]-zPudh[0]) - (zPsdh[1]-zPsdh[0]) - zPldh[1] << G4endl;
  G4cout << "Height of lower det housing cap: " << zPldh[1] - zPldh[0] << G4endl;
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
  //G4PVPlacement* physcoax1 = new G4PVPlacement(0,G4ThreeVector(0.*cm,-(sdcx_thicknessH+rOsdh[0]),coaxpos),"coax1",logiccoax1,physicalTower,false,0);

  //G4PVPlacement* physcoax1 = new G4PVPlacement(0,G4ThreeVector(0.*cm,-4.18*cm,coaxpos),"coax1",logiccoax1,physicalTower,false,0);

  coaxpos = (Tower_zPcut[0]+zPldh[1] + (6-0.5-1)*Zip_Househeight + sdcx_lenH[1]);
 // G4PVPlacement* physcoax2 = new G4PVPlacement(coaxrotation1,G4ThreeVector(0+3.62*cm,0-2.09*cm,coaxpos),"coax2",logiccoax2,physicalTower,false,0);

  coaxpos = (Tower_zPcut[0]+zPldh[1] + (6-0.5-2)*Zip_Househeight + sdcx_lenH[2]);
  //G4PVPlacement* physcoax3 = new G4PVPlacement(coaxrotation2,G4ThreeVector(0+3.62*cm,0+2.09*cm,coaxpos),"coax3",logiccoax3,physicalTower,false,0);

  coaxpos = (Tower_zPcut[0]+zPldh[1] + (6-0.5-5)*Zip_Househeight + sdcx_lenH[3]);
  //G4PVPlacement* physcoax4 = new G4PVPlacement(0,G4ThreeVector(0,sdcx_thicknessH+rOsdh[0],coaxpos),"coax4",logiccoax4,physicalTower,false,0);  

  //G4PVPlacement* physcoax4 = new G4PVPlacement(0,G4ThreeVector(0,0+4.18*cm,coaxpos),"coax4",logiccoax4,physicalTower,false,0);

  coaxpos = (Tower_zPcut[0]+zPldh[1] + (6-0.5-4)*Zip_Househeight + sdcx_lenH[4]);
  //G4PVPlacement* physcoax5 = new G4PVPlacement(coaxrotation1,G4ThreeVector(0-3.62*cm,0+2.09*cm,coaxpos),"coax5",logiccoax5,physicalTower,false,0);

  coaxpos = (Tower_zPcut[0]+zPldh[1] + (6-0.5-5)*Zip_Househeight + sdcx_lenH[5]);
  //G4PVPlacement* physcoax6 = new G4PVPlacement(coaxrotation2,G4ThreeVector(-(sdcx_thicknessH+rOsdh[0])*sin(60*deg),-(sdcx_thicknessH+rOsdh[0])*cos(60*deg),coaxpos),"coax6",logiccoax6,physicalTower,false,0);

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

  if(FirstDetGe)
    DetMaterials[0]  = 1; //Ge  ///Ge = 1 Si = 0
  else{
    DetMaterials[0]  = 0; //Si  ///Ge = 1 Si = 0
  }
  if(ConstructVetoBool)    {ConstructVeto(logicalWorld);}
  if(ConstructShieldsBool) {ConstructShields(logicalWorld);}
  if(ConstructIceBoxBool)  {ConstructIceBox(logicalWorld);}
  if(ConstructFloorBool)  {ConstructFloor(physicalWorld);}
  if(ConstructWallsBool)  {ConstructWalls(physicalWorld);}
  if(ConstructCeilingBool)  {ConstructCeiling(physicalWorld);}
  if(ConstructWestReflectorBool)  {ConstructWestReflector(physicalWorld);}
  if(ConstructFrameBool)  {ConstructFrame(physicalWorld);}
  if(ConstructPuBeSourceAndShieldBool)  {ConstructPuBeSourceAndShield(physicalWorld);}
  if(ConstructThermalNeutronBoxBool)  {ConstructThermalNeutronBox(physicalWorld);}
  if(ConstructNaIArrayBool) {ConstructNaIArray(logicalWorld);} //spandey

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
  polyVis->SetForceSolid(true);
  
  // base panel setup
  G4ThreeVector panelPosition = G4ThreeVector(frame_x+10*2.54*cm,frame_y+12*2.54*cm,frame_z-23.5*2.54*cm);
  G4Box* baseShield = new G4Box("baseShield", 11*2.54*cm, 11*2.54*cm, 4*2.54*cm);
  G4LogicalVolume* logicBase = new G4LogicalVolume(baseShield, polyMat, "logicBase",0,0,0);
  if(shieldParams.addBasePoly) //usually dont put poly panels on base 
    new G4PVPlacement(0,panelPosition,"physicBase",logicBase,physicalWorld,false,0);
  logicBase->SetVisAttributes(polyVis);

  // square side panels
  // initial setup
  G4Box* squareBase = new G4Box("squareBase",11*2.54*cm,4*2.54*cm,27.5*2.54*cm);
  //G4Tubs* hole = new G4Tubs("hole", 0, .55*2.54*cm, 2.1*2.54*cm, 0, 2*pi);
  G4Tubs* hole = new G4Tubs("hole", 0, .55*2.54*cm, 4.1*2.54*cm, 0, 2*pi);
  G4ThreeVector off(1*m,1*m,1*m);
  G4RotationMatrix holeRot;
  holeRot.rotateX(pi/2*rad);
  G4Transform3D offset(holeRot,off);
  //G4SubtractionSolid* oldSquare = new G4SubtractionSolid("oldSquare", squareBase, hole, offset);
  G4SubtractionSolid* oldSquare = (G4SubtractionSolid*) squareBase; //not good to put a subtraction outside original solid, especially not just for convenience, try casting
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

  //can make this sensitive to count flux for Pu/Be runs
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  if(pubeNaIParams.westPolySensitivity){
    G4String SDname = "WestPoly";
    G4int collID = -1; collID = SDman->GetCollectionID(SDname);
    k100_StdSD* polyFlux;
    ConstructGenericSensitiveInt=2; //?FIXME I actually forgot what role this is supposed to play 

    k100CollName[SDname] = 0;
    polyFlux = new k100_StdSD(SDname,k100CollName[SDname]);
    G4cout << "West poly walll flux counter " << k100CollName[SDname] << G4endl;
    k100CollPointStd[SDname] = polyFlux;
    SDman->AddNewDetector(polyFlux);
    logicSquare->SetSensitiveDetector(polyFlux);

  }

  panelPosition=G4ThreeVector(frame_x+10*2.54*cm,frame_y-5*2.54*cm,frame_z);
  G4cout << "Is Pu/Be Mod? " << shieldParams.mod << G4endl;
  if((shieldParams.mod!=2) || (ConstructPuBeSourceAndShieldBool==false)) //mod 2 is simply without this poly
    new G4PVPlacement(0,panelPosition,"physicSquare",logicSquare,physicalWorld,false,0);
  panelPosition=G4ThreeVector(frame_x+10*2.54*cm,frame_y+29*2.54*cm,frame_z);
  if(((shieldParams.mod!=1)&&(shieldParams.mod!=2)&&(!pubeNaIParams.doOrb)) || (ConstructPuBeSourceAndShieldBool==false)) //mod 1 and 2 is simply without this poly
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
  if(!shieldParams.addNaISouth&!shieldParams.HPGeboron&!ConstructNaIArrayBool) //only place south panels when NOT using NaI detectors or boron HPGe
    new G4PVPlacement(0,panelPosition,"physicRect",logicRect,physicalWorld,false,0);
  panelPosition=G4ThreeVector(frame_x+25*2.54*cm,frame_y+12*2.54*cm,frame_z);
  new G4PVPlacement(0,panelPosition,"physicRect1",logicRect,physicalWorld,false,1);
  // //spandey updating north
  // if(!ConstructNaIArrayBool)
  //   new G4PVPlacement(0,panelPosition,"physicRect1",logicRect,physicalWorld,false,1);
  logicRect->SetVisAttributes(polyVis);

 // -------------- Make NaI detector frame (and place upon request) ------------

 //variables
 G4double plateSpace=3.0*2.54*cm;

 // lead visuals
 G4VisAttributes* detVis = new G4VisAttributes(G4Colour(255/255.,0/255.,255/255.));
 detVis->SetForceSolid(false);
 detVis->SetForceWireframe(true);  //I want a Wireframe
 G4VisAttributes* naiAndLightGuideVis = new G4VisAttributes(G4Colour(255/255.,255/255.,0/255.));
 naiAndLightGuideVis->SetForceSolid(true);

 if(shieldParams.addNaISouth){
 
	 
   //FIXME using simple supports, could use 8020 
   G4Box* naiSupport = new G4Box("naiSupport",0.5*1.5*2.54*cm,0.5*1.5*2.54*cm,30.0*2.54*cm);
   G4LogicalVolume* logicNaISupport = new G4LogicalVolume(naiSupport,G4NISTAl,"logicNaISupport",0,0,0);

   G4ThreeVector supportPosition;
   supportPosition=G4ThreeVector(frame_x-(13.25-0.25-0.5*1.5)*2.54*cm,frame_y+(12-0.5*(18-1.5))*2.54*cm,frame_z+(30-27)*2.54*cm); //pole shifted to sit on platform
   new G4PVPlacement(0,supportPosition,"physicNaISupport0",logicNaISupport,physicalWorld,false,0);
   supportPosition=G4ThreeVector(frame_x-(13.25-0.25-0.5*1.5)*2.54*cm,frame_y+(12+0.5*(18-1.5))*2.54*cm,frame_z+(30-27)*2.54*cm); //pole shifted to sit on platform
   new G4PVPlacement(0,supportPosition,"physicNaISupport0",logicNaISupport,physicalWorld,false,0);
   logicNaISupport->SetVisAttributes(frameVis);

   //plates
   G4VSolid* naiPlate = new G4Box("naiPlate",0.5*0.25*2.54*cm,0.5*18.0*2.54*cm,0.5*18.0*2.54*cm);
   G4Tubs* large_hole = new G4Tubs("large_hole",0.0,0.5*8.0*2.54*cm,2*2.54*cm,0.0,2*pi);
   //have to make a rotation for the corner plate cut 
   G4RotationMatrix *holerot = new G4RotationMatrix;
   holerot->rotateY(-M_PI/2.0*rad);
   naiPlate = new G4SubtractionSolid("naiPlate_Cut",naiPlate,large_hole,holerot,G4ThreeVector(0,0,0));
   G4LogicalVolume* logicNaIPlate = new G4LogicalVolume(naiPlate,G4NISTAl,"logicNaIPlate",0,0,0);

   G4ThreeVector platePosition;
   platePosition=G4ThreeVector(frame_x-(13.25-0.25-1.5-0.5*0.25)*2.54*cm,frame_y+(12)*2.54*cm,frame_z-(0.5*18)*2.54*cm-0.5*plateSpace); 
   new G4PVPlacement(0,platePosition,"physicNaIPlate0",logicNaIPlate,physicalWorld,false,0);
   platePosition=G4ThreeVector(frame_x-(13.25-0.25-1.5-0.5*0.25)*2.54*cm,frame_y+(12)*2.54*cm,frame_z+(0.5*18)*2.54*cm+0.5*plateSpace); 
   new G4PVPlacement(0,platePosition,"physicNaIPlate1",logicNaIPlate,physicalWorld,false,0);
   logicNaIPlate->SetVisAttributes(frameVis);

   //detectors
   G4VSolid* collet = new G4Tubs("collet",0,0.5*10.0*2.54*cm,0.5*1.0*2.54*cm,0.0,2*pi);
   G4VSolid* casing = new G4Tubs("casing",0,0.5*8.0*2.54*cm,0.5*8.0*2.54*cm,0.0,2*pi);
   G4VSolid* naiHouse = new G4UnionSolid("naiHouse",collet,casing,0,G4ThreeVector(0,0,+0.5*(8.0+1.0)*2.54*cm));
   G4LogicalVolume* logicNaIHousing = new G4LogicalVolume(naiHouse,G4NISTAl,"logicNaIHousing",0,0,0);
   
   G4ThreeVector detPosition;
   detPosition=G4ThreeVector(frame_x-(13.25-0.25-1.5+0.5*1.0)*2.54*cm,frame_y+(12)*2.54*cm,frame_z-(0.5*18)*2.54*cm-0.5*plateSpace); 
   G4PVPlacement *det0 = new G4PVPlacement(holerot,detPosition,"physicNaIDet0",logicNaIHousing,physicalWorld,false,0);
   detPosition=G4ThreeVector(frame_x-(13.25-0.25-1.5+0.5*1.0)*2.54*cm,frame_y+(12)*2.54*cm,frame_z+(0.5*18)*2.54*cm+0.5*plateSpace); 
   G4PVPlacement *det1 = new G4PVPlacement(holerot,detPosition,"physicNaIDet1",logicNaIHousing,physicalWorld,false,0);
   logicNaIHousing->SetVisAttributes(detVis);

   //detectors
   G4double naiThk = 6.0; //in inches
   G4double sp = 0.0; //in inches, theoretical space between LG and NaI, nominally this is zero but use it for visualization
   G4double lgThk = (9.0-6.0-0.125-sp); // in inches
   G4VSolid* nai = new G4Tubs("nai",0,0.5*(8.0-0.25)*2.54*cm,0.5*(naiThk)*2.54*cm,0.0,2*pi);
   G4VSolid* lg = new G4Tubs("lg",0,0.5*(8.0-0.25)*2.54*cm,0.5*(lgThk)*2.54*cm,0.0,2*pi);
   G4LogicalVolume* logicNaIB = new G4LogicalVolume(nai,G4NISTNaI,"logicNaIB",0,0,0);
   G4LogicalVolume* logicNaIC = new G4LogicalVolume(nai,G4NISTNaI,"logicNaIC",0,0,0);
   G4LogicalVolume* logicLG = new G4LogicalVolume(lg,G4NISTlucite,"logicLG",0,0,0);

   G4ThreeVector matPosition;
   matPosition=G4ThreeVector(0,0,(0.5+8.0-(naiThk/2.0)-0.125)*2.54*cm); 
   new G4PVPlacement(0,matPosition,"physicNaI0",logicNaIB,det0,false,0);
   matPosition=G4ThreeVector(0,0,(-0.5+(lgThk/2.0))*2.54*cm); 
   new G4PVPlacement(0,matPosition,"physicLG0",logicLG,det0,false,0);
   matPosition=G4ThreeVector(0,0,(0.5+8.0-(naiThk/2.0)-0.125)*2.54*cm); 
   new G4PVPlacement(0,matPosition,"physicNaI1",logicNaIC,det1,false,0);
   matPosition=G4ThreeVector(0,0,(-0.5+(lgThk/2.0))*2.54*cm); 
   new G4PVPlacement(0,matPosition,"physicLG1",logicLG,det1,false,0);
   logicNaIB->SetVisAttributes(naiAndLightGuideVis);
   logicNaIC->SetVisAttributes(naiAndLightGuideVis);
   logicLG->SetVisAttributes(naiAndLightGuideVis);

   //------------------------------------------------ 
   // Coincident Sensitive detectors
   //------------------------------------------------ 
      
   // Prepare to declare sensitive detectors
   //G4SDManager* SDman = G4SDManager::GetSDMpointer();

   if(pubeNaIParams.NaIsensitivity){
     G4String SDname = "NaIB";
     G4int collID = -1; collID = SDman->GetCollectionID(SDname);
     k100_StdSD* naiSD0;
     ConstructGenericSensitiveInt=2; //?FIXME I actually forgot what role this is supposed to play 

     k100CollName[SDname] = 1; //was 7
     naiSD0 = new k100_StdSD(SDname,k100CollName[SDname]);
     G4cout << "NaI B is Detector " << k100CollName[SDname] << G4endl;
     k100CollPointStd[SDname] = naiSD0;
     SDman->AddNewDetector(naiSD0);
     logicNaIB->SetSensitiveDetector(naiSD0);

     SDname = "NaIC";
     collID = -1; collID = SDman->GetCollectionID(SDname);
     k100_StdSD* naiSD1;
     ConstructGenericSensitiveInt=2; //?FIXME I actually forgot what role this is supposed to play 

     k100CollName[SDname] = 2; //was 8
     naiSD1 = new k100_StdSD(SDname,k100CollName[SDname]);
     G4cout << "NaI C is Detector " << k100CollName[SDname] << G4endl;
     k100CollPointStd[SDname] = naiSD1;
     SDman->AddNewDetector(naiSD1);
     logicNaIC->SetSensitiveDetector(naiSD1);
   }

  }// end addNaISouth if statement

  if(shieldParams.HPGeboron){ //this should be mutually exclusive with addNaISouth
 
   G4Tubs* dewar = new G4Tubs("dewar", ((17.0/2.0)-(1/2.0))*2.54*cm
       ,(17.0/2.0)*2.54*cm, (16.0/2.0)*2.54*cm, 0, 2*pi);
   G4LogicalVolume* logicDewar = new G4LogicalVolume(dewar,G4NISTstainless,"logicDewar",0,0,0);

   G4ThreeVector dewarPosition;
   //dewarPosition=G4ThreeVector(frame_x-(13.25-0.25-0.5*1.5)*2.54*cm,frame_y+(12-0.5*(18-1.5))*2.54*cm,frame_z+((16.0/2.0)-27)*2.54*cm); //dewar shifted to sit on platform
   dewarPosition=G4ThreeVector(-53*cm,-5*cm,frame_z+((16.0/2.0)-27)*2.54*cm); //dewar shifted to sit on platform
   new G4PVPlacement(0,dewarPosition,"physicDewar",logicDewar,physicalWorld,false,0);
   // frame visuals
   G4VisAttributes* dewarVis = new G4VisAttributes(G4Colour(255/255.,0/255.,0/255.));
   //hpgeVis->SetForceSolid(true);
   dewarVis->SetForceWireframe(true);
   logicDewar->SetVisAttributes(dewarVis);
	 

   //Ok, let's build this pain-in-the-ass thing
   G4double habove_dewar_inches = (16+7+2+3.25); //keep in inches
   G4double casethk = 1/4.0; //in inches
   G4double rhpge = ((3-(2*casethk))/2)*2.54*cm; //3" minus 1/8" casing Al
   G4double rhpge_bore = (1/7.0)*rhpge; //measured from Mirion/canberra brochure assuming to-scale
   G4double lhpge = (3+(5/8))*2.54*cm; //anyone's guess on the length, this number was in the manual not clear if it's our model
   G4double lhpge_bore = 0.75*lhpge; //measured from Mirion/canberra brochure.
   //start with Ge:
   //zipGeMat = Germanium;
   G4Tubs* hpge_core = new G4Tubs("hpge_core", rhpge_bore ,rhpge, lhpge_bore/2.0, 0, 2*pi);
   G4Tubs* hpge_top = new G4Tubs("hpge_top", 0 ,rhpge, (lhpge-lhpge_bore)/2.0, 0, 2*pi);
   G4VSolid* hpge = new G4UnionSolid("hpge",hpge_core,hpge_top,0,G4ThreeVector(0,0,(lhpge_bore/2.0)+((lhpge-lhpge_bore)/2.0)));
   G4LogicalVolume* logicHPGe = new G4LogicalVolume(hpge,zipGeMat,"logicHPGe",0,0,0);

   //now do the casing
   //G4NISTAl
   G4double rhpge_case = (3.0/2)*2.54*cm; //3" 
   G4double lhpge_case = (8)*2.54*cm; //casing is def 8"
   G4Tubs* hpge_case = new G4Tubs("hpge_case", 0 ,rhpge_case, lhpge_case/2.0, 0, 2*pi);
   G4LogicalVolume* logicHPGe_case = new G4LogicalVolume(hpge_case,G4NISTAl,"logicHPGe_case",0,0,0);

   //now do some vacuum
   //defaultMat = Vacuum;
   G4double rhpge_vac = ((3-(2*casethk))/2)*2.54*cm; //same radius of Ge block 
   G4double lhpge_vac = (8-(casethk))*2.54*cm; //casing is def 8", make the casing 1/8 thick 
   G4Tubs* hpge_vac = new G4Tubs("hpge_vac", 0 ,rhpge_vac, lhpge_vac/2.0, 0, 2*pi);
   G4LogicalVolume* logicHPGe_vac = new G4LogicalVolume(hpge_vac,defaultMat,"logicHPGe_vac",0,0,0);
  

   //get the rotation for the casing
   G4RotationMatrix *crot = new G4RotationMatrix;
   crot->rotateX(90.*deg);
   crot->rotateY(-45.*deg);
   //G4Transform3D r1(crot, G4ThreeVector(0,0,0));

   G4ThreeVector hpgePosition;
   G4double xdel = (3.5*2.54*cm+1.375*2.54*cm+(lhpge_case/2.0))*sin(45*deg);
   G4double ydel = (3.5*2.54*cm+1.375*2.54*cm+(lhpge_case/2.0))*cos(45*deg);
   hpgePosition=G4ThreeVector(-53*cm+xdel,-5*cm+ydel,frame_z+(lhpge_case/2.0)+(-27+habove_dewar_inches)*2.54*cm); //hpge shifted to be on top of dewar 
   G4PVPlacement* casing = new G4PVPlacement(crot,hpgePosition,"physicHPGe_package",logicHPGe_case,physicalWorld,false,0);
   G4VisAttributes* hpgeVis = new G4VisAttributes(G4Colour(255/255.,0/255.,255/255.));
   hpgeVis->SetForceSolid(true);
   //hpgeVis->SetForceWireframe(true);
   logicHPGe->SetVisAttributes(hpgeVis);
   G4VisAttributes* hpgecaseVis = new G4VisAttributes(G4Colour(255/255.,0/255.,0/255.));
   hpgecaseVis->SetForceWireframe(true);
   //hpgecaseVis->SetForceSolid(true);
   logicHPGe_case->SetVisAttributes(hpgecaseVis);
   G4VisAttributes* hpgevacVis = new G4VisAttributes(G4Colour(255/255.,255/255.,153/255.));
   hpgevacVis->SetForceWireframe(true);
   logicHPGe_vac->SetVisAttributes(hpgevacVis);

   //place everything inside the casing
   G4ThreeVector relPosition;
   relPosition=G4ThreeVector(0,0,(0.5*casethk)*2.54*cm); // 
   G4PVPlacement* vac = new G4PVPlacement(0,relPosition,"vac",logicHPGe_vac,casing,false,0);
   relPosition=G4ThreeVector(0,0,(lhpge_vac/2.0)-(lhpge_bore/2)-(lhpge-lhpge_bore)); // shift up by half-length of vac, down by half-length of bore, down by full length of top 
   G4PVPlacement* hpge_in_vac = new G4PVPlacement(0,relPosition,"hpge_in_vac",logicHPGe,vac,false,0);

   //now make the HPGe sensitive
   G4String SDname = "HPGe";
   G4int collID = -1; collID = SDman->GetCollectionID(SDname);
   k100_StdSD* HPGe_sensitive;
   ConstructGenericSensitiveInt=2; //?FIXME I actually forgot what role this is supposed to play 

   k100CollName[SDname] = 3; 
   HPGe_sensitive = new k100_StdSD(SDname,k100CollName[SDname]);
   k100CollPointStd[SDname] = HPGe_sensitive;
   SDman->AddNewDetector(HPGe_sensitive);
   logicHPGe->SetSensitiveDetector(HPGe_sensitive);

   if(shieldParams.HPGeboron_wshield){
     //make the boron shield
     G4double container_thk = 0.5*cm;
     G4double rcontainer = (5.0/2.0)*2.54*cm;
     G4double rhole = (3.0/2.0)*2.54*cm;
     G4double lcontainer = (6.0)*2.54*cm;
     G4double boron_thk = 1.5*cm; 
     G4double lcontainer_hole = (6.0)*2.54*cm - 2*container_thk - boron_thk;
     G4Tubs* bshield  = new G4Tubs("bshield", 0 ,rcontainer, lcontainer/2.0, 0, 2*pi);
     G4Tubs* bshieldHole = new G4Tubs("bshieldHole", 0,rhole,lcontainer-1.0*2.54*cm,0,2*pi);
    
     G4ThreeVector bshieldHole_shift(0,0,(lcontainer/2.0)+(lcontainer_hole/2.0));
     G4Transform3D off(noRot,bshieldHole_shift);
     G4SubtractionSolid *new_bshield = new G4SubtractionSolid("bshield",bshield,bshieldHole,off);
     G4LogicalVolume* logic_bcontainer = new G4LogicalVolume(new_bshield,polyMat,"logic_bcontainer",0,0,0); //FIXME poly is probably way too dense

     G4double rpowder_in = rhole+container_thk;
     G4double rpowder_out = rpowder_in+boron_thk;
     G4Tubs* bpowder_base  = new G4Tubs("bpowder_base", rpowder_in ,rpowder_out, 
         (lcontainer-2*container_thk-boron_thk)/2.0, 0, 2*pi);
     G4Tubs* bpowder_top  = new G4Tubs("bpowder_top", 0 ,rpowder_out, 
         boron_thk/2.0, 0, 2*pi);
     G4VSolid* bpowder = new G4UnionSolid("bpowder",bpowder_base,bpowder_top,0,G4ThreeVector(0,0,(lcontainer-2*container_thk-boron_thk)/2.0+((boron_thk)/2.0)));
     G4LogicalVolume* logic_bpowder = new G4LogicalVolume(bpowder,sodium_borate_anhydrous,"logic_bpowder",0,0,0); 

     G4ThreeVector bcontainerPosition;
     xdel = (3.5*2.54*cm+1.375*2.54*cm+(lhpge_case/2.0)+(lhpge_case/2.0)+(lcontainer/2.0)-(lcontainer_hole))*sin(45*deg);
     ydel = (3.5*2.54*cm+1.375*2.54*cm+(lhpge_case/2.0)+(lhpge_case/2.0)+(lcontainer/2.0)-(lcontainer_hole))*cos(45*deg);
     //xdel = (3.5*2.54*cm+1.375*2.54*cm+(lhpge_case/2.0))*sin(45*deg);
     //ydel = (3.5*2.54*cm+1.375*2.54*cm+(lhpge_case/2.0))*cos(45*deg);
     bcontainerPosition=G4ThreeVector(-53*cm+xdel,-5*cm+ydel,frame_z+(lhpge_case/2.0)+(-27+habove_dewar_inches)*2.54*cm); //lots of shifts 
     G4PVPlacement* bcontainer = new G4PVPlacement(crot,bcontainerPosition,"physic_bcontainer",logic_bcontainer,physicalWorld,false,0);
     G4VisAttributes* bcontainerVis = new G4VisAttributes(G4Colour(255/255.,255/255.,255/255.));
     //bcontainerVis->SetForceSolid(true);
     bcontainerVis->SetForceWireframe(true);
     logic_bcontainer->SetVisAttributes(bcontainerVis);
    
     //place the boron
     G4PVPlacement* bpowder_phys = new G4PVPlacement(0,G4ThreeVector(0,0,((lcontainer-2*container_thk-boron_thk)/2.0)-lcontainer/2.0+container_thk),"physic_bpowder",logic_bpowder,bcontainer,false,0);
     G4VisAttributes* bpowderVis = new G4VisAttributes(G4Colour(255/255.,0/255.,0/255.));
     bpowderVis->SetForceSolid(true);
     //bpowderVis->SetForceWireframe(true);
     logic_bpowder->SetVisAttributes(bpowderVis);

     //now make the boron stuff sensitive
     SDname = "boron";
     collID = -1; collID = SDman->GetCollectionID(SDname);
     k100_StdSD* boron_sensitive;

     k100CollName[SDname] = 4; 
     boron_sensitive = new k100_StdSD(SDname,k100CollName[SDname]);
     k100CollPointStd[SDname] = boron_sensitive;
     SDman->AddNewDetector(boron_sensitive);
     logic_bpowder->SetSensitiveDetector(boron_sensitive);

     SDname = "boron_case";
     collID = -1; collID = SDman->GetCollectionID(SDname);
     k100_StdSD* boron_case_sensitive;

     k100CollName[SDname] = 5; 
     boron_case_sensitive = new k100_StdSD(SDname,k100CollName[SDname]);
     k100CollPointStd[SDname] = boron_case_sensitive;
     SDman->AddNewDetector(boron_case_sensitive);
     logic_bcontainer->SetSensitiveDetector(boron_case_sensitive);

   }//end conditional for HPGeboron_shield casing

  }// end HPGeboron if statement
 // --------------------- Lead Frame Panels --------------------------
  // This section contains the aluminum fram that surrounds the lead shield.

  // lead frame visuals
  G4VisAttributes* lframeVis = new G4VisAttributes(G4Colour(128/255.,128/255.,128/255.));
  lframeVis->SetForceSolid(false);
 
 // --------------------- Lead Panels --------------------------

  // lead visuals
  //G4VisAttributes* leadVis = new G4VisAttributes(G4Colour(102/255.,0/255.,187/255.));
  G4VisAttributes* leadVis = new G4VisAttributes(G4Colour(255/255.,0/255.,0/255.));
  leadVis->SetForceSolid(true);

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
  if(shieldParams.addBaseLead)//usually do NOT add base lead sheets
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
  if(!pubeNaIParams.doOrb) //don't place this lead if doing orb
    new G4PVPlacement(0,leadPosition,"physicLeadSquare1",logicLeadSquare,physicalWorld,false,1);
  logicLeadSquare->SetVisAttributes(leadVis);


  // rect side panels
  // initial setup

  G4Box* rectLeadBaseY = new G4Box("rectLeadBase",0.25*2.54*cm,0.5*45*2.54*cm,27.5*2.54*cm);
  G4Box* rectLeadBaseY_forNaI = new G4Box("rectLeadBase",0.25*2.54*cm,0.5*19.0*2.54*cm,27.5*2.54*cm);
  G4Box* rectLeadBaseX = new G4Box("rectLeadBase",0.5*45*2.54*cm,0.25*2.54*cm,27.5*2.54*cm);
  
  //place squares
  G4LogicalVolume* logicLeadRectY = new G4LogicalVolume(rectLeadBaseY,shieldPbMat,"logicLeadRectY",0,0,0);
  G4LogicalVolume* logicLeadRectY_forNaI = new G4LogicalVolume(rectLeadBaseY_forNaI,shieldPbMat,"logicLeadRectY_forNaI",0,0,0);

  //spandey
  //if(!shieldParams.addNaISouth&!shieldParams.HPGeboron){//change lead on side with NaI detectors
  if(!shieldParams.addNaISouth&!shieldParams.HPGeboron){//change lead on side with NaI detectors
    
    leadPosition=G4ThreeVector(frame_x-(11.25)*2.54*cm,frame_y+12*2.54*cm,frame_z);
    new G4PVPlacement(0,leadPosition,"physicLeadRectY0",logicLeadRectY,physicalWorld,false,0);
  }
  else if(shieldParams.addNaISouth&!shieldParams.HPGeboron&!ConstructNaIArrayBool){
    
    leadPosition=G4ThreeVector(frame_x-(13.25)*2.54*cm,frame_y+(12+22.5-(19.0/2.0))*2.54*cm,frame_z);
    new G4PVPlacement(0,leadPosition,"physicLeadRectY_forNaI0",logicLeadRectY_forNaI,physicalWorld,false,0);
    leadPosition=G4ThreeVector(frame_x-(13.25)*2.54*cm,frame_y+(12-22.5+(19.0/2.0))*2.54*cm,frame_z);
    new G4PVPlacement(0,leadPosition,"physicLeadRectY_forNaI1",logicLeadRectY_forNaI,physicalWorld,false,0);
  }
  leadPosition=G4ThreeVector(frame_x+(33.25)*2.54*cm,frame_y+12*2.54*cm,frame_z);
  new G4PVPlacement(0,leadPosition,"physicLeadRectY1",logicLeadRectY,physicalWorld,false,1);
  logicLeadRectY->SetVisAttributes(leadVis);
  logicLeadRectY_forNaI->SetVisAttributes(leadVis);

  G4LogicalVolume* logicLeadRectX = new G4LogicalVolume(rectLeadBaseX,shieldPbMat,"logicLeadRectX",0,0,0);
  leadPosition=G4ThreeVector(frame_x+(10.0)*2.54*cm,frame_y+(12+25)*2.54*cm,frame_z);
  if(!pubeNaIParams.doOrb) //don't place this lead if doing orb
    new G4PVPlacement(0,leadPosition,"physicLeadRectX0",logicLeadRectX,physicalWorld,false,0);
  leadPosition=G4ThreeVector(frame_x+(10.0)*2.54*cm,frame_y+(12-25)*2.54*cm,frame_z);
  new G4PVPlacement(0,leadPosition,"physicLeadRectX1",logicLeadRectX,physicalWorld,false,1);
  logicLeadRectX->SetVisAttributes(leadVis);
                          
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
  G4double r30 = 0.5*5.7*2.54*cm;
  G4double thk30 = 0.5*0.12*2.54*cm;
  //G4Tubs* mK30Solid = new G4Tubs("mK30Solid", 5*4.438*2.54*cm,.5*4.5*2.54*cm, .5*9.75*2.54*cm, 0, 2*pi);
  G4Tubs* mK30Solid = new G4Tubs("mK30Solid", r30-thk30,r30, .5*9.75*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicmK30 = new G4LogicalVolume(mK30Solid,towerMat,"logicmK30",0,0,0);
  new G4PVPlacement(0,mK30Pos,"physicmK30",logicmK30,physicalWorld,false,0);
  logicmK30->SetVisAttributes(copperVis);

  // upper endcap
  G4ThreeVector mK30topPos(fridge_x,fridge_y,fridge_z+5.*2.54*cm);
  //G4Tubs* mK30top = new G4Tubs("mK30top", 0, .5*4.5*2.54*cm, .5*.25*2.54*cm,0,2*pi);
  G4Tubs* mK30top = new G4Tubs("mK30top", 0, r30, .5*.25*2.54*cm,0,2*pi);
  G4LogicalVolume* logicmK30top = new G4LogicalVolume(mK30top,towerMat,"logicmK30top",0,0,0);
  new G4PVPlacement(0,mK30topPos,"physicmK30top",logicmK30top,physicalWorld,false,0);
  logicmK30top->SetVisAttributes(copperVis);

  // lower endcap
  G4ThreeVector mK30lowPos(fridge_x,fridge_y,fridge_z-5.*2.54*cm);
  //G4Tubs* mK30lowbase = new G4Tubs("mK30lowbase",0,.5*4.5*2.54*cm,.5*.25*2.54*cm,0,2*pi);
  G4Tubs* mK30lowbase = new G4Tubs("mK30lowbase",0,r30,.5*.25*2.54*cm,0,2*pi);
  //G4Polyhedra* towerHole = new G4Polyhedra("towerHole",0.*deg,360.*deg,6,Tower_nZcut,Tower_zPcut,Tower_rIcut,Tower_rOcut);
  G4Tubs* towerHole = new G4Tubs("towerHole",0,Tower_rOcut[0]+0.5*cm,400*cm,0,2*pi); //FIXME temporary slopply but large hole
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
  G4ThreeVector uCuPos(fridge_x,fridge_y,fridge_z+9.06*2.54*cm);
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
  G4ThreeVector ivcPos(fridge_x,fridge_y,fridge_z-1.3035*2.54*cm+1.0*2.54*cm);
  G4Tubs* ivc = new G4Tubs("ivc", .5*8.87*2.54*cm, .5*9.*2.54*cm, .5*13.483*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicivc = new G4LogicalVolume(ivc,steel,"logicivc",0,0,0);
  new G4PVPlacement(0,ivcPos,"physicivc",logicivc,physicalWorld,false,0);
  logicivc->SetVisAttributes(steelVis); 
  
  // upper endcap
  G4ThreeVector ivcTopPos(fridge_x,fridge_y,fridge_z+5.888*2.54*cm+1.0*2.54*cm);
  G4Tubs* ivcTop = new G4Tubs("ivcTop", .5*4.*2.54*cm, .5*9.*2.54*cm, .5*.9*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicivcTop = new G4LogicalVolume(ivcTop,steel,"logicivcTop",0,0,0);
  new G4PVPlacement(0,ivcTopPos,"physicivcTop",logicivcTop,physicalWorld,false,0);
  logicivcTop->SetVisAttributes(steelVis); 

  // upper pipe
  G4ThreeVector uppPos(fridge_x,fridge_y,fridge_z+14.191*2.54*cm+1.0*2.54*cm);
  G4Tubs* upp = new G4Tubs("upp", .5*3.834*2.54*cm, .5*4.*2.54*cm, .5*17.506*2.54*cm,0,2*pi);
  G4LogicalVolume* logicupp = new G4LogicalVolume(upp,steel,"logicupp",0,0,0);
  new G4PVPlacement(0,uppPos,"physicupp",logicupp,physicalWorld,false,0);
  logicupp->SetVisAttributes(steelVis);

  // upper pipe cap
  G4ThreeVector upcPos(fridge_x,fridge_y,fridge_z+23.144*2.54*cm+1.0*2.54*cm);
  G4Tubs* upc = new G4Tubs("upc", 0,.5*4.*2.54*cm,.5*.4*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicupc = new G4LogicalVolume(upc,steel,"logicupc",0,0,0);
  new G4PVPlacement(0,upcPos,"physicupc",logicupc,physicalWorld,false,0);
  logicupc->SetVisAttributes(steelVis);

  // lower endcap
  G4ThreeVector ivcLowPos(fridge_x,fridge_y,fridge_z-8.495*2.54*cm+1.0*2.54*cm);
  G4Tubs* ivcLowBase = new G4Tubs("ivcLowBase", 0, .5*9.*2.54*cm, .5*.9*2.54*cm,0,2*pi);
  G4SubtractionSolid* ivcLow = new G4SubtractionSolid("ivcLow",ivcLowBase,towerHole);
  G4LogicalVolume* logicivcLow = new G4LogicalVolume(ivcLow,steel,"logicivcLow",0,0,0);
  new G4PVPlacement(0,ivcLowPos,"physicivcLow",logicivcLow,physicalWorld,false,0);
  logicivcLow->SetVisAttributes(steelVis); 

  // electronics box
  G4ThreeVector eboxPos(fridge_x,fridge_y,fridge_z-10.785*2.54*cm+1.0*2.54*cm);
  G4Tubs* ebox = new G4Tubs("ebox", .5*7.87*2.54*cm, .5*8.*2.54*cm, .5*3.68*2.54*cm,0,2*pi);
  G4LogicalVolume* logicebox = new G4LogicalVolume(ebox,steel,"logicebox",0,0,0);
  new G4PVPlacement(0,eboxPos,"physicebox",logicebox,physicalWorld,false,0);
  logicebox->SetVisAttributes(steelVis);

  // lower cap for electronics box
  G4ThreeVector ecapPos(fridge_x,fridge_y,fridge_z-12.75*2.54*cm+1.0*2.54*cm);
  G4Tubs* ecap = new G4Tubs("ecap", 0, .5*8.*2.54*cm, .5*.25*2.54*cm,0,2*pi);
  G4LogicalVolume* logicecap=new G4LogicalVolume(ecap,steel,"logicecap",0,0,0);
  new G4PVPlacement(0,ecapPos,"physicebox",logicecap,physicalWorld,false,0);
  logicecap->SetVisAttributes(steelVis);

  // ---------------------------- Liquid Helium --------------------------
  //

  G4Material *thishelium;
  if(fridgeParams.pure3HeBath)
    thishelium = stillHe;
  else
    thishelium = helium;
  // top piece
  G4ThreeVector topHePos(fridge_x,fridge_y,fridge_z+35.9075*2.54*cm+1.0*2.54*cm);
  G4Tubs* topHe = new G4Tubs("topHe", 0, .5*10.*2.54*cm, .5*18.563*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logictopHe = new G4LogicalVolume(topHe,thishelium,"logictopHe",0,0,0);
  new G4PVPlacement(0,topHePos,"physictopHe",logictopHe,physicalWorld,false,0);
  logictopHe->SetVisAttributes(heliumVis);

  // top center piece
  G4ThreeVector tcHePos(fridge_x,fridge_y,fridge_z+24.985*2.54*cm+1.0*2.54*cm);
  G4Tubs* tcHe = new G4Tubs("tcHe", 0, .5*14.*2.54*cm, .5*3.282*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logictcHe = new G4LogicalVolume(tcHe,thishelium,"logictcHe",0,0,0);
  new G4PVPlacement(0,tcHePos,"physictcHe",logictcHe,physicalWorld,false,0);
  logictcHe->SetVisAttributes(heliumVis);

  // center piece
  G4ThreeVector ceHePos(fridge_x,fridge_y,fridge_z+14.841*2.54*cm+1.0*2.54*cm);
  G4Tubs* ceHe = new G4Tubs("ceHe", .5*4.*2.54*cm, .5*14.*2.54*cm, .5*17.006*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logiceHe = new G4LogicalVolume(ceHe,thishelium,"logiceHe",0,0,0);
  new G4PVPlacement(0,ceHePos,"physiceHe",logiceHe,physicalWorld,false,0);
  logiceHe->SetVisAttributes(heliumVis);
  
  // lower center piece
  G4ThreeVector lcHePos(fridge_x,fridge_y,fridge_z+5.8255*2.54*cm+1.0*2.54*cm);
  G4Tubs* lcHe = new G4Tubs("lcHe", .5*9.*2.54*cm, .5*14.*2.54*cm, .5*1.025*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logiclcHe = new G4LogicalVolume(lcHe,thishelium,"logiclcHe",0,0,0);
  new G4PVPlacement(0,lcHePos,"physiclcHe",logiclcHe,physicalWorld,0,0,0);
  logiclcHe->SetVisAttributes(heliumVis);
  
  // upper lower piece
  G4ThreeVector ulHePos(fridge_x,fridge_y,fridge_z-1.816*2.54*cm+1.0*2.54*cm);
  G4Tubs* ulHe = new G4Tubs("ulHe", .5*9.*2.54*cm, .5*10.*2.54*cm, .5*14.258*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logiculHe = new G4LogicalVolume(ulHe,thishelium,"logiculHe",0,0,0);
  new G4PVPlacement(0,ulHePos,"physiculHe",logiculHe,physicalWorld,0,0,0);
  logiculHe->SetVisAttributes(heliumVis);

  // center lower piece
  G4ThreeVector clHePos(fridge_x,fridge_y,fridge_z-10.91*2.54*cm+1.0*2.54*cm);
  G4Tubs* clHe = new G4Tubs("clHe", .5*8.*2.54*cm, .5*10.*2.54*cm, .5*3.93*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logiclHe = new G4LogicalVolume(clHe,thishelium,"logiclHe",0,0,0);
  new G4PVPlacement(0,clHePos,"physiclHe",logiclHe,physicalWorld,0,0,0);
  logiclHe->SetVisAttributes(heliumVis);

  // bottom piece
  G4ThreeVector botHePos(fridge_x,fridge_y,fridge_z-13.*2.54*cm+1.0*2.54*cm);
  G4Tubs* botHe = new G4Tubs("botHe", 0, .5*10.*2.54*cm, .5*.25*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicbotHe = new G4LogicalVolume(botHe,thishelium,"logicbotHe",0,0,0);
  new G4PVPlacement(0,botHePos,"physicbotHe",logicbotHe,physicalWorld,0,0,0);
  logicbotHe->SetVisAttributes(heliumVis);

  // ------------------------- Shell ----------------------
  // Outer Edge
  G4ThreeVector oeShellPos(fridge_x,fridge_y,fridge_z+12.9045*2.54*cm+2.0*2.54*cm);
  G4Tubs* oeShell = new G4Tubs("oeShell", .5*21.5*2.54*cm, .5*22.*2.54*cm, .5*64.067*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicoeShell = new G4LogicalVolume(oeShell,steel,"logicoeShell",0,0,0);
  new G4PVPlacement(0,oeShellPos,"physicoeShell",logicoeShell,physicalWorld,false,0);
  logicoeShell->SetVisAttributes(steelVis);

  // Top Edge
  G4ThreeVector teShellPos(fridge_x,fridge_y,fridge_z+45.063*2.54*cm+2.0*2.54*cm);
  G4Tubs* teShell = new G4Tubs("teShell", .5*10.*2.54*cm, .5*22.*2.54*cm, .5*.25*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicteShell = new G4LogicalVolume(teShell,steel,"logicteShell",0,0,0);
  new G4PVPlacement(0,teShellPos,"physicteShell",logicteShell,physicalWorld,false,0);
  logicteShell->SetVisAttributes(steelVis);

  // Bottom Edge
  G4ThreeVector btShellPos(fridge_x,fridge_y,fridge_z-19.254*2.54*cm+2.0*2.54*cm);
  G4Tubs* btShell = new G4Tubs("btShell", 0, .5*22.*2.54*cm, .5*.25*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicbtShell = new G4LogicalVolume(btShell,steel,"logicbtShell",0,0,0);
  new G4PVPlacement(0,btShellPos,"physicbtShell",logicbtShell,physicalWorld,false,0);
  logicbtShell->SetVisAttributes(steelVis);

  // Upper Inner (core)
  G4ThreeVector ucShellPos(fridge_x,fridge_y,fridge_z+35.9075*2.54*cm+1.5*2.54*cm);
  G4Tubs* ucShell = new G4Tubs("ucShell", .5*10.1*2.54*cm, .5*10.2*2.54*cm, .5*18.063*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicucShell = new G4LogicalVolume(ucShell,steel,"logicucShell",0,0,0);
  new G4PVPlacement(0,ucShellPos,"physicucShell",logicucShell,physicalWorld,false,0);
  logicucShell->SetVisAttributes(steelVis);

  // Upper Interface
  G4ThreeVector uiShellPos(fridge_x,fridge_y,fridge_z+26.751*2.54*cm+1.5*2.54*cm);
  G4Tubs* uiShell = new G4Tubs("uiShell", .5*10.1*2.54*cm, .5*14.5*2.54*cm, .5*.25*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicuiShell = new G4LogicalVolume(uiShell,steel,"logicuiShell",0,0,0);
  //new G4PVPlacement(0,uiShellPos,"physicuiShell",logicuiShell,physicalWorld,false,0);
  logicuiShell->SetVisAttributes(steelVis);

  // Middle Inner (core)
  G4ThreeVector mcShellPos(fridge_x,fridge_y,fridge_z+15.9695*2.54*cm+1.0*2.54*cm);
  G4Tubs* mcShell = new G4Tubs("mcShell", .5*14.1*2.54*cm, .5*14.5*2.54*cm, .5*21.313*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicmcShell = new G4LogicalVolume(mcShell,steel,"logicmcShell",0,0,0);
  new G4PVPlacement(0,mcShellPos,"physicmcShell",logicmcShell,physicalWorld,false,0);
  logicmcShell->SetVisAttributes(steelVis);

  // Lower Interface
  G4ThreeVector liShellPos(fridge_x,fridge_y,fridge_z+5.188*2.54*cm+1.0*2.54*cm);
  G4Tubs* liShell = new G4Tubs("liShell", .5*10.1*2.54*cm, .5*14.5*2.54*cm, .5*.25*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicliShell = new G4LogicalVolume(liShell,steel,"logicliShell",0,0,0);
  new G4PVPlacement(0,liShellPos,"physicliShell",logicliShell,physicalWorld,false,0);
  logicliShell->SetVisAttributes(steelVis);

  // Lower Inner (core)
  G4ThreeVector lcShellPos(fridge_x,fridge_y,fridge_z-4.031*2.54*cm+1.0*2.54*cm);
  G4Tubs* lcShell = new G4Tubs("lcShell", .5*10.1*2.54*cm, .5*10.5*2.54*cm, .5*18.188*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logiclcShell = new G4LogicalVolume(lcShell,steel,"logiclcShell",0,0,0);
  new G4PVPlacement(0,lcShellPos,"physiclcShell",logiclcShell,physicalWorld,false,0);
  logiclcShell->SetVisAttributes(steelVis);

  // Lower Cap
  G4ThreeVector lpShellPos(fridge_x,fridge_y,fridge_z-13.25*2.54*cm+1.0*2.54*cm);
  G4Tubs* lpShell = new G4Tubs("lpShell", 0, .5*10.5*2.54*cm, .5*.25*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logiclpShell = new G4LogicalVolume(lpShell,steel,"logiclpShell",0,0,0);
  new G4PVPlacement(0,lpShellPos,"physiclpShell",logiclpShell,physicalWorld,false,0);
  logiclpShell->SetVisAttributes(steelVis);

  // ------------------------ Insulation --------------------
  // Upper piece
  G4ThreeVector upInPos(fridge_x,fridge_y,fridge_z+35.9075*2.54*cm+1.0*2.54*cm);
  G4Tubs* upIn = new G4Tubs("upIn", .5*10.5*2.54*cm, .5*21.5*2.54*cm, .5*18.063*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicupIn = new G4LogicalVolume(upIn,super,"logicupIn",0,0,0);
  new G4PVPlacement(0,upInPos,"physicupIn",logicupIn,physicalWorld,false,0);
  logicupIn->SetVisAttributes(superVis);

  // Middle piece
  G4ThreeVector miInPos(fridge_x,fridge_y,fridge_z+15.9695*2.54*cm+1.0*2.54*cm);
  G4Tubs* miIn = new G4Tubs("miIn", .5*14.5*2.54*cm, .5*21.5*2.54*cm, .5*21.813*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicmiIn = new G4LogicalVolume(miIn,super,"logicmiIn",0,0,0);
  new G4PVPlacement(0,miInPos,"physicmiIn",logicmiIn,physicalWorld,false,0);
  logicmiIn->SetVisAttributes(superVis);

  // Lower piece
  G4ThreeVector loInPos(fridge_x,fridge_y,fridge_z-4.156*2.54*cm+1.0*2.54*cm);
  G4Tubs* loIn = new G4Tubs("loIn", .5*10.5*2.54*cm, .5*21.5*2.54*cm, .5*18.438*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicloIn = new G4LogicalVolume(loIn,super,"logicloIn",0,0,0);
  new G4PVPlacement(0,loInPos,"physicloIn",logicloIn,physicalWorld,false,0);
  logicloIn->SetVisAttributes(superVis);

  // Bottom piece
  G4ThreeVector btInPos(fridge_x,fridge_y,fridge_z-16.252*2.54*cm+1.5*2.54*cm);
  G4Tubs* btIn = new G4Tubs("btIn", 0, .5*21.5*2.54*cm, .5*4.754*2.54*cm, 0, 2*pi);
  G4LogicalVolume* logicbtIn = new G4LogicalVolume(btIn,super,"logicbtIn",0,0,0);
  new G4PVPlacement(0,btInPos,"physicbtIn",logicbtIn,physicalWorld,false,0);
  logicbtIn->SetVisAttributes(superVis);
} // ends IceBox Construction

void k100_DetectorConstruction::ConstructFloor(G4VPhysicalVolume*  world)
{
  // THIS AREA FOR THE K100 Floor around the fridge in PAN 43 
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
  G4VisAttributes* concreteVis = new G4VisAttributes(G4Colour(127/255.,128/255.,118/255.));
  concreteVis->SetForceSolid(true);
  G4VisAttributes* holeVis = new G4VisAttributes(G4Colour(0/255.,0/255.,0/255.));
  holeVis->SetForceSolid(true);
  G4VisAttributes* holeCapVis = new G4VisAttributes(G4Colour(255/255.,0/255.,0/255.));
  holeCapVis->SetForceSolid(true);
  G4VisAttributes* polyVis = new G4VisAttributes(G4Colour(255/255.,255/255.,255/255.));
  polyVis->SetForceSolid(true);
  G4VisAttributes* leadVis = new G4VisAttributes(G4Colour(255/255.,0/255.,0/255.));
  leadVis->SetForceSolid(true);


  // ------------------------- Solid Floor ----------------------
  // Outer Edge
  //G4ThreeVector oeShellPos(fridge_x,fridge_y,fridge_z+12.9045*2.54*cm);
  //G4Tubs* oeShell = new G4Tubs("oeShell", .5*21.5*2.54*cm, .5*22.*2.54*cm, .5*64.067*2.54*cm, 0, 2*pi);
  //G4LogicalVolume* logicoeShell = new G4LogicalVolume(oeShell,steel,"logicoeShell",0,0,0);
  //new G4PVPlacement(0,oeShellPos,"physicoeShell",logicoeShell,physicalWorld,false,0);
  //logicoeShell->SetVisAttributes(steelVis);

  //Floor block
  //G4double fridgeHalfHeightToBottomPlate = (12.9045+13.25+0.25)*2.54*cm;
  G4double fridgeHalfHeightToBottomPlate = (12.9045+19.254+0.25)*2.54*cm; //modified 1/1/18 to get floor height right
  //G4double distanceToFloorZ = fridge_z+12.9045*2.54*cm - fridgeHalfHeightToBottomPlate - 21.0*2.54*cm;
  //G4double distanceCenterToFloor = fridgeHalfHeightToBottomPlate + 21.0*2.54*cm;
  G4double distanceCenterToFloor = fridgeHalfHeightToBottomPlate + 21.0*2.54*cm -70.86*mm; //compensate for 70.86mm discrepancy in floor distance 1/1/18
  G4double floorZ = fridge_z+12.9045*2.54*cm - distanceCenterToFloor;
  G4double floorThk = floorZ - (-0.5*world_z);  //need minus sign because it is COORDINATE of bottom of world
  //FIXME world_z should be checked as being big enough to fit a reasonable floor!
  G4double floorMiddleZ = floorZ - 0.5*floorThk;
  G4cout << "fridgeHalfHeightToBottomPlate: " << fridgeHalfHeightToBottomPlate << G4endl;
  G4cout << "distanceCenterToFloor: " << distanceCenterToFloor << G4endl;
  G4cout << "world_z: " << world_z << G4endl;
  G4cout << "floorZ: " << floorZ << G4endl;
  G4cout << "floorThk: " << floorThk << G4endl;
  G4cout << "floorMiddleZ: " << floorMiddleZ << G4endl;
  G4ThreeVector floorBlockPos(0,0,floorMiddleZ);
  G4Box* floorBlockBox = new G4Box("floorBlock",0.5*world_x,0.5*world_y,0.5*floorThk);
  G4LogicalVolume* logicalFloorBlock = new G4LogicalVolume(floorBlockBox,G4NISTconcrete,"logicalFloorBlock",0,0,0);
  G4PVPlacement *floor = new G4PVPlacement(0,floorBlockPos,"physicalFloorBlock",logicalFloorBlock,world,false,0);
  logicalFloorBlock->SetVisAttributes(concreteVis);

  //hole 6ft deep 42" square
  G4double lipThk = 1.5*2.54*cm;
  G4double holeDepth = 6*12.0*2.54*cm; 
  G4double holeWidth = 42.0*2.54*cm;
  G4ThreeVector holePos(0,0,0.5*(floorThk-holeDepth-lipThk));
  G4cout << "holeDepth: " << holeDepth << G4endl;
  G4cout << "holePosZ: " << holePos.z() << G4endl;
  //FIXME world_z should be checked as being big enough!
  G4Box* holeBox = new G4Box("hole",0.5*holeWidth,0.5*holeWidth,0.5*(holeDepth-lipThk)); //minus the lip
  G4LogicalVolume* logicalHole = new G4LogicalVolume(holeBox,G4NISTair,"logicalHoleBox",0,0,0);
  G4PVPlacement *pitHole = new G4PVPlacement(0,holePos,"physicalHoleBox",logicalHole,floor,false,0);
  logicalHole->SetVisAttributes(holeVis);

  //hole cap
  G4ThreeVector holeCapPos(0,0,0.5*(floorThk-lipThk));
  G4Box* holeCap = new G4Box("holeCap",0.5*(holeWidth+2*lipThk),0.5*(holeWidth+2*lipThk),0.5*lipThk); //minus the lip
  G4LogicalVolume* logicalHoleCap = new G4LogicalVolume(holeCap,G4NISTPVC,"logicalHoleCap",0,0,0);
  new G4PVPlacement(0,holeCapPos,"physicalHoleCap",logicalHoleCap,floor,false,0);
  logicalHoleCap->SetVisAttributes(holeCapVis);

  //PuBe source modifications mod=1
  if((pubeNaIParams.mod==1) && (ConstructPuBeSourceAndShieldBool==true)){

    //poly
    G4double distanceFromOrigin = 1.6*m; 
    G4double x = distanceFromOrigin - lipThk + floorZ; //floorz is negative, this is intermediate var
    G4double polyReflectorThk = 4.0*2.54*cm;
    G4Box* polyReflectorBox = new G4Box("polyReflector",0.5*holeWidth,0.5*holeWidth,0.5*polyReflectorThk); 
    G4LogicalVolume* logicalPolyReflector = new G4LogicalVolume(polyReflectorBox,polyMat,"logicalPolyReflector",0,0,0);
    G4ThreeVector reflectorPos(0,0,0.5*(holeDepth-lipThk)-x);
    new G4PVPlacement(0,reflectorPos,"physicalPolyReflector",logicalPolyReflector,pitHole,false,0);
    logicalPolyReflector->SetVisAttributes(polyVis);

    //lead
    G4double leadThk = 4.0*2.54*cm;
    G4double leadWidth = 12.0*2.54*cm;
    G4Box* leadBox0 = new G4Box("lead0",0.5*holeWidth,0.5*leadWidth,0.5*leadThk); 
    G4LogicalVolume* logicalLead0 = new G4LogicalVolume(leadBox0,shieldPbMat,"logicalLead",0,0,0);
    G4ThreeVector leadPos(0,0,floorZ+0.5*leadThk);
    new G4PVPlacement(0,leadPos,"physicalLead0",logicalLead0,world,false,0);
    logicalLead0->SetVisAttributes(leadVis);
    
    leadThk = 8.0*2.54*cm;
    G4Box* leadBox1 = new G4Box("lead0",0.5*holeWidth,0.5*leadWidth,0.5*leadThk); 
    G4LogicalVolume* logicalLead1 = new G4LogicalVolume(leadBox1,shieldPbMat,"logicalLead",0,0,0);
    leadPos = G4ThreeVector(0,0,0.5*(holeDepth-lipThk)-0.5*leadThk);
    new G4PVPlacement(0,leadPos,"physicalLead1",logicalLead1,pitHole,false,0);
    logicalLead1->SetVisAttributes(leadVis);
   

  }//end Source mod 1

} // ends Floor Construction
void k100_DetectorConstruction::ConstructWalls(G4VPhysicalVolume*  world)
{
  // THIS AREA FOR THE K100 South and West walls around the fridge in PAN 43 
  // just to test...some zero crosshairs:
  
  // --------------------------- visuals ------------------------------
  G4VisAttributes* wallVis = new G4VisAttributes(G4Colour(215/255.,127/255.,0/255.));
  wallVis->SetForceSolid(true);
  G4VisAttributes* airVis = new G4VisAttributes(G4Colour(0/255.,191/255.,255/255.));
  airVis->SetForceSolid(true);


  // ------------------------- Wall Placement ----------------------
  //Floor block
  G4double fridgeHalfHeightToBottomPlate = (12.9045+19.254+0.25)*2.54*cm; //modified 1/1/18 to get floor height right
  G4double distanceCenterToFloor = fridgeHalfHeightToBottomPlate + 21.0*2.54*cm -70.86*mm; //compensate for 70.86mm discrepancy in floor distance 1/1/18
  G4double floorZ = fridge_z+12.9045*2.54*cm - distanceCenterToFloor;

  //ceiling height
  G4double fridgeFrameHeight = (2.0+0.25+(7.0/8.0)+8+(7.0/8.0)+4+0.75)*2.54*cm+1.735*m; //FIXME taken from hard-coded other places in code
  G4double wallHeight = fridgeFrameHeight + 2.5*m;
  G4double ceilingHeight = wallHeight + 0.5*m; //distance into hollow part of ceiling
  G4double topSlabThk = 23.0*cm; //thickness of topmost concrete slab above concrete voids
  G4double totalCeilingThk = 73.0*cm; //thickness of topmost concrete slab above concrete voids


  //wall thickness
  G4double dwallthk = (3.0/4.0)*2.54*cm; //3/4" assumption
  G4double wallthk = 5.5*2.54*cm; //walls measured as 5.5" thick approx 
  
 
  //West wall
  //G4Box* westWallBox = new G4Box("westWall",0.5*world_x,0.5*wallthk,-floorZ);
  G4Box* westWallBox = new G4Box("westWall",0.5*world_x,0.5*wallthk,0.5*wallHeight);
  G4LogicalVolume* logicalWestWall = new G4LogicalVolume(westWallBox,G4NISTair,"logicalWestWall",0,0,0);
  G4ThreeVector westWallPos(0,0.5*world_y-0.5*wallthk,0.5*wallHeight-(-floorZ)); //floorZ is negative
  G4PVPlacement *westWall = new G4PVPlacement(0,westWallPos,"physicalWestWall",logicalWestWall,world,false,0);
  logicalWestWall->SetVisAttributes(airVis);

  //G4Box* westDryWallBox = new G4Box("westDryWall",0.5*world_x,0.5*dwallthk,-floorZ);
  G4Box* westDryWallBox = new G4Box("westDryWall",0.5*world_x,0.5*dwallthk,0.5*wallHeight);
  G4LogicalVolume* logicalWestDryWall = new G4LogicalVolume(westDryWallBox,G4NISTGypsum,"logicalWestDryWall",0,0,0);
  G4ThreeVector westDryWallPos(0,0.5*wallthk-0.5*dwallthk,0); //position inside of wall
  G4PVPlacement *westDryWall0 = new G4PVPlacement(0,westDryWallPos,"physicalWestDryWall0",logicalWestDryWall,westWall,false,0);
  westDryWallPos = G4ThreeVector(0,-0.5*wallthk+0.5*dwallthk,0); //position inside of wall
  G4PVPlacement *westDryWall1 = new G4PVPlacement(0,westDryWallPos,"physicalWestDryWall1",logicalWestDryWall,westWall,false,0);
  logicalWestDryWall->SetVisAttributes(wallVis);

  //South wall
  //G4Box* southWallBox = new G4Box("southWall",0.5*wallthk,0.5*(world_y-wallthk),-floorZ); //subtract thickness for other wall
  G4Box* southWallBox = new G4Box("southWall",0.5*wallthk,0.5*(world_y-wallthk),0.5*wallHeight); //subtract thickness for other wall
  G4LogicalVolume* logicalSouthWall = new G4LogicalVolume(southWallBox,G4NISTair,"logicalSouthWall",0,0,0);
  G4ThreeVector southWallPos(-0.5*world_x+0.5*wallthk,-0.5*wallthk,0.5*wallHeight-(-floorZ)); //floorZ is negative
  G4PVPlacement *southWall = new G4PVPlacement(0,southWallPos,"physicalSouthWall",logicalSouthWall,world,false,0);
  logicalSouthWall->SetVisAttributes(airVis);

  //G4Box* southDryWallBox = new G4Box("southDryWall",0.5*dwallthk,0.5*(world_y-wallthk),-floorZ);
  G4Box* southDryWallBox = new G4Box("southDryWall",0.5*dwallthk,0.5*(world_y-wallthk),0.5*wallHeight);
  G4LogicalVolume* logicalSouthDryWall = new G4LogicalVolume(southDryWallBox,G4NISTGypsum,"logicalSouthDryWall",0,0,0);
  G4ThreeVector southDryWallPos(0.5*wallthk-0.5*dwallthk,0,0); //position inside of wall
  G4PVPlacement *southDryWall0 = new G4PVPlacement(0,southDryWallPos,"physicalSouthDryWall0",logicalSouthDryWall,southWall,false,0);
  southDryWallPos = G4ThreeVector(-0.5*wallthk+0.5*dwallthk,0,0); //position inside of wall
  G4PVPlacement *southDryWall1 = new G4PVPlacement(0,southDryWallPos,"physicalSouthDryWall1",logicalSouthDryWall,southWall,false,0);
  logicalSouthDryWall->SetVisAttributes(wallVis);

} // ends Walls Construction
void k100_DetectorConstruction::ConstructCeiling(G4VPhysicalVolume*  world)
{
  // THIS AREA FOR THE K100 South and West walls around the fridge in PAN 43 
  // just to test...some zero crosshairs:
  
  // --------------------------- visuals ------------------------------
  G4VisAttributes* ceilVis = new G4VisAttributes(G4Colour(127/255.,128/255.,118/255.));
  ceilVis->SetForceSolid(true);
  G4VisAttributes* airVis = new G4VisAttributes(G4Colour(0/255.,191/255.,255/255.));
  airVis->SetForceSolid(true);

  // ------------------------- Ceiling Placement ----------------------
  //Floor block
  G4double fridgeHalfHeightToBottomPlate = (12.9045+19.254+0.25)*2.54*cm; //modified 1/1/18 to get floor height right
  G4double distanceCenterToFloor = fridgeHalfHeightToBottomPlate + 21.0*2.54*cm -70.86*mm; //compensate for 70.86mm discrepancy in floor distance 1/1/18
  G4double floorZ = fridge_z+12.9045*2.54*cm - distanceCenterToFloor;

  //ceiling height
  G4double fridgeFrameHeight = (2.0+0.25+(7.0/8.0)+8+(7.0/8.0)+4+0.75)*2.54*cm+1.735*m; //FIXME taken from hard-coded other places in code
  G4double wallHeight = fridgeFrameHeight + 2.5*m;
  G4double ceilingHeight = wallHeight + 0.5*m; //distance into hollow part of ceiling
  G4double topSlabThk = 23.0*cm; //thickness of topmost concrete slab above concrete voids
  G4double totalCeilingThk = 73.0*cm; //thickness of topmost concrete slab above concrete voids

  //slab
  G4Box* slabBox = new G4Box("slab",0.5*world_x,0.5*world_y,0.5*totalCeilingThk);
  G4LogicalVolume* logicalSlab = new G4LogicalVolume(slabBox,G4NISTconcrete,"logicalSlab",0,0,0);
  G4ThreeVector slabPos(0,0,wallHeight-(-floorZ)+0.5*totalCeilingThk);//floorZ is negative
  G4PVPlacement *slab = new G4PVPlacement(0,slabPos,"physicalSlab",logicalSlab,world,false,0);
  logicalSlab->SetVisAttributes(ceilVis);

  //west voids
  G4Box* westVoidBox = new G4Box("westVoid",0.5*0.5*(world_x-30*cm-5.5*2.54*cm),0.5*2.5*m,0.5*(totalCeilingThk-topSlabThk));
  G4LogicalVolume* logicalWestVoid = new G4LogicalVolume(westVoidBox,G4NISTair,"logicalWestVoid",0,0,0);
  G4ThreeVector voidPos(-0.5*world_x+0.5*(0.5*(world_x-30*cm-5.5*2.54*cm))+5.5*2.54*cm,0.5*world_y-0.5*2.5*m-5.5*2.54*cm,-0.5*totalCeilingThk+0.5*(totalCeilingThk-topSlabThk));//floorZ is negative
  new G4PVPlacement(0,voidPos,"physicalWestVoid0",logicalWestVoid,slab,false,0);
  voidPos = G4ThreeVector(0.5*world_x-0.5*(0.5*(world_x-30*cm-5.5*2.54*cm)),0.5*world_y-0.5*2.5*m-5.5*2.54*cm,-0.5*totalCeilingThk+0.5*(totalCeilingThk-topSlabThk));//floorZ is negative
  new G4PVPlacement(0,voidPos,"physicalWestVoid1",logicalWestVoid,slab,false,0);
  logicalWestVoid->SetVisAttributes(airVis);

  //east voids
  G4Box* eastVoidBox = new G4Box("eastVoid",0.5*0.5*(world_x-30*cm-5.5*2.54*cm),0.5*(world_y-2.5*m-30*cm-5.5*2.54*cm),0.5*(totalCeilingThk-topSlabThk));
  G4LogicalVolume* logicalEastVoid = new G4LogicalVolume(eastVoidBox,G4NISTair,"logicalEastVoid",0,0,0);
  voidPos = G4ThreeVector(-0.5*world_x+0.5*(0.5*(world_x-30*cm-5.5*2.54*cm))+5.5*2.54*cm,-0.5*world_y+0.5*(world_y-2.5*m-30*cm-5.5*2.54*cm),-0.5*totalCeilingThk+0.5*(totalCeilingThk-topSlabThk));//floorZ is negative
  new G4PVPlacement(0,voidPos,"physicalEastVoid",logicalEastVoid,slab,false,0);
  voidPos = G4ThreeVector(0.5*world_x-0.5*(0.5*(world_x-30*cm-5.5*2.54*cm)),-0.5*world_y+0.5*(world_y-2.5*m-30*cm-5.5*2.54*cm),-0.5*totalCeilingThk+0.5*(totalCeilingThk-topSlabThk));//floorZ is negative
  new G4PVPlacement(0,voidPos,"physicalEastVoid",logicalEastVoid,slab,false,0);
  logicalEastVoid->SetVisAttributes(airVis);


} // ends Ceiling Construction

void k100_DetectorConstruction::ConstructWestReflector(G4VPhysicalVolume*  world)
{
  // THIS AREA FOR THE K100 West wall reflector 
  // just to test...some zero crosshairs:
  
  // --------------------------- visuals ------------------------------
  G4VisAttributes* polyVis = new G4VisAttributes(G4Colour(255/255.,255/255.,255/255.));
  polyVis->SetForceSolid(true);

  // ------------------------- Panel Placement ----------------------
  //Floor block
  G4double fridgeHalfHeightToBottomPlate = (12.9045+19.254+0.25)*2.54*cm; //modified 1/1/18 to get floor height right
  G4double distanceCenterToFloor = fridgeHalfHeightToBottomPlate + 21.0*2.54*cm -70.86*mm; //compensate for 70.86mm discrepancy in floor distance 1/1/18
  G4double floorZ = fridge_z+12.9045*2.54*cm - distanceCenterToFloor;

  //FIXME hard-coded post locations
  std::vector<double> postPointsX,postPointsY; //X is north/south south toward positive Y is east/west west toward positive
  //anti-clockwise facing south
  postPointsY.push_back((-19.75-16.0)*2.54*cm); 
  postPointsX.push_back(-22.0*2.54*cm); 
  postPointsY.push_back((-19.75-16.0+73.6)*2.54*cm); 
  postPointsX.push_back(-22.0*2.54*cm);
  postPointsY.push_back((-19.75-16.0+73.6)*2.54*cm); 
  postPointsX.push_back((-22.0+44.0)*2.54*cm);
  postPointsY.push_back((-19.75-16.0)*2.54*cm); 
  postPointsX.push_back((-22.0+44.0)*2.54*cm);

  //ceiling height
  G4double fridgeFrameHeight = (2.0+0.25+(7.0/8.0)+8+(7.0/8.0)+4+0.75)*2.54*cm+1.735*m; //FIXME taken from hard-coded other places in code
  G4double wallHeight = fridgeFrameHeight + 2.5*m;
  G4double ceilingHeight = wallHeight + 0.5*m; //distance into hollow part of ceiling
  G4double topSlabThk = 23.0*cm; //thickness of topmost concrete slab above concrete voids
  G4double totalCeilingThk = 73.0*cm; //thickness of topmost concrete slab above concrete voids

  //useful locations
  //G4ThreeVector barrelPos = G4ThreeVector(0,postPointsY[1]-0.5*16*2.54*cm+0.5*27.5*cm-2.25*2.54*cm+0.5*12*2.54*cm+0.5*16.25*2.54*cm,floorZ+0.5*21.0*2.54*cm);
  //G4Tubs* steelBarrel = new G4Tubs("steelBarrel",0,0.5*(16.25)*2.54*cm,0.5*(21.0)*2.54*cm,0,2*pi);
  G4double westmostEdge = postPointsY[1]-0.5*16*2.54*cm+0.5*27.5*cm-2.25*2.54*cm+0.5*12*2.54*cm+0.5*16.25*2.54*cm + 0.5*(16.25)*2.54*cm; //barrel center plus radius

  //setable parameters
  G4double slabReflectorThk = 4.0*2.54*cm; //4 inch thick
  G4double slabReflectorH = 48.0*2.54*cm; //4ft tall 
  G4double slabReflectorW = 36.0*2.54*cm; //3ft wide 

  //slab Reflector
  G4Box* slabReflectorBox = new G4Box("slabReflector",0.5*slabReflectorW,0.5*slabReflectorThk,0.5*slabReflectorH);
  G4LogicalVolume* logicalSlabReflector = new G4LogicalVolume(slabReflectorBox,polyMat,"logicalSlabReflector",0,0,0);
  G4ThreeVector slabPos(0,westmostEdge+0.5*slabReflectorThk,floorZ+0.5*slabReflectorH);//floorZ is negative
  G4PVPlacement *slabReflector = new G4PVPlacement(0,slabPos,"physicalSlabReflector",logicalSlabReflector,world,false,0);
  logicalSlabReflector->SetVisAttributes(polyVis);


} // ends Reflector Construction

void k100_DetectorConstruction::ConstructFrame(G4VPhysicalVolume*  world)
{
  // THIS AREA FOR THE K100 Floor around the fridge in PAN 43 
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
  G4VisAttributes* concreteVis = new G4VisAttributes(G4Colour(127/255.,128/255.,118/255.));
  concreteVis->SetForceSolid(true);
  G4VisAttributes* alVis = new G4VisAttributes(G4Colour(255/255.,0/255.,0/255.));
  alVis->SetForceSolid(true);
  G4VisAttributes* sandVis = new G4VisAttributes(G4Colour(255/255.,255/255.,0/255.));
  sandVis->SetForceSolid(true);


  // ------------------------- Floor Feet and Posts ----------------------
  //G4double fridgeHalfHeightToBottomPlate = (12.9045+13.25+0.25)*2.54*cm;
  G4double fridgeHalfHeightToBottomPlate = (12.9045+19.254+0.25)*2.54*cm; //modified 1/1/18 to get floor height right
  //G4double distanceCenterToFloor = fridgeHalfHeightToBottomPlate + 21.0*2.54*cm;
  G4double distanceCenterToFloor = fridgeHalfHeightToBottomPlate + 21.0*2.54*cm -70.86*mm; //compensate for 70.86mm discrepancy in floor distance 1/1/18
  G4double floorZ = fridge_z+12.9045*2.54*cm - distanceCenterToFloor;
  std::vector<double> postPointsX,postPointsY; //X is north/south south toward positive Y is east/west west toward positive
  //anti-clockwise facing south
  postPointsY.push_back((-19.75-16.0)*2.54*cm); 
  postPointsX.push_back(-22.0*2.54*cm); 
  postPointsY.push_back((-19.75-16.0+73.6)*2.54*cm); 
  postPointsX.push_back(-22.0*2.54*cm);
  postPointsY.push_back((-19.75-16.0+73.6)*2.54*cm); 
  postPointsX.push_back((-22.0+44.0)*2.54*cm);
  postPointsY.push_back((-19.75-16.0)*2.54*cm); 
  postPointsX.push_back((-22.0+44.0)*2.54*cm);

  //Floor feet 
  //FIXME check world_x and world_y are large enough
  G4Box* floorFootC = new G4Box("floorFootC",0.5*16*2.54*cm,0.5*16*2.54*cm,0.5*2.0*2.54*cm);
  G4LogicalVolume* logicalFootC = new G4LogicalVolume(floorFootC,G4NISTconcrete,"logicalFloorFootC",0,0,0);
  for(int i=0;i<postPointsY.size();i++){
    G4ThreeVector footPos_C(postPointsX[i],postPointsY[i],floorZ+1.0*2.54*cm);
    std::ostringstream name;
    name << "physicalFloorFootC_" << i;
    new G4PVPlacement(0,footPos_C,name.str().c_str(),logicalFootC,world,false,0);
  }
  logicalFootC->SetVisAttributes(concreteVis);

  //Post plates 
  G4Box* postPlate = new G4Box("postPlate",0.5*14*2.54*cm,0.5*14*2.54*cm,0.5*0.25*2.54*cm);
  G4LogicalVolume* logicalPostPlate = new G4LogicalVolume(postPlate,G4NISTAl,"logicalPostPlate",0,0,0);
  for(int i=0;i<postPointsY.size();i++){
    G4ThreeVector postPos(postPointsX[i],postPointsY[i],floorZ+(2.0+0.125)*2.54*cm);
    std::ostringstream name;
    name << "physicalPostPlate_" << i;
    new G4PVPlacement(0,postPos,name.str().c_str(),logicalPostPlate,world,false,0);
  }
  logicalPostPlate->SetVisAttributes(alVis);

  //Posts 
  //G4Tubs* post = new G4Tubs("post",0.5*27*cm-0.5*2.54*cm,0.5*27*cm,0.5*1.735*m,0,2*pi); //hollow
  G4Tubs* post = new G4Tubs("post",0,0.5*27*cm,0.5*1.735*m,0,2*pi); //with sand
  G4Tubs* sand = new G4Tubs("post",0,0.5*27*cm-0.5*2.54*cm,0.5*1.735*m,0,2*pi);
  G4LogicalVolume* logicalPost = new G4LogicalVolume(post,G4NISTAl,"logicalPost",0,0,0);
  G4LogicalVolume* logicalSand = new G4LogicalVolume(sand,blastsand,"logicalSand",0,0,0);
  for(int i=0;i<postPointsY.size();i++){
    G4ThreeVector postPos(postPointsX[i],postPointsY[i],floorZ+(2.0+0.25)*2.54*cm+0.5*1.735*m);
    std::ostringstream name,sand_name;
    name << "physicalPost_" << i;
    sand_name << "physicalPostSand_" << i;
    G4PVPlacement *thispost = new G4PVPlacement(0,postPos,name.str().c_str(),logicalPost,world,false,0);
    new G4PVPlacement(0,G4ThreeVector(0,0,0),sand_name.str().c_str(),logicalSand,thispost,false,0);
  }
  logicalPost->SetVisAttributes(alVis);
  logicalSand->SetVisAttributes(sandVis);

  //Post top plates 
  G4Box* postTopPlate = new G4Box("postTopPlate",0.5*12*2.54*cm,0.5*12*2.54*cm,0.5*(7.0/8.0)*2.54*cm);
  G4LogicalVolume* logicalPostTopPlate = new G4LogicalVolume(postTopPlate,G4NISTAl,"logicalPostTopPlate",0,0,0);
  for(int i=0;i<postPointsY.size();i++){
    G4ThreeVector postPos(postPointsX[i],postPointsY[i],floorZ+(2.0+0.25+(7.0/16.0))*2.54*cm+1.735*m);
    std::ostringstream name;
    name << "physicalPostTopPlate_" << i;
    new G4PVPlacement(0,postPos,name.str().c_str(),logicalPostTopPlate,world,false,0);
  }
  logicalPostTopPlate->SetVisAttributes(alVis);

  //Post top assemblies 
  G4Box* postBoxedIBeam = new G4Box("postBoxedIBeam",0.5*6*2.54*cm,0.5*2.75*2.54*cm,0.5*8.0*2.54*cm);
  G4Box* postBoxedIBeam_AirPocket = new G4Box("postBoxedIBeam_AirPocket",0.5*(6-1)*2.54*cm,0.5*0.5*(2.75-0.25)*2.54*cm,0.5*(8.0-0.5)*2.54*cm);
  G4LogicalVolume* logicalPostBoxedIBeam = new G4LogicalVolume(postBoxedIBeam,G4NISTAl,"logicalPostBoxedIBeam",0,0,0);
  G4LogicalVolume* logicalPostBoxedIBeam_AirPocket = new G4LogicalVolume(postBoxedIBeam_AirPocket,G4NISTair,"logicalPostBoxedIBeam_AirPocket",0,0,0);
  for(int i=0;i<postPointsY.size();i++){
    G4ThreeVector postPosWest(postPointsX[i],postPointsY[i]-(6.0-0.5*2.75)*2.54*cm,floorZ+(2.0+0.25+(7.0/8.0)+4.0)*2.54*cm+1.735*m);
    G4ThreeVector postPosEast(postPointsX[i],postPointsY[i]+(6.0-0.5*2.75)*2.54*cm,floorZ+(2.0+0.25+(7.0/8.0)+4.0)*2.54*cm+1.735*m);
    G4ThreeVector airShiftWest(0,(0.5*0.25+0.5*0.5*(2.75-0.25))*2.54*cm,0);
    G4ThreeVector airShiftEast(0,-(0.5*0.25+0.5*0.5*(2.75-0.25))*2.54*cm,0);
    std::ostringstream nameW,nameE,nameWa,nameEa;
    nameW << "physicalPostBoxedIBeamW_" << i;
    nameE << "physicalPostBoxedIBeamE_" << i;
    nameWa << "physicalPostBoxedIBeamW_AirPocket_" << i;
    nameEa << "physicalPostBoxedIBeamE_AirPocket_" << i;
    G4PVPlacement * w = new G4PVPlacement(0,postPosWest,nameW.str().c_str(),logicalPostBoxedIBeam,world,false,0);
    G4PVPlacement * e = new G4PVPlacement(0,postPosEast,nameE.str().c_str(),logicalPostBoxedIBeam,world,false,0);
    new G4PVPlacement(0,airShiftWest,nameWa.str().c_str(),logicalPostBoxedIBeam_AirPocket,w,false,0);
    new G4PVPlacement(0,airShiftEast,nameEa.str().c_str(),logicalPostBoxedIBeam_AirPocket,e,false,0);
  }
  logicalPostBoxedIBeam->SetVisAttributes(alVis);
  logicalPostBoxedIBeam_AirPocket->SetVisAttributes(sandVis);

  //Post top assembly plates 
  G4Box* postTopAPlate = new G4Box("postTopAPlate",0.5*10*2.54*cm,0.5*12*2.54*cm,0.5*(7.0/8.0)*2.54*cm);
  G4LogicalVolume* logicalPostTopAPlate = new G4LogicalVolume(postTopAPlate,G4NISTAl,"logicalPostTopAPlate",0,0,0);
  for(int i=0;i<postPointsY.size();i++){
    G4ThreeVector postPos(postPointsX[i],postPointsY[i],floorZ+(2.0+0.25+(7.0/8.0)+8+(7.0/16.0))*2.54*cm+1.735*m);
    std::ostringstream name;
    name << "physicalPostTopAPlate_" << i;
    new G4PVPlacement(0,postPos,name.str().c_str(),logicalPostTopAPlate,world,false,0);
  }
  logicalPostTopAPlate->SetVisAttributes(alVis);

  //Top Struts
  G4Box* postTopStrutY = new G4Box("postTopStrutY",0.5*4*2.54*cm,0.5*(73.6+8.0)*2.54*cm,0.5*4*2.54*cm);
  G4LogicalVolume* logicalTopStrutY = new G4LogicalVolume(postTopStrutY,G4NISTAl,"logicalPostTopStrutY",0,0,0);
  G4Box* postTopStrutY_AirPocket = new G4Box("postTopStrutY_AirPocket",0.5*(4-0.5)*2.54*cm,0.5*(73.6+8.0)*2.54*cm,0.5*(4-0.5)*2.54*cm);
  G4LogicalVolume* logicalTopStrutY_AirPocket = new G4LogicalVolume(postTopStrutY_AirPocket,G4NISTair,"logicalPostTopStrutY_AirPocket",0,0,0);
  G4Box* postTopStrutX = new G4Box("postTopStrutX",0.5*(44-4)*2.54*cm,0.5*4*2.54*cm,0.5*4*2.54*cm);
  G4LogicalVolume* logicalTopStrutX = new G4LogicalVolume(postTopStrutX,G4NISTAl,"logicalPostTopStrutX",0,0,0);
  G4Box* postTopStrutX_AirPocket = new G4Box("postTopStrutX_AirPocket",0.5*(44-4)*2.54*cm,0.5*(4-0.5)*2.54*cm,0.5*(4-0.5)*2.54*cm);
  G4LogicalVolume* logicalTopStrutX_AirPocket = new G4LogicalVolume(postTopStrutX_AirPocket,G4NISTair,"logicalPostTopStrutX_AirPocket",0,0,0);
    
  G4ThreeVector northYStrutPos(postPointsX[0],postPointsY[0]+(postPointsY[1]-postPointsY[0])/2.0,floorZ+(2.0+0.25+(7.0/8.0)+8+(7.0/8.0)+2)*2.54*cm+1.735*m);
  G4PVPlacement *northYStrut = new G4PVPlacement(0,northYStrutPos,"physicalNorthYStrut",logicalTopStrutY,world,false,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),"physicalNorthYStrut_AirPocket",logicalTopStrutY_AirPocket,northYStrut,false,0);

  G4ThreeVector southYStrutPos(postPointsX[2],postPointsY[0]+(postPointsY[1]-postPointsY[0])/2.0,floorZ+(2.0+0.25+(7.0/8.0)+8+(7.0/8.0)+2)*2.54*cm+1.735*m);
  G4PVPlacement *southYStrut = new G4PVPlacement(0,southYStrutPos,"physicalSouthYStrut",logicalTopStrutY,world,false,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),"physicalSouthYStrut_AirPocket",logicalTopStrutY_AirPocket,southYStrut,false,0);

  G4ThreeVector eastXStrutPos(postPointsX[1]+(postPointsX[2]-postPointsX[1])/2.0,postPointsY[0],floorZ+(2.0+0.25+(7.0/8.0)+8+(7.0/8.0)+2)*2.54*cm+1.735*m);
  G4PVPlacement *eastXStrut = new G4PVPlacement(0,eastXStrutPos,"physicalEastXStrut",logicalTopStrutX,world,false,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),"physicalEastXStrut_AirPocket",logicalTopStrutX_AirPocket,eastXStrut,false,0);

  G4ThreeVector westXStrutPos(postPointsX[1]+(postPointsX[2]-postPointsX[1])/2.0,postPointsY[1],floorZ+(2.0+0.25+(7.0/8.0)+8+(7.0/8.0)+2)*2.54*cm+1.735*m);
  G4PVPlacement *westXStrut = new G4PVPlacement(0,westXStrutPos,"physicalWestXStrut",logicalTopStrutX,world,false,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),"physicalWestXStrut_AirPocket",logicalTopStrutX_AirPocket,westXStrut,false,0);

  logicalTopStrutY->SetVisAttributes(alVis);
  logicalTopStrutY_AirPocket->SetVisAttributes(sandVis);
  logicalTopStrutX->SetVisAttributes(alVis);
  logicalTopStrutX_AirPocket->SetVisAttributes(sandVis);

  //frame plates
  G4VSolid* framePlateMiddle = new G4Box("framePlateMiddle",0.5*(44+4)*2.54*cm,0.5*(33.75)*2.54*cm,0.5*1.5*2.54*cm);
  G4Tubs* large_hole = new G4Tubs("large_hole",0.0,0.5*22.1*2.54*cm,2*2.54*cm,0.0,2*pi);
  framePlateMiddle = new G4SubtractionSolid("framePlateMiddle_Cut",framePlateMiddle,large_hole,NULL,G4ThreeVector(0,0,0));
  // using "light" aluminum to simulate the existence of a lattice of 1/4-20 holes with a pitch of 2" about 99% of area is filled and same fraction
  // of volume becaus they are through-holes
  G4LogicalVolume* logicalFramePlateMiddle = new G4LogicalVolume(framePlateMiddle,lightaluminum,"logicalFramePlateMiddle",0,0,0);

  G4ThreeVector framePlateMiddlePos(0,0,floorZ+(2.0+0.25+(7.0/8.0)+8+(7.0/8.0)+4+0.75)*2.54*cm+1.735*m);
  G4PVPlacement *framePlateMiddlePV = new G4PVPlacement(0,framePlateMiddlePos,"physicalFramePlateMiddle",logicalFramePlateMiddle,world,false,0);
  logicalFramePlateMiddle->SetVisAttributes(alVis);

  G4Box* framePlateEast = new G4Box("framePlateEast",0.5*(44+4)*2.54*cm,0.5*(23.75)*2.54*cm,0.5*0.25*2.54*cm);
  G4LogicalVolume* logicalFramePlateEast = new G4LogicalVolume(framePlateEast,G4NISTAl,"logicalFramePlateEast",0,0,0);

  G4ThreeVector framePlateEastPos(0,-0.5*(33.75+23.75)*2.54*cm,floorZ+(2.0+0.25+(7.0/8.0)+8+(7.0/8.0)+4+0.125)*2.54*cm+1.735*m);
  G4PVPlacement *framePlateEastPV = new G4PVPlacement(0,framePlateEastPos,"physicalFramePlateEast",logicalFramePlateEast,world,false,0);
  logicalFramePlateEast->SetVisAttributes(alVis);

} // ends Frame Construction

void k100_DetectorConstruction::ConstructPuBeSourceAndShield(G4VPhysicalVolume*  world)
{
  // THIS AREA FOR THE K100 Floor around the fridge in PAN 43 
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
  G4VisAttributes* woodVis = new G4VisAttributes(G4Colour(209/255.,179/255.,127/255.));
  woodVis->SetForceSolid(true);
  G4VisAttributes* leadVis = new G4VisAttributes(G4Colour(255/255.,0/255.,0/255.));
  leadVis->SetForceSolid(true);
  G4VisAttributes* alVis = new G4VisAttributes(G4Colour(0/255.,255/255.,0/255.));
  alVis->SetForceSolid(true);
  //alVis->SetForceWireframe(true);
  G4VisAttributes* airVis = new G4VisAttributes(G4Colour(255/255.,255/255.,0/255.));
  airVis->SetForceSolid(true);
  //airVis->SetForceWireframe(true);
  G4VisAttributes* steelVis = new G4VisAttributes(G4Colour(255/255.,0/255.,255/255.));
  steelVis->SetForceSolid(true);
  steelVis->SetForceWireframe(false);
  G4VisAttributes* barrelVis = new G4VisAttributes(G4Colour(255/255.,0/255.,255/255.));
  barrelVis->SetForceSolid(true);
  //barrelVis->SetForceWireframe(false);
  G4VisAttributes* paraffinVis = new G4VisAttributes(G4Colour(0/255.,255/255.,0/255.));
  paraffinVis->SetForceSolid(true);
  //paraffinVis->SetForceWireframe(true);
  G4VisAttributes* luciteVis = new G4VisAttributes(G4Colour(255/255.,165/255.,0/255.));
  luciteVis->SetForceSolid(true);


  // ------------------------- Floor Feet and Posts ----------------------
  //G4double fridgeHalfHeightToBottomPlate = (12.9045+13.25+0.25)*2.54*cm;
  G4double fridgeHalfHeightToBottomPlate = (12.9045+19.254+0.25)*2.54*cm; //modified 1/1/18 to get floor height right
  //G4double distanceCenterToFloor = fridgeHalfHeightToBottomPlate + 21.0*2.54*cm;
  G4double distanceCenterToFloor = fridgeHalfHeightToBottomPlate + 21.0*2.54*cm -70.86*mm; //compensate for 70.86mm discrepancy in floor distance 1/1/18
  G4double floorZ = fridge_z+12.9045*2.54*cm - distanceCenterToFloor;
  std::vector<double> postPointsX,postPointsY; //X is north/south south toward positive Y is east/west west toward positive
  //anti-clockwise facing south
  postPointsY.push_back((-19.75-16.0)*2.54*cm); 
  postPointsX.push_back(-22.0*2.54*cm); 
  postPointsY.push_back((-19.75-16.0+73.6)*2.54*cm); 
  postPointsX.push_back(-22.0*2.54*cm);
  postPointsY.push_back((-19.75-16.0+73.6)*2.54*cm); 
  postPointsX.push_back((-22.0+44.0)*2.54*cm);
  postPointsY.push_back((-19.75-16.0)*2.54*cm); 
  postPointsX.push_back((-22.0+44.0)*2.54*cm);

  if(pubeNaIParams.doR66){
    //Wood
    G4Box* woodBase = new G4Box("woodBase",0.5*(postPointsX[2]-postPointsX[1]) -0.5*16.0*2.54*cm,0.5*10.5*2.54*cm,0.5*1.5*2.54*cm);
    G4LogicalVolume* logicalWoodBase = new G4LogicalVolume(woodBase,wood,"logicalWoodBase",0,0,0);
    G4ThreeVector basePos(0,postPointsY[1]-(0.5*(16.0-10.5)+3.25)*2.54*cm,floorZ+0.5*1.5*2.54*cm);
    new G4PVPlacement(0,basePos,"physicalWoodBase",logicalWoodBase,world,false,0);
    logicalWoodBase->SetVisAttributes(woodVis);

    //lead
    G4Box* leadBase = new G4Box("leadBase",0.5*(postPointsX[2]-postPointsX[1]) -0.5*16.0*2.54*cm - 0.5*2.0*2.54*cm,0.5*8.0*2.54*cm,0.5*0.5*2.54*cm + 0.5*cm);
    G4LogicalVolume* logicalLeadBase = new G4LogicalVolume(leadBase,shieldPbMat,"logicalLeadBase",0,0,0);
    basePos =  G4ThreeVector(0,postPointsY[1]-(0.5*(16.0-10.5)+3.25)*2.54*cm,floorZ+1.5*2.54*cm+0.5*(0.5*2.54*cm + 1*cm));
    new G4PVPlacement(0,basePos,"physicalleadBase",logicalLeadBase,world,false,0);
    logicalLeadBase->SetVisAttributes(leadVis);

    //aluminum box
    G4VSolid* alBox = new G4Box("alBox",0.5*(postPointsX[2]-postPointsX[1]) -0.5*16.0*2.54*cm,0.5*25.5*cm,0.5*25.5*cm);
    G4Box* alCap = new G4Box("alCap",0.5*1.0*cm,0.5*27.5*cm,0.5*27.5*cm);
    alBox = new G4UnionSolid("alBox0",alBox,alCap,0,G4ThreeVector(0.5*(postPointsX[2]-postPointsX[1]) -0.5*16*2.54*cm +0.5*1.0*cm,0,0));
    alBox = new G4UnionSolid("alBox1",alBox,alCap,0,G4ThreeVector(-0.5*(postPointsX[2]-postPointsX[1])+0.5*16*2.54*cm -0.5*1.0*cm,0,0));
    G4LogicalVolume* logicalAlBox = new G4LogicalVolume(alBox,G4NISTAl,"logicalAlBox",0,0,0);
    //G4ThreeVector boxPos(0,postPointsY[1]-0.5*16*2.54*cm+0.5*27.5*cm-2.25*2.54*cm,floorZ+1.5*2.54*cm+0.5*2.54*cm+1*cm+0.5*25.5*cm);
    G4ThreeVector boxPos(0,postPointsY[1]-0.5*16*2.54*cm+0.5*27.5*cm-2.25*2.54*cm,floorZ+2.0*2.54*cm+0.5*27.5*cm);
    G4PVPlacement *alBoxWorld = new G4PVPlacement(0,boxPos,"physicalAlBox",logicalAlBox,world,false,0);
    logicalAlBox->SetVisAttributes(alVis);

    G4Box* airInAlBox = new G4Box("airInAlBox",0.5*28.0*2.54*cm,0.5*(25.5-0.6)*cm,0.5*(25.5-0.6)*cm);
    G4LogicalVolume* logicalAirInAlBox = new G4LogicalVolume(airInAlBox,G4NISTair,"logicalAirInAlBox",0,0,0);
    G4PVPlacement *airBoxWorld = new G4PVPlacement(0,G4ThreeVector(0,0,0),"physicalAirInAlBox",logicalAirInAlBox,alBoxWorld,false,0);
    logicalAirInAlBox->SetVisAttributes(airVis);

    G4Box* alIntPlate = new G4Box("alInt",0.5*(25.5-0.6)*cm + 1.0*2.54*cm,0.5*1.0*2.54*cm,0.5*(25.5-0.6)*cm);
    G4LogicalVolume* logicalAlIntPlate = new G4LogicalVolume(alIntPlate,G4NISTAl,"logicalAlIntPlate",0,0,0);
    new G4PVPlacement(0,G4ThreeVector(0,0,0),"physicalAlIntPlate",logicalAlIntPlate,airBoxWorld,false,0);
    logicalAlIntPlate->SetVisAttributes(alVis);

    //plate lead sits on
    G4Box* alLeadPlate = new G4Box("alLeadPlate",0.5*(24-1.5)*2.54*cm,0.5*(12)*2.54*cm,0.5*(3.0/8.0)*2.54*cm);
    G4LogicalVolume* logicalAlLeadPlate = new G4LogicalVolume(alLeadPlate,G4NISTAl,"logicalAlLeadPlate",0,0,0);
    G4ThreeVector alLeadPlatePos(0,postPointsY[1]-0.5*16*2.54*cm+0.5*27.5*cm-2.25*2.54*cm,floorZ+1.5*2.54*cm+0.5*2.54*cm+1*cm+25.5*cm+0.5*(3.0/8.0)*2.54*cm);
    new G4PVPlacement(0,alLeadPlatePos,"physicalAlLeadPlate",logicalAlLeadPlate,world,false,0);
    logicalAlLeadPlate->SetVisAttributes(alVis);

    //lead stack
    G4Box* leadStack = new G4Box("leadStack",0.5*(24)*2.54*cm,0.5*(12)*2.54*cm,0.5*(16)*2.54*cm);
    G4LogicalVolume* logicalLeadStack = new G4LogicalVolume(leadStack,shieldPbMat,"logicalLeadStack",0,0,0);
    G4ThreeVector leadStackPos(0,postPointsY[1]-0.5*16*2.54*cm+0.5*27.5*cm-2.25*2.54*cm,floorZ+1.5*2.54*cm+0.5*2.54*cm+1*cm+25.5*cm+(3.0/8.0)*2.54*cm+0.5*16*2.54*cm);
    new G4PVPlacement(0,leadStackPos,"physicalLeadStack",logicalLeadStack,world,false,0);
    logicalLeadStack->SetVisAttributes(leadVis);

    G4Box* leadStack1 = new G4Box("leadStack1",0.5*(8)*2.54*cm,0.5*(8)*2.54*cm,0.5*(6)*2.54*cm);
    G4LogicalVolume* logicalLeadStack1 = new G4LogicalVolume(leadStack1,shieldPbMat,"logicalLeadStack1",0,0,0);
    leadStackPos = G4ThreeVector(-(12-4.0)*2.54*cm,postPointsY[1]-0.5*16*2.54*cm+0.5*27.5*cm-2.25*2.54*cm-2.0*2.54*cm,floorZ+1.5*2.54*cm+0.5*2.54*cm+1*cm+25.5*cm+(3.0/8.0)*2.54*cm+16*2.54*cm+0.5*6.0*2.54*cm);
    new G4PVPlacement(0,leadStackPos,"physicalLeadStack1",logicalLeadStack1,world,false,0);
    logicalLeadStack1->SetVisAttributes(airVis);

    G4Box* leadStack2 = new G4Box("leadStack2",0.5*(8)*2.54*cm,0.5*(4)*2.54*cm,0.5*(2)*2.54*cm);
    G4LogicalVolume* logicalLeadStack2 = new G4LogicalVolume(leadStack2,shieldPbMat,"logicalLeadStack2",0,0,0);
    leadStackPos = G4ThreeVector(-(12-4.0)*2.54*cm,postPointsY[1]-0.5*16*2.54*cm+0.5*27.5*cm-2.25*2.54*cm-2.0*2.54*cm,floorZ+1.5*2.54*cm+0.5*2.54*cm+1*cm+25.5*cm+(3.0/8.0)*2.54*cm+16*2.54*cm+6.0*2.54*cm+0.5*2.0*2.54*cm);
    new G4PVPlacement(0,leadStackPos,"physicalLeadStack2",logicalLeadStack2,world,false,0);
    logicalLeadStack2->SetVisAttributes(airVis);

    //barrel
    if(pubeNaIParams.addBarrel){
      G4Tubs* steelBarrel = new G4Tubs("steelBarrel",0,0.5*(16.25)*2.54*cm,0.5*(21.0)*2.54*cm,0,2*pi);
      G4LogicalVolume* logicalSteelBarrel = new G4LogicalVolume(steelBarrel,carbonsteel,"logicalSteelBarrel",0,0,0);
      G4ThreeVector barrelPos = G4ThreeVector(0,postPointsY[1]-0.5*16*2.54*cm+0.5*27.5*cm-2.25*2.54*cm+0.5*12*2.54*cm+0.5*16.25*2.54*cm,floorZ+0.5*21.0*2.54*cm);
      G4PVPlacement *barrelWorld = new G4PVPlacement(0,barrelPos,"physicalBarrel",logicalSteelBarrel,world,false,0);
      logicalSteelBarrel->SetVisAttributes(steelVis);

      G4Tubs* paraffinInsert = new G4Tubs("paraffinInsert",0,0.5*(16.25-0.125)*2.54*cm,0.5*(21.0-(1.0/16.0))*2.54*cm,0,2*pi);
      G4LogicalVolume* logicalParaffinInsert = new G4LogicalVolume(paraffinInsert,G4NISTparaffin,"logicalParaffinInsert",0,0,0);
      G4ThreeVector paraffinPos = G4ThreeVector(0,0,(1.0/32.0)*2.54*cm);
      G4PVPlacement *paraffinWorld = new G4PVPlacement(0,paraffinPos,"physicalParaffin",logicalParaffinInsert,barrelWorld,false,0);
      logicalParaffinInsert->SetVisAttributes(paraffinVis);

      G4cout << "floorZ: " << floorZ << G4endl;
      G4Tubs* pipeInner = new G4Tubs("pipeInner",0,0.5*(2.5)*2.54*cm,0.5*(9.0)*2.54*cm,0,2*pi);
      G4LogicalVolume* logicalPipeInner = new G4LogicalVolume(pipeInner,carbonsteel,"logicalPipeInner",0,0,0);
      G4ThreeVector pipePos = G4ThreeVector(0,0,0.5*(21.0-(1/16.0)-9.0)*2.54*cm-2*2.54*cm);
      G4PVPlacement *pipeWorld = new G4PVPlacement(0,pipePos,"physicalPipe",logicalPipeInner,paraffinWorld,false,0);
      logicalPipeInner->SetVisAttributes(steelVis);

      G4Tubs* pipeAir = new G4Tubs("pipeAir",0,0.5*(2.5-0.25)*2.54*cm,0.5*(9.0-0.25)*2.54*cm,0,2*pi);
      G4LogicalVolume* logicalPipeAir = new G4LogicalVolume(pipeAir,G4NISTair,"logicalPipeAir",0,0,0);
      G4ThreeVector airPos = G4ThreeVector(0,0,0.5*(1/8.0)*2.54*cm);
      G4PVPlacement *airWorld = new G4PVPlacement(0,airPos,"physicalPipeAir",logicalPipeAir,pipeWorld,false,0);
      logicalPipeAir->SetVisAttributes(airVis);

      G4Tubs* barrelAir0 = new G4Tubs("barrelAir0",0,0.5*(16.25-0.125)*2.54*cm,0.5*(2.0)*2.54*cm,0,2*pi);
      G4LogicalVolume* logicalBarrelAir0 = new G4LogicalVolume(barrelAir0,G4NISTair,"logicalBarrelAir0",0,0,0);
      airPos = G4ThreeVector(0,0,0.5*(21.0 - (1/16.0))*2.54*cm - 0.5*(2.0)*2.54*cm);
      G4PVPlacement *abovePipeAirWorld = new G4PVPlacement(0,airPos,"physicalBarrelAir0",logicalBarrelAir0,paraffinWorld,false,0);
      logicalBarrelAir0->SetVisAttributes(airVis);

      G4Tubs* barrelAir1 = new G4Tubs("barrelAir1",0.5*2.5*2.54*cm,0.5*(16.25-0.125)*2.54*cm,0.5*(1.5)*2.54*cm,0,2*pi);
      G4LogicalVolume* logicalBarrelAir1 = new G4LogicalVolume(barrelAir1,G4NISTair,"logicalBarrelAir1",0,0,0);
      airPos = G4ThreeVector(0,0,0.5*(21.0 - (1/16.0))*2.54*cm - (2.0)*2.54*cm - 0.5*(1.5)*2.54*cm );
      new G4PVPlacement(0,airPos,"physicalBarrelAir1",logicalBarrelAir1,paraffinWorld,false,0);
      logicalBarrelAir1->SetVisAttributes(airVis);

      G4Tubs* luciteRodInPipe = new G4Tubs("luciteRodInPipe",0.,0.5*(2.0)*2.54*cm,0.5*(3.5)*2.54*cm,0,2*pi);
      G4LogicalVolume* logicalLuciteRodInPipe = new G4LogicalVolume(luciteRodInPipe,G4NISTlucite,"logicalLuciteRodInPipe",0,0,0);
      G4ThreeVector lucitePos = G4ThreeVector(0,0,0.5*(9.0-0.25)*2.54*cm - 0.5*3.5*2.54*cm);
      new G4PVPlacement(0,lucitePos,"physicalLuciteRodInPipe",logicalLuciteRodInPipe,airWorld,false,0);
      logicalLuciteRodInPipe->SetVisAttributes(luciteVis);

      G4Tubs* luciteRodAbovePipe = new G4Tubs("luciteRodAbovePipe",0.,0.5*(2.0)*2.54*cm,0.5*(2.0)*2.54*cm,0,2*pi);
      G4LogicalVolume* logicalLuciteRodAbovePipe = new G4LogicalVolume(luciteRodAbovePipe,G4NISTlucite,"logicalLuciteRodAbovePipe",0,0,0);
      lucitePos = G4ThreeVector(0,0,0);
      new G4PVPlacement(0,lucitePos,"physicalLuciteRodAbovePipe",logicalLuciteRodAbovePipe,abovePipeAirWorld,false,0);
      logicalLuciteRodAbovePipe->SetVisAttributes(luciteVis);

      G4Tubs* luciteRodAboveBarrel = new G4Tubs("luciteRodAboveBarrel",0.,0.5*(2.0)*2.54*cm,0.5*(9.75)*2.54*cm,0,2*pi);
      G4LogicalVolume* logicalLuciteRodAboveBarrel = new G4LogicalVolume(luciteRodAboveBarrel,G4NISTlucite,"logicalLuciteRodAboveBarrel",0,0,0);
      lucitePos = G4ThreeVector(0,postPointsY[1]-0.5*16*2.54*cm+0.5*27.5*cm-2.25*2.54*cm+0.5*12*2.54*cm+0.5*16.25*2.54*cm,floorZ+(21)*2.54*cm+0.5*(9.75)*2.54*cm);
      new G4PVPlacement(0,lucitePos,"physicalLuciteRodAboveBarrel",logicalLuciteRodAboveBarrel,world,false,0);
      logicalLuciteRodAboveBarrel->SetVisAttributes(luciteVis);
    
    } //end pubeNaIParams.addBarrel if statement
    G4cout << "Position of Source: " << 0 << "," << postPointsY[1]-0.5*16*2.54*cm+0.5*27.5*cm-2.25*2.54*cm+0.5*12*2.54*cm+0.5*16.25*2.54*cm << "," << floorZ+(21-2.0-9.0+2.25)*2.54*cm << G4endl;
  }//end do R66

  if(pubeNaIParams.doOrb){

    //lead
    G4Orb* leadOrb = new G4Orb("leadOrb",pubeNaIParams.OrbRad);
    G4LogicalVolume* logicalLeadOrb = new G4LogicalVolume(leadOrb,shieldPbMat,"logicalLeadOrb",0,0,0);
    new G4PVPlacement(0,pubeNaIParams.OrbPos,"physicalleadOrb",logicalLeadOrb,world,false,0);
    logicalLeadOrb->SetVisAttributes(leadVis);

  } //end do Orb

} // ends PuBe SourceAndShield 

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
	G4double towerassy_shift = -6.5*cm; //FIXME hard-coded?
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
								point - (relative.unit()*shieldTestParams.sizethk/2.0), //put source at edge of shielding
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
void k100_DetectorConstruction::ConstructSimpleGammaCoin(G4VPhysicalVolume *world)
{
        //Get the coordinates of copies and source point and things
	G4ThreeVector stackrel_detorigin = zipParam->GetCoordinates(0);  //is the first one is 0
	G4double towerassy_shift = -6.5*cm; //FIXME hard-coded?
	G4double zipstack_shift = -(Tower_zPcut[0] - (zPldh[1] - zPldh[0]) - zPsdh[1] -Zip_z/2.0); //negative because tower is flipped
	G4double instack_shift = -stackrel_detorigin.z(); //negative because tower is flipped
	G4ThreeVector detorigin = G4ThreeVector(0,0,towerassy_shift + zipstack_shift +instack_shift);
	G4ThreeVector point(gammaCoinParams.xcntr,gammaCoinParams.ycntr,gammaCoinParams.zcntr);
	G4ThreeVector relative = point - detorigin;

	G4cout << "xdet: " << detorigin.x() << " ydet: " << detorigin.y() << " zdet: " << detorigin.z() << G4endl;
	G4cout << "towerassy_shift: " << towerassy_shift << " zipstack_shift: " << zipstack_shift << " instack_shift: " << instack_shift << G4endl;

	//Do the appropriate rotations
	G4RotationMatrix *gammadetrot = new G4RotationMatrix;
	gammadetrot->rotateZ(-relative.phi());
	gammadetrot->rotateY(-relative.theta());

	//create the simple rectangular shield
	G4Tubs* GeGammaCyl = new G4Tubs("GeGammaCyl_S",0.0,gammaCoinParams.sizer,gammaCoinParams.sizethk/2.0,0,2*pi);
	G4LogicalVolume* logicalGeGammaCyl;
        logicalGeGammaCyl = new G4LogicalVolume(GeGammaCyl,gammaCoinParams.coinmaterial,"GeGammaCyl_L",0,0,0);
	G4VPhysicalVolume* GeGammaCylWorld = new G4PVPlacement(gammadetrot, 
								point - (relative.unit()*gammaCoinParams.sizethk/2.0), //put source at edge of shielding
								"shieldBox_P",
								logicalGeGammaCyl,
								world,
								false,
								0);

	// Visualization attributes
	//G4VisAttributes* VisAttGeGammaCyl = new G4VisAttributes(G4Colour(128.0/255.0,0/255.,0/255.));
	G4cout << "Gamma Coincidence Detector is: " << gammaCoinParams.coinmaterial->GetName() << G4endl;
        G4VisAttributes* VisAttGeGammaCyl = new G4VisAttributes(G4Colour(255/255.,0/255.,0/255.));
	VisAttGeGammaCyl->SetForceWireframe(false);  //I want a Wireframe of the me
	logicalGeGammaCyl->SetVisAttributes(VisAttGeGammaCyl);  

	      //------------------------------------------------ 
        // Sensitive detectors
        //------------------------------------------------ 
    
        // Prepare to declare sensitive detectors
        G4SDManager* SDman = G4SDManager::GetSDMpointer();

        G4String detectorZipSDname = "gammaCoin1";
        G4int collID = -1; collID = SDman->GetCollectionID(detectorZipSDname);
        k100_ZipSD* azipSD1;
        ConstructGenericSensitiveInt=2; //?FIXME I actually forgot what role this is supposed to play 
        azipSD1 = new k100_ZipSD(detectorZipSDname, k100CollName.size()+1);
        k100CollName[detectorZipSDname] = k100CollName.size()+1;
        SDman->AddNewDetector(azipSD1);
        logicalGeGammaCyl->SetSensitiveDetector(azipSD1);
  
}
void k100_DetectorConstruction::ConstructPuBeNaI(G4VPhysicalVolume *world)
{
        //Get the coordinates of copies and source point and things
	G4ThreeVector stackrel_detorigin = zipParam->GetCoordinates(0);  //is the first one is 0
	G4double towerassy_shift = -6.5*cm; //FIXME hard-coded?
	G4double zipstack_shift = -(Tower_zPcut[0] - (zPldh[1] - zPldh[0]) - zPsdh[1] -Zip_z/2.0); //negative because tower is flipped
	G4double instack_shift = -stackrel_detorigin.z(); //negative because tower is flipped
	G4ThreeVector detorigin = G4ThreeVector(0,0,towerassy_shift + zipstack_shift +instack_shift);
	G4ThreeVector point(gammaCoinParams.xcntr,gammaCoinParams.ycntr,gammaCoinParams.zcntr);
	G4ThreeVector relative = point - detorigin;

	G4cout << "xdet: " << detorigin.x() << " ydet: " << detorigin.y() << " zdet: " << detorigin.z() << G4endl;
	G4cout << "towerassy_shift: " << towerassy_shift << " zipstack_shift: " << zipstack_shift << " instack_shift: " << instack_shift << G4endl;

	//Do the appropriate rotations
	G4RotationMatrix *gammadetrot = new G4RotationMatrix;
	gammadetrot->rotateZ(-relative.phi());
	gammadetrot->rotateY(-relative.theta());

	//create the simple rectangular shield
	G4Tubs* GeGammaCyl = new G4Tubs("GeGammaCyl_S",0.0,gammaCoinParams.sizer,gammaCoinParams.sizethk/2.0,0,2*pi);
	G4LogicalVolume* logicalGeGammaCyl;
        logicalGeGammaCyl = new G4LogicalVolume(GeGammaCyl,gammaCoinParams.coinmaterial,"GeGammaCyl_L",0,0,0);
	G4VPhysicalVolume* GeGammaCylWorld = new G4PVPlacement(gammadetrot, 
								point - (relative.unit()*gammaCoinParams.sizethk/2.0), //put source at edge of shielding
								"shieldBox_P",
								logicalGeGammaCyl,
								world,
								false,
								0);

	// Visualization attributes
	//G4VisAttributes* VisAttGeGammaCyl = new G4VisAttributes(G4Colour(128.0/255.0,0/255.,0/255.));
	G4cout << "Gamma Coincidence Detector is: " << gammaCoinParams.coinmaterial->GetName() << G4endl;
        G4VisAttributes* VisAttGeGammaCyl = new G4VisAttributes(G4Colour(255/255.,0/255.,0/255.));
	VisAttGeGammaCyl->SetForceWireframe(false);  //I want a Wireframe of the me
	logicalGeGammaCyl->SetVisAttributes(VisAttGeGammaCyl);  

	//------------------------------------------------ 
        // Sensitive detectors
        //------------------------------------------------ 
    
        // Prepare to declare sensitive detectors
        G4SDManager* SDman = G4SDManager::GetSDMpointer();

        G4String detectorZipSDname = "gammaCoin1";
        G4int collID = -1; collID = SDman->GetCollectionID(detectorZipSDname);
        k100_ZipSD* azipSD1;
        ConstructGenericSensitiveInt=2; //?FIXME I actually forgot what role this is supposed to play 
        azipSD1 = new k100_ZipSD(detectorZipSDname, k100CollName.size()+1);
        k100CollName[detectorZipSDname] = k100CollName.size()+1;
        SDman->AddNewDetector(azipSD1);
        logicalGeGammaCyl->SetSensitiveDetector(azipSD1);
  
}
void k100_DetectorConstruction::DetSizeMod(G4double R, G4double thk)
{
  Zip_Rout = R; 
  Zip_z = thk;
  return;
}

//spandey
void k100_DetectorConstruction::ConstructNaIArray(G4LogicalVolume*  logicalWorld) {

  G4VisAttributes* superVis = new G4VisAttributes(G4Colour(0./255.,255./255.,0/255.));
  superVis->SetForceSolid(false);
  G4VisAttributes* mother_NaI = new G4VisAttributes(G4Colour(0./255.,255./255.,0/255.,0.2));
  mother_NaI->SetForceSolid(false);
  G4double NaI_tile_length = 406 *mm;
  G4double NaI_tile_height = 57 *mm;
  G4double NaI_tile_width = 102 *mm;

  //Floor block
  //G4double fridgeHalfHeightToBottomPlate = (12.9045+13.25+0.25)*2.54*cm;
  G4double fridgeHalfHeightToBottomPlate = (12.9045+19.254+0.25)*2.54*cm; //modified 1/1/18 to get floor height right
  //G4double distanceToFloorZ = fridge_z+12.9045*2.54*cm - fridgeHalfHeightToBottomPlate - 21.0*2.54*cm;
  //G4double distanceCenterToFloor = fridgeHalfHeightToBottomPlate + 21.0*2.54*cm;
  G4double distanceCenterToFloor = fridgeHalfHeightToBottomPlate + 21.0*2.54*cm -70.86*mm; //compensate for 70.86mm discrepancy in floor distance 1/1/18
  G4double floorZ = fridge_z+12.9045*2.54*cm - distanceCenterToFloor;
  
  G4Box* NaI_tile = new G4Box("NaI_tile", 0.5 * NaI_tile_width, 0.5* NaI_tile_length, 0.5* NaI_tile_height);
  G4LogicalVolume* NaI_tile_LV = new G4LogicalVolume(NaI_tile,naiMat,"NaI_tile_LV",0,0,0);
  NaI_tile_LV->SetVisAttributes(superVis);

  /////////////////////////
  ///   Vertical Array  ///
  ////////////////////////

  G4Box* NaI_Mother_vertical = new G4Box("NaI_Mother_vertical", 0.5*(NaI_tile_width + 2*cm), 0.5*(NaI_tile_length + 2*cm) , 0.5*(1255 + 2*cm));
  G4LogicalVolume* NaI_Mother_vertical_LV = new G4LogicalVolume(NaI_Mother_vertical,G4NISTair,"NaI_Mother_vertical_LV",0,0,0);
  NaI_Mother_vertical_LV->SetVisAttributes(mother_NaI);
  //G4ThreeVector tilePos(fridge_x + 175.1*mm, fridge_y, fridge_z);
  G4double xposition = (-1*frame_x) + (0.5 * NaI_tile_width) + 10*cm; 
  G4double zposition = fridge_z-19.254*2.54*cm+2.0*2.54*cm; //floorZ + 89 + (1255 .0+ 2*cm)/2;
  zposition = zposition + 1255.0/2 - 89.0;
  //Zposition = floorZ + 89 + (fTileHeight)/2;
  G4ThreeVector NaI_Mother_vertical_pos(-1.0*xposition , 0, zposition);
  new G4PVPlacement(0,NaI_Mother_vertical_pos,"NaI_Mother_vertical_placement",NaI_Mother_vertical_LV,physicalWorld,false,0);


  G4double stack_height = 292;
  G4double stack_width = NaI_tile_width;
  G4double stack_length = NaI_tile_length;
  G4int nTiles = 20;
  k100_NaIParametrization *NaIParam = new k100_NaIParametrization(stack_height, stack_length, stack_width,
                                                                  NaI_tile_height, NaI_tile_length, NaI_tile_width,
                                                                  3,nTiles);
  
  

  new G4PVParameterised("NaIArray_stack1",NaI_tile_LV,NaI_Mother_vertical_LV,kZAxis, nTiles, NaIParam);
  //std::cout<<"Mother VOLUME x position = "<<xposition<<std::endl;
  //G4cout<<"xxxx NaI cooredinates of first copy = "<<NaIParam->GetCoordinates(0)<<G4endl;

  

  /////////////////////////
  ///   Bottom Array  ///
  ////////////////////////

  G4Box* NaI_Mother_bottom = new G4Box("NaI_Mother_bottom", 0.5*(324.0 + 2*cm), 0.5*(NaI_tile_length + 2*cm) , 
    0.5*(NaI_tile_height + 2*cm));
  G4LogicalVolume* NaI_Mother_bottom_LV = new G4LogicalVolume(NaI_Mother_bottom,G4NISTair,"NaI_Mother_bottom_LV",0,0,0);
  NaI_Mother_bottom_LV->SetVisAttributes(mother_NaI);
  //G4ThreeVector tilePos(fridge_x + 175.1*mm, fridge_y, fridge_z);
  xposition = 0.0; 
  zposition = fridge_z-19.254*2.54*cm+2.0*2.54*cm; //floorZ + 89 + (1255 .0+ 2*cm)/2;
  zposition = zposition  - 166.0;
  //std::cout<<"zposition for bottom NaI array = "<<zposition<<std::endl;
  //Zposition = floorZ + 89 + (fTileHeight)/2;
  G4ThreeVector NaI_Mother_bottom_pos(xposition , 0, zposition);
  G4RotationMatrix r90Rotation;    // Rotate bottom tiles
  r90Rotation.rotateZ(90.*deg);
  G4Transform3D BottomRotate(r90Rotation, NaI_Mother_bottom_pos);
  //new G4PVPlacement(0,NaI_Mother_bottom_pos,"NaI_Mother_bottom_placement",NaI_Mother_bottom_LV,physicalWorld,false,0);
  new G4PVPlacement(BottomRotate,"NaI_Mother_bottom_placement",NaI_Mother_bottom_LV,physicalWorld,false,0);


  
  stack_height = NaI_tile_height;
  stack_width = 324;
  stack_length = NaI_tile_length;
  nTiles = 3;
  k100_NaIParametrization *NaIParam_bottom = new k100_NaIParametrization(stack_height, stack_length, stack_width,
                                                                  NaI_tile_height, NaI_tile_length, NaI_tile_width,
                                                                  1,nTiles);



  new G4PVParameterised("NaIArray_stack2",NaI_tile_LV,NaI_Mother_bottom_LV,kXAxis, nTiles, NaIParam_bottom);



  //------------------------------------------------ 
  // Sensitive detectors
  //------------------------------------------------ 

  
  // Prepare to declare sensitive detectors
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String detectorZipSDname = "NaITile";
  G4int collID = -1; collID = SDman->GetCollectionID(detectorZipSDname);
  k100_ZipSD* azipSD1;
  ConstructGenericSensitiveInt=2; //?FIXME I actually forgot what role this is supposed to play 
  azipSD1 = new k100_ZipSD(detectorZipSDname, k100CollName.size()+1);
  k100CollName[detectorZipSDname] = k100CollName.size()+1;
  SDman->AddNewDetector(azipSD1);
  NaI_tile_LV->SetSensitiveDetector(azipSD1);

  

}



































