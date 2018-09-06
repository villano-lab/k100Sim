// ------------------------------------------------
//
// k100_NR_DetectorParameterDef.hh 
//
// ------------------------------------------------

// -- Detector box --
G4double microHair;

// "Copper Upper Tower" Parameters
static const G4int nZcut=2;
G4double zPcut[nZcut], rIcut[nZcut], rOcut[nZcut];

// Lower cap of the "copper upper tower" Parameters
static const G4int nZclc=2;
G4double zPclc[nZclc], rIclc[nZclc], rOclc[nZclc]; 

// Connector Tube Parameters
G4double rIctu, rOctu, zHctu;

// Base of Connector Tube Parameters
G4double rIrb, rOrb, zHrb;

// Upper cap of the detector housing Parameters
static  const G4int nZudh=2;
G4double zPudh[nZudh], rIudh[nZudh], rOudh[nZudh];

// The side detector housing
static  const G4int nZsdh=2;
G4double zPsdh[nZsdh], rIsdh[nZsdh], rOsdh[nZsdh];

// The lower cap of the detector housing
static  const G4int nZldh=2;
G4double zPldh[nZldh], rIldh[nZldh], rOldh[nZldh];

// Side Coax
G4double sdcx_widthH, sdcx_thicknessH;
G4double sdcx_lenH[6];
G4double sdcx_radius;

// Tower Parameters
static const G4int Tower_nZcut = 2;
G4double Tower_zPcut[Tower_nZcut];
G4double Tower_rIcut[Tower_nZcut];
G4double Tower_rOcut[Tower_nZcut];

G4double DetBoxShim;
G4double DetBox_Rout, DetBox_z;

G4double Tower_Center_Z; 
G4double Zip_Househeight;
G4double Zip_Space;
///////Tower_Bottom_Z, 

G4double xtow[5], ytow[5], ztow[5];

// Zip Parameters
G4double Zip_Rout, Zip_z;
G4double Zip_Flat_R1, Zip_Flat_R2;

