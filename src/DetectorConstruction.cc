#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

DetectorConstruction::DetectorConstruction ()
    : air ( 0 )
    , argongas ( 0 )
    , scintillator ( 0 )
    , C3H6 ( 0 )
    , C2H4 ( 0 )
    , CO2 ( 0 )
    , ArCO2 ( 0 )
    , CsI ( 0 )
    , aluminium ( 0 )
    , copper ( 0 )
    , lead ( 0 )
    , detSide1 ( 160.0 * cm )
    , detSide2 ( 90.0 * cm )
    , detSide3 ( 160.0 * cm )
    , detSide4 ( 60.0 * cm )
    , detHeight ( 5.0 * cm )
    , deltaDet ( 3.0 * cm )
    , deltaStacks ( 45.0 * cm )
    , objLen ( 10.0 * cm )
    , objWid ( 10.0 * cm )
    , objHt ( 10.0cm )
    , numDet ( 6 )
    , numChanPerDet ( 32 )
    , worldVisAtt ( 0 )
    , hodoscopeVisAtt ( 0 )
    , chamberVisAtt ( 0 )
{
}

DetectorConstruction::~DetectorConstruction ()
{
    DestroyMaterials ();

    delete worldVisAtt;
    delete hodoscopeVisAtt;
    delete chamberVisAtt;
}

void DetectorConstruction::ConstructMaterials ()
{
    G4double a;
    G4double z;
    G4double density;
    G4double weightRatio;
    G4String name;
    G4String symbol;
    G4int nElem;

    // Argon gas
    a = 39.95 * g / mole;
    argongas = new G4Material ( name = "ArgonGas", z = 18, a, density );

    // Elements for micture and compounds
    a = 1.01 * g / mole;
    G4Element* elH = new G4Element ( name = "Hydrogen", symbol = "H", z = 1.0, a );

    a = 12.01 * g / mole;
    G4Element* elC = new G4Elelment ( name = "Carbon", symbol = "C", z = 6.0, a );

    a = 14.01 * g / mole;
    G4Elelment* elN = new G4Element ( name = "Nitrogen", symbol = "N", z = 7.0, a );

    a = 16.00 * g / mole;
    G4Elelment* elO = new G4Element ( name = "Oxygen", symbol = "O", z = 8.0, a );

    // Air
    density = 1.032 * g / cm3;
    air = new G4Material ( name = "Air", density, nElem = 2 );
    air->AddElement ( elN, fracMass = 70 * percent );
    air->AddElement ( elO, fracMass = 30 * percent );

    // CO2
    density = 1.8714e-3 * g / cm3;
    CO2 = new G4Material ( name = "CO2", density, nElem = 2 );
    CO2->AddElement ( elC, 1 );
    CO2->AddElement ( elO, 2 );

    // Ar-CO2 mixture
    density = 1.77e-3 * gm / cm3;
    ArCO2 = new G4Material ( name = "Ar-CO2", density, ncomp = 2 );
    ArCO2->AddElement ( argongas, fracMass = 90 * percent );
    ArCO2->AddMaterial ( CO2, fracMass = 10 * percent );

    // scintillator
    density = 1.032 * g / cm3;
    scintillator = new G4Material ( name = "Scintillator", density, nElem = 2 );
    scintillator->AddElement ( elC, 9 );
    scintillator->AddElement ( elH, 10 );

    // polypropelene
    density = 0.946 * g / cm3;
    C3H6 = new G4Material ( name = "Polypropelene", density, nElem = 2 );
    C3H6->AddElement ( elC, 3 );
    C3H6->AddElement ( elH, 6 );

    // polyethylene
    density = 0.94 * g / cm3;
    C2H4 = new G4Material ( name = "Polyethylene", density, nElem = 2 );
    C2H4->AddElement ( elC, 2 );
    C2H4->AddElement ( elH, 4 );

    // Aluminium
    a = 26.982 * g / mole;
    density = 2.7 * g / cm3;
    aluminium = new G4Material ( name = "Al", z = 13.0, a, density );

    // Copper
    a = 63.546 * g / mole;
    density = 8.94 * g / cm3;
    copper = new G4Material ( name = "Copper", z = 29, a, density );

    // Lead
    a = 207.19 * g / mole;
    density = 11.35 * g / cm3;
    lead = new G4Material ( name = "Lead", z = 82.0, a, density );

    // CsI
    a = 126.9 * g / mole;
    G4Elelment* elI = new G4Element ( name = "Iodine", symbol = "I", z = 53.0, a );
    a = 132.9 * g / mole;
    G4Elelment* elCs = new G4Element ( name = "Cesium", symbol = "Cs", z = 55.0, a );
    density = 4.51 * g / cm3;
    CsI = new G4Material ( name = "CsI", density, nElem = 2 );
    CsI->AddElement ( elI, weightRatio = 0.5 );
    CsI->AddElement ( elCs, weightRatio = 0.5 );

    G4cout << G4endl << "The materials defined are: " << G4endl << G4endl;
    G4cout << *( G4Material::GetMaterialTable () ) << G4endl;
}

void DetectorConstruction::DestroyMaterials ()
{
    // Destroy all allocated elements and materials
    size_t index;
    G4MaterialTable* matTable = (G4MaterialTable)G4Material::GetMaterialTable ();
    for ( index = 0; index < matTable->size (); index++ ) {
        delete (*matTable))[index];
    }
    matTable->clear ();
    G4ElementTable* elemTable = (G4ElementTable*)G4Element::GetElementTable ();
    for ( index = 0; index < elemTable->size (); index++ ) {
        delete ( *( elemTable ) )[index];
    }
    elemTable->clear ();
}

G4VPhysicalVolume* DetectorConstruction::Construct ()
{
    G4VSensitiveDetector* hodoscope;
    G4VSensitiveDetector* chamber1;

    ConstructMaterials ();

    // geometry ******
    G4double worldx = 2.0 * m; // effectively world dim = 4m X 4m X 10m
    G4double worldy = 2.0 * m;
    G4double worldz = 5.0 * m;

    G4double hodoLen = 1.5 * max ( detSide1, detSide3 );
    G4double hodoWid = 1.5 * max ( detSide3, detSide4 );
    G4double hodoHt = deltaStacks + numDet * deltaDet + numDet * detHeight;

    // Experimental Hall (world volume)
    G4VSolid* worldSolid = new G4Box ( "worldBox", worldx, worldy, worldz );
    G4LogicalVolume* worldLogical = new G4LogicalVolume ( worldSolid, air, "worldLogical", 0, 0, 0 );
    G4VPhysicalVolume* worldPhysical
        = new G4PVPlacement ( 0, G4ThreeVector (), worldLogical, "worldPhysical", 0, 0, 0 );

    // Place the hodoscope in the world volume
    G4VSolid* hodoscopeSolid = new G4Box ( "hodoscopeBox", 0.5 * hodoLen, 0.5 * hodoWid, 0.5 * hodoHt );
    G4LogicalVolume* hodoscopeLogical = new G4LogicalVolume ( hodoscopeSolid, air, "hodoscopeogical", 0, 0, 0 );
    new G4PVPlacement ( 0, G4ThreeVector ( 0., 0., 0. ), hodoscopeLogical, "hodoscopePhysical", worldLogical, 0, 0 );

    // GEM one layer Cu - CH2 - Cu
    G4VSolid* gemCu1 = new G4Box ( "Cu1", )
}
