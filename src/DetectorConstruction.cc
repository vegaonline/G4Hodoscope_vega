#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4UserLimits.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

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
    , detSide1 (0 )
    , detSide2 ( 0 )
    , detSide3 ( 0 )
    , detSide4 ( 0 )
    , detHeight ( 0 )
    , deltaDet ( 0 )
    , deltaStacks ( 0 )
    , objLen ( 0 )
    , objWid ( 0 )
    , objHt ( 0 )
    , numDet ( 6 )
    , numChanPerDet ( 32 )
    , worldVisAtt ( 0 )
    , hodoVisAtt ( 0 )
    , GEMVisAtt ( 0 )
{
//  ConstructMaterials();
}

DetectorConstruction::~DetectorConstruction ()
{
    DestroyMaterials();

    delete worldVisAtt;
    delete hodoVisAtt;
    delete GEMVisAtt;
}

void DetectorConstruction::ConstructMaterials ()
{
    G4double a;
    G4double z;
    G4double density;
    G4double weightRatio;
    G4double fracMass;
    G4String name;
    G4String symbol;
    G4int nElem;
    G4int nComp;

    // Argon gas
    a = 39.95*g/mole;
    density = 1.782e-03*g/cm3;
    argongas = new G4Material ( name = "ArgonGas", z = 18, a, density );

    // Elements for micture and compounds
    a = 1.01 * g / mole;
    G4Element* elH = new G4Element ( name = "Hydrogen", symbol = "H", z = 1.0, a );

    a = 12.01 * g / mole;
    G4Element* elC = new G4Element ( name = "Carbon", symbol = "C", z = 6.0, a );

    a = 14.01 * g / mole;
    G4Element* elN = new G4Element ( name = "Nitrogen", symbol = "N", z = 7.0, a );

    a = 16.00 * g / mole;
    G4Element* elO = new G4Element ( name = "Oxygen", symbol = "O", z = 8.0, a );

    // Air
    density = 1.032 * g / cm3;
    air = new G4Material ( name = "Air", density, nElem = 2 );
    air->AddElement ( elN, fracMass = 70.0 * perCent );
    air->AddElement ( elO, fracMass = 30.0 * perCent );

    // CO2
    density = 1.8714e-3 * g / cm3;
    CO2 = new G4Material ( name = "CO2", density, nElem = 2 );
    CO2->AddElement ( elC, 1 );
    CO2->AddElement ( elO, 2 );

    // Ar-CO2 mixture
    density = 1.77e-3 * g / cm3;
    ArCO2 = new G4Material ( name = "ArCO2", density, nComp = 2 );
    ArCO2->AddMaterial ( argongas, fracMass = 90.0 * perCent );
    ArCO2->AddMaterial ( CO2, fracMass = 10.0 * perCent );

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
    G4Element* elI = new G4Element ( name = "Iodine", symbol = "I", z = 53.0, a );
    a = 132.9 * g / mole;
    G4Element* elCs = new G4Element ( name = "Cesium", symbol = "Cs", z = 55.0, a );
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
    G4MaterialTable* matTable = (G4MaterialTable*)G4Material::GetMaterialTable ();
    for ( index = 0; index < matTable->size (); index++ ) {
      delete (*(matTable))[index];
    }
    matTable->clear ();
    G4ElementTable* elemTable = (G4ElementTable*)G4Element::GetElementTable ();
    for ( index = 0; index < elemTable->size (); index++ ) {
        delete ( *( elemTable ) )[index];
    }
    elemTable->clear ();
}

void DetectorConstruction::DumpGeometricalTree(G4VPhysicalVolume* aVolume,G4int depth)
{
  for(int isp=0;isp<depth;isp++)
  { G4cout << "  "; }
  G4cout << aVolume->GetName() << "[" << aVolume->GetCopyNo() << "] "
         << aVolume->GetLogicalVolume()->GetName() << " "
         << aVolume->GetLogicalVolume()->GetNoDaughters() << " "
         << aVolume->GetLogicalVolume()->GetMaterial()->GetName();
  if(aVolume->GetLogicalVolume()->GetSensitiveDetector())
  {
    G4cout << " " << aVolume->GetLogicalVolume()->GetSensitiveDetector()
                            ->GetFullPathName();
  }
  G4cout << G4endl;
  for(int i=0;i<aVolume->GetLogicalVolume()->GetNoDaughters();i++)
  { DumpGeometricalTree(aVolume->GetLogicalVolume()->GetDaughter(i),depth+1); }
}


G4VPhysicalVolume* DetectorConstruction::Construct ()
{
    G4VSensitiveDetector* hodoscope;
    
    G4bool checkOverlaps = true;    

     ConstructMaterials();

    // geometry ******
    detSide1 = 160.0*cm; detSide2 = 90.0*cm; detSide3 = 160.0*cm; detSide4 = 60.0*cm;
    deltaDet = 3.0*cm; deltaStacks = 45.0*cm; 

    G4double GEMCathodeHt = 0.3*cm;   G4double GEMDriftHt = 0.3*cm; 
    G4double G1G2gap = 0.1*cm;  G4double G2G3gap = 0.2*cm; G4double G3readergap = 0.1*cm;
    G4double GEMreaderHt = 0.3*cm;
    G4double GEMCuHt = 0.5e-3*cm; G4double GEMPolyHt = 5.0e-3*cm;  G4double GEMHt = GEMPolyHt + 2.0 * GEMCuHt;
    detHeight = GEMCathodeHt + GEMDriftHt + G1G2gap + G2G3gap + G3readergap + GEMreaderHt + GEMHt;

    objLen = 10.0*cm; objWid = 10.0*cm; objHt = 10.0*cm;

    G4double worldx = 2.0 * m; // effectively world dim = 4m X 4m X 10m
    G4double worldy = 2.0 * m;
    G4double worldz = 5.0 * m;

    G4double hodoLen = 1.5 * std::max ( detSide1, detSide3 );
    G4double hodoWid = 1.5 * std::max ( detSide3, detSide4 );
    G4double hodoHt = deltaStacks + numDet * deltaDet + numDet * detHeight;

    // Experimental Hall (world volume)
    G4VSolid* worldSolid = new G4Box ( "worldBox", worldx, worldy, worldz );
    G4LogicalVolume* worldLogical = new G4LogicalVolume ( worldSolid, air, "worldLogical", 0, 0, 0 );
    G4VPhysicalVolume* worldPhysical
      = new G4PVPlacement ( 0, G4ThreeVector(), worldLogical, "worldPhysical", 0, 0, 0, checkOverlaps );
 
    // Place the hodoscope in the world volume
    G4VSolid* hodoscopeSolid = new G4Box ( "hodoscopeBox", 0.5 * hodoLen, 0.5 * hodoWid, 0.5 * hodoHt );
    G4LogicalVolume* hodoscopeLV = new G4LogicalVolume ( hodoscopeSolid, air, "hodoscopeogical", 0, 0, 0 );
    new G4PVPlacement ( 0, G4ThreeVector ( 0., 0., 0. ), hodoscopeLV, "hodoscopePhysical", worldLogical, 0, 0, checkOverlaps );

    // GEM one layer Cu - C2H4 - Cu
    G4double dx1Trd = 0.5 * detSide2;
    G4double dx2Trd = 0.5 * detSide4;
    G4double dzTrd = 0.5 * detSide1;
    
    G4VSolid* gemCuTop = new G4Trd ( "CuTop", dx1Trd, dx2Trd, 0.5 * GEMCuHt, 0.5 * GEMCuHt, dzTrd );
    G4LogicalVolume* gemCuTopLV = new G4LogicalVolume ( gemCuTop, copper, "gemCuTopLV", 0, 0, 0 );
    G4VSolid* gemPoly = new G4Trd ( "Poly", dx1Trd, dx2Trd, 0.5 * GEMPolyHt, 0.5 * GEMPolyHt, dzTrd);
    G4LogicalVolume* gemPolyLV = new G4LogicalVolume ( gemPoly, C2H4, "gemPolyLV", 0, 0, 0 );
    G4VSolid* gemCuBot = new G4Trd ( "CuBot", dx1Trd, dx2Trd, 0.5 * GEMCuHt, 0.5 * GEMCuHt, dzTrd);
    G4LogicalVolume* gemCuBotLV = new G4LogicalVolume ( gemCuBot, copper, "gemCuBotLV", 0, 0, 0 );
    
    // Drift chamber
    G4VSolid* driftCuTop = new G4Trd("DriftCuTop", dx1Trd, dx2Trd, 0.5*GEMCuHt, 0.5*GEMCuHt, dzTrd);
    G4LogicalVolume* driftCuTopLV = new G4LogicalVolume(driftCuTop, copper, "driftCuTop", 0, 0, 0);
    
    G4VSolid* driftChamber = new G4Trd("Drift", dx1Trd, dx2Trd, 0.5 * GEMDriftHt, 0.5 * GEMDriftHt, dzTrd);
    G4LogicalVolume* driftChamberLV = new G4LogicalVolume(driftChamber, ArCO2, "Drift", 0, 0, 0);
    
    G4VSolid* readerCuBot = new G4Trd("readerCuBot", dx1Trd, dx2Trd, 0.5*GEMCuHt, 0.5*GEMCuHt, dzTrd);
    G4LogicalVolume* readerCuBotLV = new G4LogicalVolume(readerCuBot, copper, "readerCuBot", 0, 0, 0);

    // Construction of single GEM
    G4RotationMatrix Ritem;    
    G4ThreeVector Titem;

    G4VSolid* GEMDet = new G4Trd("GemDet", dx1Trd, dx2Trd, 0.5 * detHeight, 0.5 * detHeight, dzTrd);
    G4LogicalVolume* GEMDetLV = new G4LogicalVolume(GEMDet, ArCO2, "GEMDetector");
    Ritem.rotateX(90.0*deg); Ritem.rotateZ(90.0*deg); 
    Titem.setX(0.0); Titem.setY(0.0); Titem.setZ(0.5*hodoHt);
    new G4PVPlacement(G4Transform3D(Ritem, Titem), GEMDetLV, "GEMdetector", hodoscopeLV, 0, 0, checkOverlaps);
    
    // Place drift board, GEM foils and readout in one GEM
    G4RotationMatrix Ra; G4ThreeVector Ta;
    Ta.setX(0.0); Ta.setY(-0.5*detHeight); Ta.setZ(0.0);
    new G4PVPlacement(0, Ta, driftCuTopLV, "driftPhysical", GEMDetLV, 0, 0, checkOverlaps);
    

    // Setting Visual Attribute
    worldVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    worldVisAtt->SetVisibility(false);
    worldLogical->SetVisAttributes(worldVisAtt);
    
    hodoVisAtt = new G4VisAttributes(G4Colour(0.8888, 0.0, 0.0));
    hodoVisAtt->SetForceWireframe(true);
    hodoscopeLV->SetVisAttributes(hodoVisAtt);

    

    // return the world physical volume ----------------------------------------

    G4cout << G4endl << "The geometrical tree defined are : " << G4endl << G4endl;
    DumpGeometricalTree(worldPhysical);

    return worldPhysical;



}
