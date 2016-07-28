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













DetectorConstruction::DetectorConstruction()
	: air(0), argongas(0), scintillator(0), CsI(0), lead(0), 
	detSide1(160.0*cm), detSide2(90.0*cm), detSide3(160.0*cm), detSide4(60.0*cm),
	detHeight(5.0*cm), deltaDet(3.0*cm), deltaStacks(45.0*cm), 
	objLen(10.0*cm), objWid(10.0*cm), objHt(10.0cm),
	  worldVisAtt(0), hodoscopeVisAtt(0), chamberVisAtt(0)
{


}

DetectorConstruction::~DetectorConstruction() {
	DestroyMaterials();

	delete worldVisAtt;
	delete hodoscopeVisAtt;
	delete chamberVisAtt;
}

G4VPhysicalVolume* DetectorConstruction::Construct() {
	G4VSensitiveDetector* hodoscope;
	G4VSensitiveDetector* chamber1;

	ConstructMaterials();

	// geometry ******
	G4double hodoLen = 1.5 * max(detSide1, detSide3);
	G4double hodoWid = 1.5 * max(detSide3, detSide4);
	G4double hodoHt  = 6 * deltaDet + deltaStacks + 6 * detHeight; 

	
	// Experimental Hall (world volume)
	G4VSolid* worldSolid = new G4Box("worldBox", 2.0*m, 2.0*m, 5.0*m);
	G4LogicalVolume* worldLogical = new G4LogicalVolume(worldSolid, air, "worldLogical", 0,0,0);
	G4VPhysicalVolume* worldPhysical = new G4PVPlacement(0, G4ThreeVector(), worldLogical, "worldPhysical", 0, 0, 0);

	// Place the hodoscope in the world volume
	G4VSolid* hodoscopeSolid = new G4Box("hodoscopeBox", hodoLen, hodoWid, hodoHt);
	G4LogicalVolume* hodoscopeLogical = new G4LogicalVolume(hodoscopeSolid, 






}













