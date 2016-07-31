#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4UserDetectorConstruction.hh"
#include "Globals.hh"
#include "G4RotationMatrix.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4VSensitiveDetector;
class G4VisAttributes;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction ();
    virtual ~DetectorConstruction ();

  public:
    virtual G4VPhysicalVolume* Construct ();

  private:
    G4Material* air;
    G4Material* argongas;
    G4Material* aluminium;
    G4Material* copper;
    G4Material* lead;

    G4Material* scintillator;
    G4Material* C2H4;
    G4Material* C3H6;
    G4Material* CO2;
    G4Material* ArCO2;
    G4Material* CsI;

    G4Material* GEM;

    G4VisAttributes* worldVisAtt;
    G4VisAttributes* hodoscopeVisAtt;
    G4VisAttributes* chamberVisAtt;

    G4double deltaStacks; // difference between upper and lower stack where object stays
    G4double deltaDet;    // gap between each consecutive detector in either upper or lower
    G4double detSide1;    // must be length side
    G4double detSide2;    // must be breadth side
    G4double detSide3;    // must be length side
    G4double detSide4;    // must be breadth side
    G4double detHeight;

    G4double objLen;
    G4double objWid;
    G4double objHt;

    G4Int numDet;        // total number of detectors 3 top and 3 bottom
    G4Int numChanPerDet; // total number of channels per detectors
}
#endif
