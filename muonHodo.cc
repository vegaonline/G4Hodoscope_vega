//
// **************************************************************
// * Trial for making hodoscope simulation for muon
// *  - Abhijit Bhattacharyya
// **************************************************************

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "Randomize.hh"

#include "DetectorConstruction.hh"
// #include "PhysicsList.hh"
// #include "PrimaryGeneratorAction.hh"
// #include "EventAction.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

int main ( int argc, char** argv )
{
    // Detect mode and define UI session
    G4UIExecutive* ui = 0;
    if ( argc == 1 ) {
        ui = new G4UIExecutive ( argc, argv );
    }

    // Construct run manager
    G4RunManager* runManager = new G4RunManager;

    // mandatory initialization classes with runManager
    runManager->SetUserInitialization ( new DetectorConstruction() );
    // runManager->SetUserInitialization(new PhysicsList);

    // User Action Initialization classes <-- I have seen that this is optional
    // runManager->SetUserInitialization ( new ActionInitialization () );

    // Initialize Visualization
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize ();

    // Pointer to the User Interface
    G4UImanager* UImanager = G4UImanager::GetUIpointer ();

    // Starting UI or macro
    if ( !ui ) {
        G4String command = "/control/execute";
        G4String filename = argv[1];
        UImanager->ApplyCommand ( command + filename );
    } else {
        UImanager->ApplyCommand ( "/control/execute vis.mac" );
        ui->SessionStart ();
        delete ui;
    }

    // Code terminates releasing stores
    delete visManager;
    delete runManager;

    return 0;
}
