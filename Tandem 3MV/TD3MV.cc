#include "Randomize.hh"
#include <ctime>
#include <algorithm>

#include "RunInput.hh"
#include "PhysicsList.hh"
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

// Uncomment this if using multi-threading
// #ifdef G4MULTITHREADED
// #include "G4MTRunManager.hh"
// #else
#include "G4RunManager.hh"
// #endif

#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

int main(int argc, char** argv) {
    // Initialize the RunInput object to handle input configurations
    RunInput* runInput = new RunInput();
    G4bool VisualOn = runInput->IsVisualOn();

    // Set up the random number generator
    CLHEP::HepRandomEngine* randGen = new CLHEP::RanecuEngine;
    long int seed = static_cast<long int>(time(NULL));   // Get seed from current time
    G4Random::setTheEngine(randGen);                     // Choose the random engine
    G4Random::setTheSeed(seed);

    // Get the number of events to simulate
    G4int numberOfEvent = runInput->GetNumEvents();

    // Construct the default run manager
    // Uncomment this section if using multi-threading
    // #ifdef G4MULTITHREADED
    //     G4MTRunManager* runManager = new G4MTRunManager;
    //     runManager->SetNumberOfThreads(8);
    //     G4cout << "Running in multi-threading mode!" << G4endl;
    // #else
    G4RunManager* runManager = new G4RunManager;
    G4cout << "Running in sequential mode!" << G4endl;
    // #endif

    // Set mandatory initialization classes
    runManager->SetUserInitialization(new DetectorConstruction(runInput));
    runManager->SetUserInitialization(new PhysicsList());
    runManager->SetUserInitialization(new ActionInitialization(runInput, randGen));

    // Initialize the G4 kernel
    runManager->Initialize();

    // Set up the visualization manager
    G4VisManager* visManager = new G4VisExecutive();
    if (VisualOn) {
        visManager->Initialize();

        // Interactive mode: define UI session
        G4UIExecutive* UIexec = new G4UIExecutive(argc, argv);
        UIexec->SessionStart();
        delete UIexec;
    } else {
        // Uncomment this section if using batch mode
        // G4UImanager* UI = G4UImanager::GetUIpointer();
        // UI->ApplyCommand("/run/verbose 1");
        // UI->ApplyCommand("/event/verbose 1");
        // UI->ApplyCommand("/tracking/verbose 1");
        // UI->ApplyCommand("/physics/verbose 1");

        // Batch mode: run the simulation
        G4cout << G4endl << "***Running " << numberOfEvent 
               << " events with random seed " << seed << G4endl << G4endl;
        runManager->BeamOn(numberOfEvent);
    }

    // Clean up and release resources
    delete runInput;
    delete visManager;
    delete runManager;

    return 0;
}
