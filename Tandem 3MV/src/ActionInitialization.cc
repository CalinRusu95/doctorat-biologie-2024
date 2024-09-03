#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"

#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

ActionInitialization::ActionInitialization(RunInput* rinp, CLHEP::HepRandomEngine* rgen) : G4VUserActionInitialization(), runInput(rinp), randGen(rgen) {}

ActionInitialization::~ActionInitialization() {}

void ActionInitialization::Build() const {
	SetUserAction(new PrimaryGeneratorAction(runInput, randGen));
	RunAction* runAction = new RunAction(runInput);
	SetUserAction(runAction);

	EventAction* eventAction = new EventAction(runAction);
	SetUserAction(eventAction);
	SetUserAction(new SteppingAction(eventAction, runInput));
}
