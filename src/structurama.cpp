#include "interface.h"
#include "nexusfile.h"
#include "st.h"

#ifndef MAC_GUI

int main (int argc, char const *argv[]) {

	// instantiate the user interface
	NexusFile *nf = new NexusFile();
	CLIUserInterface *cliUserInterface = new CLIUserInterface();
	CompInterface compInterface(cliUserInterface, nf);

	// start the run loop
	compInterface.readLoop();

	// free up memory
	delete cliUserInterface;
	delete nf;
	
  return 0;
}

#endif
