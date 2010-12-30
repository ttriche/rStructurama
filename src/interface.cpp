#include "interface.h"
#include "st.h"
#include "nexusfile.h"
#include <sstream>
#include <ctime>
#include <string>

using namespace std;




CLIUserInterface::CLIUserInterface(void) {

	stprint("                                                                                          \n"); 
	stprint("                                            Structurama 2.0                               \n\n"); 
	stprint("                                                  by                                      \n\n"); 
	stprint("              John P. Huelsenbeck (1), Edna T. Huelsenbeck (1), and Peter Andolfatto (2)  \n\n"); 
	stprint("                                (1) Department of Integrative Biology                     \n"); 
	stprint("                                  University of California, Berkeley                      \n\n"); 
	stprint("                           (2) Department of Ecology & Evolutionary Biology               \n"); 
	stprint("                                         Princeton University                             \n\n"); 
#	if !defined(MAC_GUI)
	stprint("                                     Type \"Help\" to get started.                        \n\n"); 
#	endif
}



CLIUserInterface::~CLIUserInterface(void) {

	time_t currentTime;
	time (&currentTime);
	struct tm * ptm= localtime(&currentTime);
	int hour = ptm->tm_hour;
	string salutation = "";
	if ( hour >= 21 || hour < 3 )
		salutation = "night";
	else if ( hour >= 3 && hour < 12 )
		salutation = "morning";
	else if ( hour >= 12 && hour < 17 )
		salutation = "afternoon";
	else
		salutation = "evening";
	stprint("   Thank you for using Structurama and have a nice %s.\n\n", salutation.c_str());
}

const string CompInterface::commandNames[] = { "acknowledgments", "calcmarginal", "citations", "execute", "help", "mcmc", "model", "readsamples", "set", "showdata", "showmeanpart", "shownumpops", "showtogetherness", "quit" };
const vector<string> CompInterface::commands(commandNames,commandNames+numCommands());
int CompInterface::numCommands() { return sizeof(commandNames)/sizeof(string); }

CompInterface::CompInterface(UserInterface *ui, NexusFile *nf) { 

	userInterface = ui; 
	nexusFile = nf;
	nexusFile->setCompInterfacePtr(this);
}

void CompInterface::readLoop(void) {

	bool wasPreviousCommandBlank = false;
	for(;;) 
		{
		try 
			{
			string command;
			if (wasPreviousCommandBlank == false)
				userInterface->prompt("Structurama > ");
			if ( userInterface->readCommand(command) ) 
				{
				if ( MyString::isBlank(command) == true )
					wasPreviousCommandBlank = true;
				else
					{
					wasPreviousCommandBlank = false;
					cout << endl;
					if ( !executeCommand(command) ) 
						return;
					cout << endl;
					}
				} 
			else 
				{
				return;
				}
			}
		catch(StException &stException) 
			{
			switch (stException.getStExceptionType()) 
				{
				case StException::FATAL: return;
				case StException::ERROR: return;
				case StException::CONTROLC: return;   /* not implemented yet */
				case StException::FORMATERROR: break; /* do nothing and continue the loop */
				}
			}
		}
} 

void CompInterface::readCommandFromCommunicator(string command) {

	if ( userInterface->readCommand(command) ) 
		{
		cout << "here" << endl;
		if ( !executeCommand(command) ) 
			return;
		} 
}

bool CompInterface::executeCommand(string command) {

	istringstream buf(command);
	string word;
	
	if (MyString::isBlank(command))
		return true;
	while ( isspace(command.at(0)) )
		command.erase(0, 1);
		
	buf >> word;
	MyString::lowerCase(word);
	int key = MyString::partialFind (commands, 0, word);

	if (key == -1) 
		{
		Stmessage::warning("Unknown command \"" + word + "\"");
		return true;
		}

	if (MyString::partialFind (commands, key+1, word) != -1) 
		{
		Stmessage::warning("Command is ambiguous");
		return true;
		}

	command.replace (0, word.size(), commands[key]);
	word = commands[key];

	if (word == "quit") 
		return false;

	string fakestring = "   begin structurama;\n      " + command + ";\n   end;";
	istringstream fakestream(fakestring);
	nexusFile->interpretCmd(fakestream, userInterface);
	return true;
}










