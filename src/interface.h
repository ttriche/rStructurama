#ifndef INTERFACE_H
#define INTERFACE_H

#include <iostream>
#include <string>
#include <vector>

using namespace std;



class UserInterface {

	public:
                  virtual   ~UserInterface(void) { }
             virtual void   prompt(string p) { }
             virtual bool   readCommand(string &c)=0;
             virtual void   output(string p)=0;
};



class CLIUserInterface : public UserInterface {

	public:
                            CLIUserInterface(void);
                            ~CLIUserInterface(void);
                     void   prompt(string p) { cout << p; }
                     bool   readCommand(string &c) { getline(cin, c); return true; }
                     void   output(string p) { cout << p << endl; }
};



class GUIUserInterface : public UserInterface {

};



class RemoteUserInterface : public UserInterface {

};



class NexusFile;
class CompInterface {

	public:
                            CompInterface(UserInterface *ui, NexusFile *nf);
                     void   readLoop(void);
					 void   readCommandFromCommunicator(string command);
                     bool   executeCommand(string command);
             static const   vector<string> commands; 

	private:
             static const   string commandNames[];
               static int   numCommands(void);
                NexusFile   *nexusFile;
            UserInterface   *userInterface;
};

#endif
