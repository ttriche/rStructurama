#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
#include <string>
#include <vector>
#include <cstdarg>
#include <cstdio>

#undef	MAC_GUI
#undef THREADED_MEAN_PART

#if defined(MAC_GUI)
#include "AppController.h"
extern AppController *appControllerPtr;
#endif

using namespace std;



class StException {

	public:
		enum StExceptionType { FATAL, ERROR, CONTROLC, FORMATERROR }; 
		StException(StExceptionType m) : stExceptionType(m) {}
		StExceptionType getStExceptionType() const { return stExceptionType; }
	private:
		StExceptionType stExceptionType;
  
};

class Stmessage {

	public:
		static void error(string s) 
			{
#			if defined(MAC_GUI)
			NSString *theMsg = [[NSString alloc] initWithCString:s.c_str() encoding:NSASCIIStringEncoding];
			int theChoice = NSRunAlertPanel(@"Error", theMsg, @"OK", nil, nil);
			[theMsg release];
#			else
			cerr << "   Error: " + s + "." << endl;
#			endif
			}
		static void warning(string s)
			{
#			if defined(MAC_GUI)
			NSString *theMsg = [[NSString alloc] initWithCString:s.c_str() encoding:NSASCIIStringEncoding];
			int theChoice = NSRunAlertPanel(@"Warning", theMsg, @"OK", nil, nil);
			[theMsg release];
#			else
			cerr << "   Warning: " + s + "." << endl;
#			endif
			}
	private:
 
};

static void stprint(const char *format, ...) {

	char str[2048];
	va_list ptr;
	va_start (ptr, format);
	vsprintf (str, format, ptr);
#	if defined(MAC_GUI)
	[appControllerPtr printStatus:str];
#	else
	cout << str;
#	endif
}

class MyString {

	public:
		static string lowerCase(string &s) 
			{
			for(string::iterator p=s.begin(); p!=s.end(); p++)
				*p = tolower(*p);
			return s;
			}
			
		static bool isNumber(string &s)
			{
			bool allNumbers = true;
			for(string::iterator p=s.begin(); p!=s.end(); p++)
				if ( !isdigit(*p) && (*p) != '.' )
					allNumbers = false;
			return allNumbers;
			}

		static string trim(const string &str) 
			{
			int i = str.find_first_not_of(" \t\n");
			if (i == -1)
				return "";
			int j = str.find_last_not_of(" \t\n");
			return str.substr(i,j-i+1);
			}
			
		static string readName(istream &c) 
			{
			string name="";
			int ch;
			while ( isspace(ch=c.get()) && ch != EOF )
				{}
			if (ch == EOF)
				return "";
			if (ch == '"')
				{
				while((ch=c.get()) != '"')
					{
					if (ch == EOF)
						Stmessage::error("Unterminated string.");
					else 
						name += ch;
					}
				}
			else 
				{
				c.putback(ch);
				c >> name;
				}
			return name;
			}
			
		static void printName(ostream &c, string name) 
			{ 
			if(name.find_first_of(" \t\n") == name.npos)
				c << name; 
			else
				c << '"' << name << '"';
			}
			
		static bool isBlank(const string str) { return trim(str) == ""; }
		
		static int partialFind(vector<string> words, int offset, string partialWord) 
			{
			vector<string>::iterator p;
			int i=0;
			for (vector<string>::iterator p=words.begin()+offset;p!=words.end();p++,i++) 
				{
				if (p->find(partialWord) == 0) 
					{
					return i;
					}
				}
			return -1;
			}

		static string formatTime(float n) 
			{
			float t = n;
			char temp[10];
			string word = "";
			int nDays = (int)(t / 86400);
			t = (int)t % 86400;
			int nHours = (int)(t / 3600);
			t = (int)t % 3600;
			int nMinutes = (int)(t / 60);
			int nSeconds = (int)t % 60;
			if (nDays > 0)
				{
				sprintf (temp, "%dd:", nDays);
				word += temp;
				}
			if (nDays == 0 && nHours == 0)
				{}
			else
				{
				if (nHours < 10)
					sprintf (temp, "0%dh:", nHours);
				else
					sprintf (temp, "%dh:", nHours);
				word += temp;
				}
			if (nDays == 0 && nHours == 0 && nMinutes == 0)
				{}
			else
				{
				if (nMinutes < 10)
					sprintf (temp, "0%dm:", nMinutes);
				else
					sprintf (temp, "%dm:", nMinutes);
				word += temp;
				}
			if (nSeconds < 10)
				sprintf (temp, "0%ds", nSeconds);
			else
				sprintf (temp, "%ds", nSeconds);
			word += temp;
			return word;
			}

		static bool queryUser(string msg) 
			{
#			if defined(MAC_GUI)
			NSString *theMsg = [[NSString alloc] initWithCString:msg.c_str() encoding:NSASCIIStringEncoding];
			int theChoice = NSRunAlertPanel(@"Alert", theMsg, @"Yes", @"No", nil);
			[theMsg release];
			return (bool)theChoice;
#			else
			for(;;)
				{
				cout << msg << " (yes/no): ";
				string reply = "";
				cin >> reply;
				reply = lowerCase(reply);
				if (reply == "yes" || reply == "ye" || reply == "y")
					return true;
				else if (reply == "no" || reply == "n")
					return false;
				else
					cout << "Invalid response" << endl;
				}
#			endif
			}
};

#endif
