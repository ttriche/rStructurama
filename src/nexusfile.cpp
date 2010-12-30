#include "interface.h"
#include "iomanager.h"
#include "mcmc.h"
#include "nexusfile.h"
#include "observation.h"
#include "samples.h"
#include "st.h"
#include <iomanip>
#include <sstream>
#include <istream>
#include <fstream>
#include <map>
#include <utility>

#define MISSING 1010101



const string NexusStructuramaBlock::parmNames[] = { "ngen", "nchains", "temperature", "printfreq", "samplefreq", "outfile", "filename", "burnin", "savetofile", "meanpartition", "numpops", "concparmprior", "admixconcparmprior", "admixture", "admixtureprior", "warnuser" };
const vector<string> NexusStructuramaBlock::parms(parmNames,parmNames+numParms());
int NexusStructuramaBlock::numParms() { return sizeof(parmNames)/sizeof(string); }


NexusFile::NexusFile(void) {
	
	setAdmixture( false );
	setAdmixtureDist( "uniform" );
	setAdmixtureParm1( 0.0 );
	setAdmixtureParm2( 10.0 );
	setAdmixtureConcentrationDist( "fixed_expk" );
	setAdmixtureConcentrationParm1( 2.0 );
	setAdmixtureConcentrationParm2( -1.0 );
	setBurnInMarginalLike( 0 );
	setBurnInMeanPart( 0 );
	setBurnInTogetherness( 0 );
	setBurnInNumPops( 0 );
	setConcentrationDist( "fixed_expk" );
	setConcentrationParm1( 3.0 );
	setConcentrationParm2( -1.0 );
	setNgen( 1000000 );
	setNumChains( 1 );
	setNumPops( 2 );
	setPrintfreq( 100 );
	setSamplefreq( 25 );
	setSampleRead( false );
	setSaveToFileMeanPart( false );
	setSaveToFileNumPops( false );
	setSaveToFileTogetherness( false );
	setSumMeanpart( true );
	setTemperature( 0.2 );
	setWarnUser( true );
	
	sample = NULL;
	
	chainOut        = new IoManager();
	meanPartOut     = new IoManager();
	numPopsOut      = new IoManager();
	togethernessOut = new IoManager();
	mcmcSamples     = new McmcSamples();

	chainOut->setFileName( "strout.p" );
	meanPartOut->setFileName( "strout.mp" );
	numPopsOut->setFileName( "strout.np" );
	togethernessOut->setFileName( "strout.tg" );
}

NexusFile::~NexusFile(void) {

	delete chainOut;
	delete meanPartOut;
	delete numPopsOut;
	delete togethernessOut;
	for (vector<NexusBlock *>::iterator p=nexusBlocks.begin(); p != nexusBlocks.end(); p++)
		delete *p;
	delete mcmcSamples;
}

void NexusFile::interpretCmd(istream &c, UserInterface *ui) {

	string line;
	readNexusID(c);
	string command;
	while (!MyString::isBlank( command=NexusCommand(c) )) 
		{
		istringstream str(command);
		string word;
		str >> word;
		MyString::lowerCase(word);
		if (word != "begin") 
			{
			cerr << "word = " << word << endl;
			Stmessage::error("Ill-formated Nexus block. No begin found");
			throw StException(StException::FORMATERROR);
			}
		str >> word;
		MyString::lowerCase(word);
		if (word == "data") 
			{
			nexusBlocks.push_back( new NexusDataBlock(*this, c, ui) );
			}
		else if(word=="structurama") 
			{
			nexusBlocks.push_back( new NexusStructuramaBlock(*this, c, ui) );
			}
		line = "";
		}
}

bool NexusFile::makeNewSample() {

	sample = new Observations(getNumIndividuals(), getNumLoci());
	return true;
}

bool NexusFile::deleteSample() {

	delete sample;
	sample = NULL;
	return true;
}

void NexusFile::readNexusID (istream &c) {

	int ch;
	while (isspace(ch=c.get()) && ch!=EOF)
		{}
	if (ch == '#') 
		{
		string word;
		c >> word;
		if (word != "NEXUS")
			Stmessage::warning("Malformed file identifier " + word);
		}
	else if (ch != EOF)
		c.putback(ch);
}

NexusCommand::NexusCommand(istream &c) : string("") {

	int commentCount = 0;
	int ch;
	
	while ( (ch=c.get()) != EOF ) 
		{
		if (ch == ';' && commentCount == 0)
			return;
		else if (ch == '[')
			commentCount++;
		else if (ch == ']') 
			{
			if (--commentCount < 0)
				Stmessage::error("Unmatched right bracket");
			}
		else if (commentCount == 0)
		*this += ch;
		}
	if (commentCount != 0)
		Stmessage::warning("Unterminated comment");
	if (!MyString::isBlank(*this))
		Stmessage::warning("Ignoring extraneous characters at end of fileSt");
	return;
}

NexusDataBlock::NexusDataBlock(NexusFile &parent, istream &c, UserInterface *ui) : dimensionsFound(false) {

	string command;
	while( (command = NexusCommand(c)) != "" ) 
		{
		istringstream str(command);
		string word;
		str >> word;
		MyString::lowerCase(word);
		if (word == "dimensions")
			{
			if (parent.getSampleRead() == true)
				{
				stprint("      Deleting previously read observations\n");
				parent.deleteSample();
				parent.setSampleRead(false);
				}
			doDimensions(parent, str, ui);
			}
		else if (word == "info") 	
			{
			if (doInfo(parent, str, ui) == false)
				Stmessage::error("Problem reading data");
			else
				parent.setSampleRead(true);
			}
		else if (word=="end")
			return;
		else 
			Stmessage::warning("Command not understood " + word);
		}
}

void NexusDataBlock::doDimensions(NexusFile &parent, istream &c, UserInterface *ui) {

	int ni=0;
	int nl=0;
	
	if (dimensionsFound) 
		{
		Stmessage::warning("Only one dimensions statement allowed per block. Statement ignored");
		}
	else 
		{
		string key, val;
		while (readSetting(c, key, val)) 
			{
			if (val == "")
				{
				Stmessage::error("No parameter given to parameter \"" + key + "\".");
				return;
				}
			if (key == "nind") 
				{
				istringstream buf(val);
				int v;
				buf >> v;
				ni = v;
				}
			else if (key == "nloci") 
				{
				istringstream buf(val);
				int v;
				buf >> v;
				nl = v;
				}
			else
				{
				Stmessage::warning("Could not understand dimensions statement");
				}
			}
		}
		
	if (ni > 0 && nl > 0)
		{
		dimensionsFound = true;
		parent.setNumIndividuals(ni);
		parent.setNumLoci(nl);
		parent.makeNewSample();
		stprint("      Number of individuals               = %d\n", parent.getNumIndividuals());
		stprint("      Number of loci                      = %d\n", parent.getNumLoci());
		if (parent.getSampleRead())
			stprint("      Memory allocated for sample\n");
		}
}

bool NexusDataBlock::doInfo(NexusFile &parent, istream &c, UserInterface *ui) {

	bool readingLabel     = true;
	bool readingGps       = false;
	bool readingLocus     = false;
	bool labelRead        = false;
	int alleleNum         = 0;
	int locusNum          = 0;
	int indNum            = 0;
	int coordinateNum     = 0;
	string indName        = "";
	string allele         = "";
	string gps            = "";
	int a[2]              = {0,0};
	double coordinates[2] = {0.0,0.0};

	char ch;
	while ( (ch = c.get()) != EOF)
		{
		if ( ch == '(' )
			{
			readingLocus = true;
			readingGps   = false;
			readingLabel = false;
			}
		else if ( ch == '{' )
			{
			readingLocus = false;
			readingGps   = true;
			readingLabel = false;
			}
		else if ( ch == '}' )
			{
			double v;
			istringstream buf(gps);
			buf >> v;
			coordinates[coordinateNum++] = v;
			gps           = "";
			coordinateNum = 0;
			readingGps    = false;
			readingLabel  = false;
			readingLocus  = false;
			Individual *indPtr = parent.samplePtr()->getPtrToIndividual(indNum);
			indPtr->setLatitude(coordinates[0]);
			indPtr->setLongitude(coordinates[1]);
			}
		else if ( ch == ')' )
			{
			int v;
			if ( allele == "?" )
				v = MISSING;
			else
				{
				istringstream buf(allele);
				buf >> v;
				}
			if ( alleleNum == 0 || alleleNum == 1 )
				a[alleleNum] = v;
			else
				{
				Stmessage::error("Found too many or too few alleles at locus");
				return false;
				}
			//stprint(locusNum << " -- " << a[0] << " " << a[1]\n");
			Individual *indPtr = parent.samplePtr()->getPtrToIndividual(indNum);
			indPtr->setRawInfo(locusNum, a, alleleNum+1);
			
			allele = "";
			readingLocus = false;
			locusNum++;
			alleleNum = 0;
			}
		else if ( ch == ',' )
			{
			if ( readingLocus == false && readingGps == false )
				{
				// end of an individual
				readingLabel = true;
				labelRead = false;
				indNum++;
				locusNum = alleleNum = 0;
				}
			else if ( readingGps == true )
				{
				double v;
				istringstream buf(gps);
				buf >> v;
				coordinates[coordinateNum++] = v;
				gps = "";
				}
			else
				{
				int v;
				if ( allele == "?" )
					v = MISSING;
				else
					{
					istringstream buf(allele);
					buf >> v;
					}
				if ( alleleNum == 0 || alleleNum == 1 )
					a[alleleNum] = v;
				else
					{
					Stmessage::error("Found too many or too few alleles at locus");
					return false;
					}
				allele = "";
				alleleNum++;
				}
			}
		else if ( ch == ' ' || ch == '\t' || ch == '\n' )
			{
			if ( readingLabel == true && labelRead == true )
				{
				readingLabel = false;
				Individual *indPtr = parent.samplePtr()->getPtrToIndividual(indNum);
				indPtr->setName(indName);
				//stprint(indNum << " -> \"" << indPtr->getName() << "\""\n");
				indName = "";
				}
			}
		else
			{
			if ( readingLabel == true )
				{
				labelRead = true;
				indName += ch;
				}
			else if ( readingLocus == true )
				{
				allele += ch;
				}
			else if ( readingGps == true )
				{
				gps += ch;
cout << "gps = " << gps << endl;
				}
			}
			
		}

	/* recode the information so that alleles are labelled 0 to n-1, where n is the number
	   of unique alleles at the locus */
	parent.samplePtr()->recodeInfo();

	return true;
}

void NexusBlock::readKeyValPairs(istream &c, map<string, string> &keyValPairs) {

	/* break the stream up into tokens, and put each token into a string vector */
	vector<string> keyVals;
	string word;
	while ( c >> word )
		{
		size_t pos=word.find("=");
		if ( pos != word.npos )
			{
			string w1 = word.substr(0,pos);
			string w2 = word.substr(pos, 1);
			string w3 = word.substr(pos+1);
			if (w1 != "")
				keyVals.push_back( w1 );
			if (w2 != "")
				keyVals.push_back( w2 );
			if (w3 != "")
				keyVals.push_back( w3 );
			}
		else
			keyVals.push_back( word );
		}
		
	/* place the string tokens into a map of key/value pairs */
	vector<int> equalSignPositions;
	int i = 0;
	for (vector<string>::iterator p=keyVals.begin(); p != keyVals.end(); p++)
		{
		if ( (*p) == "=" )
			equalSignPositions.push_back( i );
		i++;
		}
		
	for (vector<int>::iterator p=equalSignPositions.begin(); p != equalSignPositions.end(); p++)
		{
		int pos = (*p);
		int numValsAfterEqual = 0;
		if ( p + 1 == equalSignPositions.end() )
			numValsAfterEqual = keyVals.size() - pos - 1;
		else
			numValsAfterEqual = *(p+1) - pos - 2;
		string key = keyVals[pos-1];
		string val = "";
		for (i=1; i<=numValsAfterEqual; i++)
			{
			if (i == 1)
				val = keyVals[pos+1];
			else
				val += (" " + keyVals[pos+i]);
			}
		MyString::lowerCase(key);
		pair<string, string> tempPair;
		tempPair.first = key;
		tempPair.second = val;
		keyValPairs.insert( tempPair );
		}
		
	
	/*for (map<string, string>::iterator p=keyValPairs.begin(); p != keyValPairs.end(); p++)
		stprint("\"" << p->first << "\"" << " -- " << "\"" << p->second << "\""\n");
	getchar();*/
}

bool NexusBlock::readSetting(istream &c, string &key, string &val) {
	
	string word;
		
	if( c >> word ) 
		{
		MyString::lowerCase(word);
		size_t pos=word.find("=");
		if ( pos == word.npos ) 
			{
			char ch;
			c >> ch;
			if(ch != '=') 
				{
				c.putback(ch);
				key = word;
				val = "";
				return true;
				}
			key = word;
			c >> val;
			MyString::lowerCase(val);
			return true;
			}
		else 
			{
			key = word.substr(0,pos);
			val = word.substr(pos+1);
			if (MyString::isBlank(val))
				c >> val;
			MyString::lowerCase(val);
			return true;
			}
		}
	return false;
}

NexusStructuramaBlock::NexusStructuramaBlock(NexusFile &parent, istream &c, UserInterface *ui) {

	string command;

	while ( (command = NexusCommand(c)) != "" ) 
		{
		istringstream str(command);
		string word;
		str >> word;
		MyString::lowerCase(word);
		
		int key = MyString::partialFind (parent.getCompInterfacePtr()->commands, 0, word);
		if (key == -1 && word != "end") 
			{
			Stmessage::warning("Unknown command \"" + word + "\"");
			return;
			}
		if (MyString::partialFind (parent.getCompInterfacePtr()->commands, key+1, word) != -1) 
			{
			Stmessage::warning("Command is ambiguous");
			return;
			}
		if (word != "end")
			word = parent.getCompInterfacePtr()->commands[key];
		
		if (word == "mcmc") 
			doMcmc(parent, str, ui);
		else if (word == "model") 
			doModel(parent, str, ui);
		else if (word == "execute") 
			doExecute(parent, str, ui);
		else if (word == "showdata")
			doShowdata(parent, str, ui);
		else if (word == "help")
			doHelp(parent, str, ui);
		else if (word == "acknowledgments")
			doAcknowledge(parent, str, ui);
		else if (word == "citations")
			doCitation(parent, str, ui);
		else if (word == "readsamples")
			doReadSamples(parent, str, ui);
		else if (word == "set")
			doSet(parent, str, ui);
		else if (word == "showmeanpart")
			doShowMeanPart(parent, str, ui);
		else if (word == "showtogetherness")
			doShowTogetherness(parent, str, ui);
		else if (word == "shownumpops")
			doShowNumPops(parent, str, ui);
		else if (word == "calcmarginal")
			doCalcMarginal(parent, str, ui);
		else if (word == "end") 
			return;
		else 
			Stmessage::warning("Unknown command "+word);
		}
}

void NexusStructuramaBlock::doAcknowledge(NexusFile &parent, istream &str, UserInterface *ui) {

	string temp;
	str >> temp;
	if (temp != "")
		{
		Stmessage::warning("Too many parameters given to Acknowlegments");
		return;
		}

	stprint("   Acknowledgments:\n\n");
	stprint("   The development of this program, and the theory behind it, were supported\n");
	stprint("   by NSF and NIH grants DEB-0445453 and GM-069801, respectively, awarded to\n");
	stprint("   J.P.H. The authors would also like to thank Youn-Ho Lee and the Kordi South\n");
	stprint("   Sea Institute for their support during the development of this program.\n");
#	if defined(MAC_GUI)
	stprint("\n");
#	endif
}

void NexusStructuramaBlock::doCalcMarginal(NexusFile &parent, istream &str, UserInterface *ui) {

	// read the possible parameters for this command
	int bi = parent.getBurnInMarginalLike();
	map<string, string> pairs;
	readKeyValPairs(str, pairs);
	for (map<string, string>::iterator p=pairs.begin(); p != pairs.end(); p++)
		{
		string key = p->first;
		string val = p->second;
		int i = MyString::partialFind (parms, 0, key);
		if (i == -1)
			{
			Stmessage::error("Parameter " + key + " not understood");
			return;
			}
		else
			{
			string parm = parms[i];
			if (parm == "burnin")
				{
				istringstream buf(val);
				int v;
				buf >> v;
				bi = v;
				}
			else
				{
				Stmessage::warning("Parameter " + key + " not understood");
				return;
				}
			}
		}
		
	// check that we have some samples
	if (parent.getMcmcSamplesPtr()->isEmpty() == true)
		{
		Stmessage::warning("No MCMC samples have been read into computer memory");
		return;
		}

	// set the parameters, in case they have changed
	parent.setBurnInMarginalLike( bi );
	
	// calculate the mean partition
	parent.getMcmcSamplesPtr()->calcMarginalLike( parent );
	
#	if defined(MAC_GUI)
	stprint("\n");
#	endif
}

void NexusStructuramaBlock::doCitation(NexusFile &parent, istream &str, UserInterface *ui) {

	string temp;
	str >> temp;
	if (temp != "")
		{
		Stmessage::warning("Too many parameters given to Citations");
		return;
		}

	stprint("   You should cite this program if you use it in a publication. The\n");
	stprint("   appropriate citation is:\n\n");
	stprint("      Huelsenbeck, J. P., E. T. Huelsenbeck, and P. Andolfatto.\n");
	stprint("         In Prep. Structurama: Bayesian inference of population\n");
	stprint("         structure. Bioinformatics.\n\n");
	stprint("   The program uses theory developed in the following papers:\n\n");
	stprint("      Pella, J., and M. Masuda. 2006. The Gibbs and split-merge\n");
	stprint("         sampler for population mixture analysis from genetic data\n");
	stprint("         with incomplete baselines. Can. J. Fish. Aquat. Sci.\n");
	stprint("         63:576-596.\n\n");
	stprint("      Pritchard, J. K., M. Stephens, and P. Donnelly. 2000. Infer-\n");
	stprint("         ence of population structure using multilocus genotype data.\n");
	stprint("         Genetics, 155:945-959.\n\n");
	stprint("      Huelsenbeck, J. P., and P. Andolfatto. 2007. Inference of\n");
	stprint("         population structure under a Dirichlet process model.\n");
	stprint("         Genetics, 175:1787-1802.\n\n");
	stprint("   which you should also cite.\n");
#	if defined(MAC_GUI)
	stprint("\n");
#	endif
}

void NexusStructuramaBlock::doExecute(NexusFile &parent, istream &str, UserInterface *ui) {

	/* get the name of the file (and perhaps the path too) */
	string fileName = "";
	int ch;
	bool leadingWhiteSpace = true;
	while ( (ch=str.get()) != EOF )
		{
		if (ch == ' ' && leadingWhiteSpace == true)
			;
		else
			{
			fileName += ch;
			leadingWhiteSpace = false;
			}
		}

	/* check that there is a file named to be opened */
	if(fileName == "") 
		{
		Stmessage::error("Missing file name");
		return;
		}

	/* open up the file manager and test for the presence
	   of the path and file */
	IoManager in(fileName);
	if ( in.testDirectory() == false )
		{
		Stmessage::error("Cannot find directory \"" + in.getFilePath() + "\"");
		return;
		}
	if ( in.testFile() == false )
		{
		Stmessage::error("Cannot open file \"" + in.getFileName() + "\"");
		return;
		}
	
	/* set output file information */
	parent.getChainOutPtr()->setFilePath( in.getFilePath() );
				
	/* open the file */
	ifstream nexusStream;
	if ( in.openFile( nexusStream ) == false )
		{
		Stmessage::error("Cannot open file \"" + in.getFileName() + "\"");
		return;
		}
	
	/* read the contents of the file */
	stprint("   Reading file \"%s\"\n", in.getFileName().c_str());
	parent.interpretCmd(nexusStream, ui);
	in.closeFile(nexusStream);
	stprint("   Successfully read in file \"%s\"\n", in.getFileName().c_str());
	
	/* set the output file name for the MCMC */
	string mcmcOutName = in.getFileName() + ".p";
	parent.getChainOutPtr()->setFileName(mcmcOutName);

#	if defined(MAC_GUI)
	stprint("\n");
#	endif
}

void NexusStructuramaBlock::doHelp(NexusFile &parent, istream &str, UserInterface *ui) {

	int pos;
	string cmdhelp;
	string temp;
	string helpType;
	string tempCmd;
	
	str >> cmdhelp;
	
	if (cmdhelp == "")
		{
		helpType = "program";
		}
	else
		{
		pos = MyString::partialFind(CompInterface::commands, 0, cmdhelp);
		if (pos == -1)
			{
			tempCmd = MyString::lowerCase(cmdhelp);
			Stmessage::warning("Structurama cannot help you with the unknown command \"" + cmdhelp + "\"");
			return;
			}
		else 
			cmdhelp = CompInterface::commands[pos];
		str >> temp;
		if (temp != "")
			{
			Stmessage::warning("Too many parameters given to Help");
			return;
			}
		helpType = cmdhelp;
		}
	
	if (helpType == "program")
		{
		stprint("   Structurama help                                                       \n\n");   
		stprint("      Commands that are available include:                                \n\n");
		stprint("      Acknowlegments   -- Shows the author's acknowledgments              \n");
		stprint("      Calcmarginal     -- Calculates the marginal likelihood              \n");
		stprint("      Citations        -- Provides suggested citations for the program    \n");
		stprint("      Execute          -- Executes a file                                 \n");
		stprint("      Help             -- Provides a detailed description of the commands \n");
		stprint("      Mcmc             -- Starts the Markov chain Monte Carlo analysis    \n");
		stprint("      Model            -- Sets the assumptions of the analysis            \n");
		stprint("      Readsamples      -- Reads MCMC samples into memory from a file      \n");
		stprint("      Set              -- Set details of program behavior                 \n");
		stprint("      Showdata         -- Shows the data currently in computer memory     \n");
		stprint("      Showmeanpart     -- Calculates the mean partition                   \n");
		stprint("      Shownumpops      -- Summarizes information on number of populations \n");
		stprint("      Showtogetherness -- Probability that pairs are from same population \n");
		stprint("      Quit             -- Quits the program                               \n\n");
		stprint("      Note that Structurama supports the use of the shortest unambiguous  \n");
		stprint("      spelling of the above commands (e.g., \"e\" instead of \"execute\").\n");
		}
	else if (helpType == "quit")
		{
		stprint("   Help Quit                                                              \n\n");
		stprint("      This command quits the program. Even the most novice user           \n");
		stprint("      should be able to easily learn this command.                        \n");
		}
	else if (helpType == "citations")
		{
		stprint("   Help Citations                                                         \n\n");
		stprint("      This command shows suggested citations for the program.             \n");
		}
	else if (helpType == "acknowledgments")
		{
		stprint("   Help Acknowledgments                                                   \n\n");
		stprint("      This command shows the author's acknowledgments.                    \n");
		}
	else if (helpType == "mcmc")
		{
		stprint("   Help Mcmc                                                              \n\n");
		stprint("      Runs the Markov chain Monte Carlo algorithm for approximating       \n");
		stprint("      the posterior probability distribution of the parameters. The       \n");
		stprint("      usage of the command is:                                            \n\n");
		stprint("         Mcmc <parameter> = <value> <parameter> = <value> ...             \n\n");
		stprint("      Valid parameters are:                                               \n\n");
		stprint("         Ngen        -- The number of cycles of the MCMC chain.           \n");
		stprint("         Nchains     -- The number of Markov chains.                      \n");
		stprint("         Temperature -- The temperature parameter for MCMCMC.             \n");
		stprint("         Printfreq   -- Frequency of output to screen.                    \n");
		stprint("         Samplefreq  -- Frequency with which the chain is sampled.        \n");
		stprint("         Outfile     -- The name of the output file.                      \n\n");
		stprint("      Parameter           Options                   Current Setting       \n");
		stprint("      -------------------------------------------------------------       \n");
		stprint("      Ngen                <number>                  %d                    \n", parent.getNgen() );
		stprint("      Nchains             <number>                  %d                    \n", parent.getNumChains() );
		stprint("      Temperature         <number>                  %1.2lf                \n", parent.getTemperature() );
		stprint("      Printfreq           <number>                  %d                    \n", parent.getPrintfreq() );
		stprint("      Samplefreq          <number>                  %d                    \n", parent.getSamplefreq() );
		stprint("      Outfile             <name>                    %s                    \n", parent.getChainOutPtr()->getFileName().c_str() );
		stprint("      -------------------------------------------------------------       \n");
		}
	else if (helpType == "showdata")
		{
		stprint("   Help Showdata                                                         \n\n");
		stprint("      Shows the data currently read into computer memory. The usage      \n");
		stprint("      of the command is                                                  \n\n");
		stprint("         Showdata                                                        \n\n");
		stprint("      and the command does not take any parameters.                      \n");
		}
	else if (helpType == "readsamples")
		{
		stprint("   Help Readsamples                                                      \n\n");
		stprint("      Reads a file containing MCMC samples into computer memory. The     \n");
		stprint("      samples can then be summarized in a variety of ways. The usage     \n");
		stprint("      is:                                                                \n\n");
		stprint("         Readsamples Filename = <name>                                   \n\n");
		stprint("      Note that a file name must be specified.                           \n\n");
		stprint("      Parameter           Options                   Current Setting      \n");
		stprint("      -------------------------------------------------------------      \n");
		stprint("      Filename             <name>                  No default value      \n");
		stprint("      -------------------------------------------------------------      \n");
		}
	else if (helpType == "set")
		{
		stprint("   Help Set                                                              \n\n");
		stprint("      This command customizes some aspects of the program's behavior.    \n");
		stprint("      The usage of the command is:                                       \n\n");
		stprint("         Set <parameter> = <value> <parameter> = <value>                 \n\n");
		stprint("      And the only valid parameter is:                                   \n\n");
		stprint("         Warnuser    -- Should the user be prompted                      \n\n");
		stprint("      Parameter           Options                   Current Setting      \n");
		stprint("      -------------------------------------------------------------      \n");
		if (parent.getWarnUser() == true)
			stprint("      Warnuser            Yes/No                    Yes              \n");
		else
			stprint("      Warnuser            Yes/No                    No               \n");
		stprint("      -------------------------------------------------------------      \n");
		}
	else if (helpType == "showmeanpart")
		{
		stprint("   Help Showmeanpart                                                    \n\n");
		stprint("      Summarizes the results of a Markov chain Monte Carlo analysis     \n");
		stprint("      as a mean partition--a partition that minimizes the squared       \n");
		stprint("      distances to all of the sampled partitions. The usage is:         \n\n");
		stprint("         Showmeanpart <parameter> = <value> <parameter> = <value>       \n\n");
		stprint("      Valid parameters are:                                             \n\n");
		stprint("         Burnin     -- The number of initial samples to discard         \n");
		stprint("         Savetofile -- Whether results should be saved to a file        \n");
		stprint("         Outfile    -- The name of the file for output                  \n\n");
		stprint("      Note that the burn-in is the number of samples in the file        \n");
		stprint("      that should be discarded, not the number of MCMC cycles to        \n");
		stprint("      discard. For example, let's say that a MCMC analysis was          \n");
		stprint("      run for one million cycles, and sampled every 100th cycle         \n");
		stprint("      (i.e., Mcmc ngen=1000000 samplefreq=100). This means that the     \n");
		stprint("      output file specified by the MCMC command will contain 10,000     \n");
		stprint("      samples. The burn-in specifies the number of these 10,000         \n");
		stprint("      samples to discard. The mean partition can take a while to        \n");
		stprint("      calculate. A heuristic search is performed to find the            \n");
		stprint("      partition that minimizes the squared distance to the sampled      \n");
		stprint("      partitions.                                                       \n\n");
		stprint("      Parameter           Options                   Current Setting     \n");
		stprint("      -------------------------------------------------------------     \n");
		stprint("      Burnin              <number>                  %d                  \n", parent.getBurnInMeanPart() );
		stprint("      Savetofile          Yes/No                    ");
		if (parent.getSaveToFileMeanPart() == true)
			stprint("Yes\n");
		else	
			stprint("No\n");
		stprint("      Outfile             <name>                    %s                 \n", parent.getMeanPartOutPtr()->getFileName().c_str() );
		stprint("      -------------------------------------------------------------    \n");
		}
	else if (helpType == "calcmarginal")
		{
		stprint("   Help Calcmarginal                                                   \n\n");
		stprint("      Calculates the marginal likelihood of the model using the        \n");
		stprint("      samples from the Markov chain Monte Carlo analysis. The usage    \n");
		stprint("      is:                                                              \n\n");
		stprint("         Calcmarginal <parameter> = <value>                            \n\n");
		stprint("      Valid parameters are:                                            \n\n");
		stprint("         Burnin     -- The number of initial samples to discard        \n\n");
		stprint("      Note that the burn-in is the number of samples in the file       \n");
		stprint("      that should be discarded, not the number of MCMC cycles to       \n");
		stprint("      discard. The marginal likelihood is calculated as the harmonic   \n");
		stprint("      mean of the likelihoods of the MCMC samples, a method first      \n");
		stprint("      described by Newton and Raftery (1994). The harmonic mean        \n");
		stprint("      is a notoriously unstable estimate of the marginal likelihood,   \n");
		stprint("      and should be interpreted with caution.                          \n\n");
		stprint("      Parameter           Options                   Current Setting    \n");
		stprint("      -------------------------------------------------------------    \n");
		stprint("      Burnin              <number>                  %d                 \n", parent.getBurnInMarginalLike() );
		stprint("      -------------------------------------------------------------    \n");
		}
	else if (helpType == "showtogetherness")
		{
		stprint("   Help Showtogetherness                                               \n\n");
		stprint("      Shows the probability that pairs of individuals are found in     \n");
		stprint("      the same population. The usage is:                               \n\n");
		stprint("         Showtogetherness <parameter> = <value> <parameter> = <value>  \n\n");
		stprint("      Valid parameters are:                                            \n\n");
		stprint("         Burnin     -- The number of initial samples to discard        \n");
		stprint("         Savetofile -- Whether results should be saved to a file       \n");
		stprint("         Outfile    -- The name of the file for output                 \n\n");
		stprint("      Note that the burn-in is the number of samples in the file       \n");
		stprint("      that should be discarded, not the number of MCMC cycles to       \n");
		stprint("      discard. For example, let's say that a MCMC analysis was         \n");
		stprint("      run for one million cycles, and sampled every 100th cycle        \n");
		stprint("      (i.e., Mcmc ngen=1000000 samplefreq=100). This means that the    \n");
		stprint("      output file specified by the MCMC command will contain 10,000    \n");
		stprint("      samples. The burn-in specifies the number of these 10,000        \n");
		stprint("      samples to discard.                                              \n\n");
		stprint("      Parameter           Options                   Current Setting    \n");
		stprint("      -------------------------------------------------------------    \n");
		stprint("      Burnin              <number>                  %d                 \n", parent.getBurnInTogetherness() );
		stprint("      Savetofile          Yes/No                    ");
		if (parent.getSaveToFileTogetherness() == true)
			stprint("Yes\n");
		else	
			stprint("No\n");
		stprint("      Outfile             <name>                    %s                 \n", parent.getTogethernessOutPtr()->getFileName().c_str() );
		stprint("      -------------------------------------------------------------    \n");
		}
	else if (helpType == "shownumpops")
		{
		stprint("   Help Shownumpops                                                    \n\n");
		stprint("      Summarizes the probability distribution of the number of         \n");
		stprint("      populations. The usage is:                                       \n\n");
		stprint("         Shownumpops <parameter> = <value> <parameter> = <value>       \n\n");
		stprint("      Valid parameters are:                                            \n\n");
		stprint("         Burnin     -- The number of initial samples to discard        \n");
		stprint("         Savetofile -- Whether results should be saved to a file       \n");
		stprint("         Outfile    -- The name of the file for output                 \n\n");
		stprint("      Note that the burn-in is the number of samples in the file       \n");
		stprint("      that should be discarded, not the number of MCMC cycles to       \n");
		stprint("      discard. For example, let's say that a MCMC analysis was         \n");
		stprint("      run for one million cycles, and sampled every 100th cycle        \n");
		stprint("      (i.e., Mcmc ngen=1000000 samplefreq=100). This means that the    \n");
		stprint("      output file specified by the MCMC command will contain 10,000    \n");
		stprint("      samples. The burn-in specifies the number of these 10,000        \n");
		stprint("      samples to discard.                                              \n\n");
		stprint("      Parameter           Options                   Current Setting    \n");
		stprint("      -------------------------------------------------------------    \n");
		stprint("      Burnin              <number>                  %d                 \n", parent.getBurnInNumPops() );
		stprint("      Savetofile          Yes/No                    ");
		if (parent.getSaveToFileNumPops() == true)
			stprint("Yes\n");
		else	
			stprint("No\n");
		stprint("      Outfile             <name>                    %s                 \n", parent.getNumPopsOutPtr()->getFileName().c_str() );
		stprint("      -------------------------------------------------------------    \n");
		}
	else if (helpType == "execute")
		{
		stprint("   Help Execute                                                        \n\n");
		stprint("      This command reads the contents of a file that contains the      \n");
		stprint("      observations. The correct usage of the command is:               \n\n");
		stprint("         Execute <file_name>                                           \n");
		}
	else if (helpType == "help")
		{
		stprint("   Help on Help                                                        \n\n");
		stprint("      The help command provides information on the correct usage of    \n");
		stprint("      the commands used in the program Structurama. The correct        \n");
		stprint("      usage of the Help command is either                              \n\n");
		stprint("         Help                                                          \n\n");
		stprint("      or                                                               \n\n");
		stprint("         Help <command name>                                           \n\n");
		stprint("      The help information for some of the commands will also show     \n");
		stprint("      the current settings of the command.                             \n");
		}
	else if (helpType == "model")
		{
		stprint("   Help Model                                                          \n\n");
		stprint("      This command sets the assumptions of the analysis. The usage     \n");
		stprint("      of the command is:                                               \n\n");
		stprint("         Model <parameter> = <value> <parameter> = <value>             \n\n");
		stprint("      Valid parameters are:                                            \n\n");
		stprint("      Numpops                                                          \n\n");
		stprint("         The number of populations. If this is set to a number, then   \n");
		stprint("         the program conditions on having Numpops populations. If,     \n");
		stprint("         on the other hand, Numpops is set to be \"Rv\", then the      \n");
		stprint("         number of populations is considered to be a random variable   \n");
		stprint("         with a Dirichlet process prior.                               \n\n");
		stprint("      Concparmprior                                                    \n\n");
		stprint("         This option specifies the prior distribution on the           \n");
		stprint("         concentration parameter of the Dirichlet process prior        \n");
		stprint("         (which is specified when Numpops=rv). The valid options       \n");
		stprint("         include:                                                      \n\n");
		stprint("            Fixed(<value>)                                             \n");
		stprint("            Fixed_expk(<value>)                                        \n");
		stprint("            Gamma(<shape>,<scale>)                                     \n\n");
		stprint("         Fixed sets the concentration parameter to be some specific    \n");
		stprint("         value. Fixed_expectedk, on the other hand, fixes the          \n");
		stprint("         concentration parameter such that the prior mean of the       \n");
		stprint("         number of populations is some value.                          \n\n");
		stprint("      Admixconcparmprior                                               \n\n");
		stprint("         This option specifies the prior distribution on the           \n");
		stprint("         concentration parameter of the Dirichlet process prior        \n");
		stprint("         specifying the number of admixture components for each        \n");
		stprint("         individual. Valid options include:                            \n\n");
		stprint("            Fixed(<value>)                                             \n");
		stprint("            Fixed_expk(<value>)                                        \n");
		stprint("            Gamma(<shape>,<scale>)                                     \n\n");
		stprint("         Fixed sets the concentration parameter to be some specific    \n");
		stprint("         value. Fixed_expectedk, on the other hand, fixes the          \n");
		stprint("         concentration parameter such that the prior mean of the       \n");
		stprint("         number of admixture components is some value.                 \n\n");
		stprint("      Admixture                                                        \n\n");
		stprint("         This option specifies whether the indiduals are potentially   \n");
		stprint("         an admixture from more than one population. The options are   \n");
		stprint("         either \"Yes\" or \"No\".                                     \n\n");
		stprint("      Admixtureprior                                                   \n\n");
		stprint("         The admixture proportions specify, for each individual,       \n");
		stprint("         the probability that an allele comes from some population.    \n");
		stprint("         These admixture probabilities are assumed to have a           \n");
		stprint("         Dirichlet distribution. The parameters of the Dirichlet are   \n");
		stprint("         specified as Dir(V/K, V/K, ..., V/K), where K is the          \n");
		stprint("         number of populations and V is a parameter that controls      \n");
		stprint("         the variance. The Admixtureprior option controls the          \n");
		stprint("         prior probability distribution on V, and can be               \n\n");
		stprint("            Fixed(<value>)                                             \n");
		stprint("            Uniform(<lower>,<upper>)                                   \n");
		stprint("            Gamma(<shape>,<scale>)                                     \n");
		stprint("            Exponential(<rate>)                                        \n\n");
		stprint("      Parameter            Options                 Current Settings    \n");
		stprint("      --------------------------------------------------------------   \n");
		if (parent.getNumPops() == 0)
			stprint("      Numpops              <number>/Rv             Random Variable \n");
		else
			stprint("      Numpops              <number>/Rv             %d              \n", parent.getNumPops() );
		if (parent.getAdmixture() == false)
			stprint("      Admixture            Yes/No                  No              \n");
		else
			stprint("      Admixture            Yes/No                  Yes             \n");

		stprint("      Concparmprior        Fixed/Fixed_expk/Gamma  ");
		if (parent.getConcentrationDist() == "fixed")
			stprint("Fixed(%1.3lf)\n", parent.getConcentrationParm1() );
		else if (parent.getConcentrationDist() == "fixed_expk")
			stprint("Fixed[E(K)=%1.3lf]\n", parent.getConcentrationParm1() );
		else if (parent.getConcentrationDist() == "gamma")
			stprint("Gamma(%1.3lf,%1.3lf)\n", parent.getConcentrationParm1(), parent.getConcentrationParm2() );

		stprint("      Admixconcparmprior   Fixed/Fixed_expk/Gamma  ");
		if (parent.getAdmixtureConcentrationDist() == "fixed")
			stprint("Fixed(%1.3lf)\n", parent.getAdmixtureConcentrationParm1());
		else if (parent.getAdmixtureConcentrationDist() == "fixed_expk")
			stprint("Fixed[E(K)=%1.3lf]\n", parent.getAdmixtureConcentrationParm1());
		else if (parent.getAdmixtureConcentrationDist() == "gamma")
			stprint("Gamma(%1.3lf,%1.3lf)\n", parent.getAdmixtureConcentrationParm1(), parent.getAdmixtureConcentrationParm2() );

		stprint("      Admixtureprior       Fixed/Uniform/                              \n");
		stprint("                           Gamma/Exponential       ");
		if (parent.getAdmixtureDist() == "fixed")
			stprint("Fixed(%1.3lf)\n", parent.getAdmixtureParm1() );
		else if (parent.getAdmixtureDist() == "uniform")
			stprint("Uniform(%1.3lf,%1.3lf)\n", parent.getAdmixtureParm1(), parent.getAdmixtureParm2());
		else if (parent.getAdmixtureDist() == "gamma")
			stprint("Gamma(%1.3lf,%1.3lf)\n", parent.getAdmixtureParm1(), parent.getAdmixtureParm2() );
		else if (parent.getAdmixtureDist() == "exponential")
			stprint("Exp(%1.3lf)\n", parent.getAdmixtureParm1() );
			
		stprint("      --------------------------------------------------------------   \n");
		}
#	if defined(MAC_GUI)
	stprint("\n");
#	endif
}

void NexusStructuramaBlock::doMcmc(NexusFile &parent, istream &str, UserInterface *ui) {

	/* set the delimiter between directory/file names ("/" or ":") */
	string delimiter = "/";

	/* check that we at least have a date file read into memory to work on */
	bool changedOutfile = false;
	if (parent.getSampleRead() == false)
		{
		Stmessage::warning("No observations have been read in");
		return;
		}
		
	/* get the current settings for all of the parameters used by MCMC */
	int ng    = parent.getNgen();
	int pf    = parent.getPrintfreq();
	int sf    = parent.getSamplefreq();
	int nc    = parent.getNumChains();
	double tp = parent.getTemperature();
	string of = parent.getChainOutPtr()->getFileName();
	string op = parent.getChainOutPtr()->getFilePath();
	
	/* get any user-supplied parameters entered with the mcmc command */
	map<string, string> pairs;
	readKeyValPairs(str, pairs);
	for (map<string, string>::iterator p=pairs.begin(); p != pairs.end(); p++)
		{
		string key = p->first;
		string val = p->second;
		string parm;
		int i = MyString::partialFind (parms, 0, key);
		if (i == -1)
			{
			Stmessage::error("Parameter " + key + " not understood");
			return;
			}
		else
			{
			parm = parms[i];
			if (val == "")
				{
				Stmessage::error("No parameter given to parameter \"" + parm + "\"");
				return;
				}
			istringstream buf(val);
			int v;
			buf >> v;
			if (parm == "ngen")
				ng = v;
			else if (parm == "printfreq")
				pf = v;
			else if (parm == "samplefreq")
				sf = v;
			else if (parm == "nchains")
				nc = v;
			else if (parm == "temperature")
				{
				double dv;
				buf >> dv;
				tp = dv;
				}
			else if (parm == "outfile")
				{
				parent.getChainOutPtr()->parsePathFileNames(val);
				if ( parent.getChainOutPtr()->getFilePath() == "" )
					parent.getChainOutPtr()->setFilePath( parent.getChainOutPtr()->getCurDirectory() );
				of = parent.getChainOutPtr()->getFileName();
				op = parent.getChainOutPtr()->getFilePath();
				changedOutfile = true;
				}
			else
				{
				Stmessage::warning("Parameter " + key + " not understood");
				return;
				}
			}
		}
		
	/* check that the parameters are reasonable */
	if (ng <= 0)
		{
		Stmessage::error("Too few MCMC cycles");
		return;
		}
	if (pf <= 0)
		{
		Stmessage::error("Incorrect print frequency");
		return;
		}
	if (sf <= 0)
		{
		Stmessage::error("Incorrect sample frequency");
		return;
		}
	if (nc < 1)
		{
		Stmessage::error("Too few MCMC chains");
		return;
		}
	if (tp < 0.0)
		{
		Stmessage::error("Temperature parameter is negative");
		return;
		}

	/* test for the presence of the directory and the output file */
	if( parent.getChainOutPtr()->getFileName() == "" ) 
		{
		Stmessage::error("Missing file name");
		return;
		}
		
	/* open up the file manager and test for the presence of the path and file */
	if ( parent.getChainOutPtr()->testDirectory() == false )
		{
		Stmessage::error("Cannot find directory \"" + parent.getChainOutPtr()->getFilePath() + "\"");
		return;
		}

	/* reset the parameters: some of them may have changed */
	parent.setNgen(ng);
	parent.setPrintfreq(pf);
	parent.setSamplefreq(sf);
	parent.setNumChains(nc);
	parent.setTemperature(tp);
			
	/* instantiate and run the Markov chain Monte Carlo analysis */
	Mcmc mcmc(parent);
	mcmc.runChain();

#	if defined(MAC_GUI)
	stprint("\n");
#	endif
}

void NexusStructuramaBlock::doModel(NexusFile &parent, istream &str, UserInterface *ui) {

	int np                           = parent.getNumPops();
	string admixturePriorDist        = parent.getAdmixtureDist();
	string concentrationPriorDist    = parent.getConcentrationDist();
	string admConcentrationPriorDist = parent.getAdmixtureConcentrationDist();
	double admixtureParm1            = parent.getAdmixtureParm1();
	double admixtureParm2            = parent.getAdmixtureParm2();
	double concentrationParm1        = parent.getConcentrationParm1();
	double concentrationParm2        = parent.getConcentrationParm2();
	double admConcentrationParm1     = parent.getAdmixtureConcentrationParm1();
	double admConcentrationParm2     = parent.getAdmixtureConcentrationParm2();

	map<string, string> pairs;
	readKeyValPairs(str, pairs);
	bool changed[5] = { false, false, false, false, false };
	for (map<string, string>::iterator p=pairs.begin(); p != pairs.end(); p++)
		{
		string key = p->first;
		string val = p->second;
		int i = MyString::partialFind (parms, 0, key);
		if (i == -1)
			{
			Stmessage::error("Parameter " + key + " not understood");
			return;
			}
		else
			{
			string parm = parms[i];
			if (val == "")
				{
				Stmessage::error("No parameter given to parameter \"" + parm + "\"");
				return;
				}
			istringstream buf(val);
			
			if (parm == "numpops")
				{
				if (val == "rv" || val == "r")
					np = 0;
				else
					{
					int vi;
					buf >> vi;
					np = vi;
					}
				changed[0] = true;
				}
			else if (parm == "admixture")
				{
				if (val == "yes" || val == "y")
					parent.	setAdmixture(true);
				else if (val == "no" || val == "n")
					parent.setAdmixture(false);
				else
					{
					Stmessage::error("Unknown option \"" + val + "\"");
					return;
					}
				changed[1] = true;
				}
			else if (parm == "concparmprior")
				{
				// get the concentration parameter (of the DPP model) prior distribution
				string dist = "";
				double val1 = 0.0, val2 = 0.0;
				int numArgs = interpretDistribution(val, dist, val1, val2);
				if (numArgs == -1)
					return;
				
				// check 
				if ( dist != "fixed" && dist != "fixed_expk" && dist != "exponential" && dist != "gamma" )
					{
					Stmessage::error("The " + dist + " distribution is not a valid option for Concparmprior");
					return;
					}
				if ( (dist == "fixed" || dist == "fixed_expk" || dist == "exponential") && numArgs != 1 )
					{
					Stmessage::error("The " + dist + " distribution requires one argument");
					return;
					}
				else if (dist == "gamma" && numArgs != 2 )
					{
					Stmessage::error("The " + dist + " distribution requires two arguments");
					return;
					}
				
				// set the distribution
				concentrationPriorDist = dist;
				concentrationParm1     = val1;
				concentrationParm2     = val2;
				changed[2] = true;
				}
			else if (parm == "admixconcparmprior")
				{
				// get the concentration parameter (of the HDPP model) prior distribution
				string dist = "";
				double val1 = 0.0, val2 = 0.0;
				int numArgs = interpretDistribution(val, dist, val1, val2);
				if (numArgs == -1)
					return;
				
				// check 
				if ( dist != "fixed" && dist != "fixed_expk" && dist != "exponential" && dist != "gamma" )
					{
					Stmessage::error("The " + dist + " distribution is not a valid option for Concparmprior");
					return;
					}
				if ( (dist == "fixed" || dist == "fixed_expk" || dist == "exponential") && numArgs != 1 )
					{
					Stmessage::error("The " + dist + " distribution requires one argument");
					return;
					}
				else if (dist == "gamma" && numArgs != 2 )
					{
					Stmessage::error("The " + dist + " distribution requires two arguments");
					return;
					}
				
				// set the distribution

				admConcentrationPriorDist = dist;
				admConcentrationParm1     = val1;
				admConcentrationParm2     = val2;
				changed[3] = true;
				}
			else if (parm == "admixtureprior")
				{
				// get the admixture prior distribution
				string dist = "";
				double val1 = 0.0, val2 = 0.0;
				int numArgs = interpretDistribution(val, dist, val1, val2);
				if (numArgs == -1)
					return;
				
				// check 
				if ( dist != "fixed" && dist != "exponential" && dist != "uniform" && dist != "gamma" )
					{
					Stmessage::error("The " + dist + " distribution is not a valid option for Admixtureprior");
					return;
					}
				if ( (dist == "fixed" || dist == "exponential") && numArgs != 1 )
					{
					Stmessage::error("The " + dist + " distribution requires one argument");
					return;
					}
				else if ( (dist == "uniform" || dist == "gamma") && numArgs != 2 )
					{
					Stmessage::error("The " + dist + " distribution requires two arguments");
					return;
					}
				
				// set the distribution
				admixturePriorDist = dist;
				admixtureParm1     = val1;
				admixtureParm2     = val2;
				changed[4] = true;
				}
			else
				{
				Stmessage::warning("Parameter " + key + " not understood");
				return;
				}
			}
		}
		
	parent.setNumPops(np);
	parent.setAdmixtureDist(admixturePriorDist);
	parent.setAdmixtureParm1(admixtureParm1);
	parent.setAdmixtureParm2(admixtureParm2);
	parent.setConcentrationDist(concentrationPriorDist);
	parent.setConcentrationParm1(concentrationParm1);
	parent.setConcentrationParm2(concentrationParm2);
	parent.setAdmixtureConcentrationDist(admConcentrationPriorDist);
	parent.setAdmixtureConcentrationParm1(admConcentrationParm1);
	parent.setAdmixtureConcentrationParm2(admConcentrationParm2);
	
	if (changed[0] == true || changed[1] == true || changed[2] == true || changed[3] == true)
		stprint("   Setting model\n");
	if (changed[0] == true || changed[2] == true)
		{
		stprint("      Number of populations ");
		if (parent.getNumPops() == 0)
			{
			stprint("is a random variable with a \n");
			if (parent.getConcentrationDist() == "gamma")
				stprint("      Dirichlet process prior [alpha ~ Gamma(%1.3lf,%1.3lf)]\n", parent.getConcentrationParm1(), parent.getConcentrationParm2() );
			else if (parent.getConcentrationDist() == "fixed")
				stprint("      Dirichlet process prior (alpha = %1.3lf)\n", parent.getConcentrationParm1() );
			else if (parent.getConcentrationDist() == "fixed_expk")
				stprint("      Dirichlet process prior [E(k) = %1.3lf]\n", parent.getConcentrationParm1() );
			}
		else
			stprint("= %d\n", parent.getNumPops() );
		}
	if (changed[1] == true)
		{
		if (parent.getAdmixture() == true)
			stprint("      Allowing individuals to be admixed\n");
		else
			stprint("      Individuals are not admixed\n");
		}
	if (changed[3] == true)
		{
		stprint("      Number of admixture classes for each individual is a random variable with a\n");
		if (parent.getAdmixtureConcentrationDist() == "gamma")
			stprint("      Dirichlet process prior [alpha ~ Gamma(%1.3lf,%1.3lf)]\n", parent.getAdmixtureConcentrationParm1(), parent.getAdmixtureConcentrationParm2() );
		else if (parent.getAdmixtureConcentrationDist() == "fixed")
			stprint("      Dirichlet process prior (alpha = %1.3lf)\n", parent.getAdmixtureConcentrationParm1() );
		else if (parent.getAdmixtureConcentrationDist() == "fixed_expk")
			stprint("      Dirichlet process prior [E(k) = %1.3lf]\n", parent.getAdmixtureConcentrationParm1() );
		}
	if (changed[4] == true)
		{
		if (parent.getAdmixtureDist() == "exponential")
			stprint("      Prior on variance parameter of admixture proportions is Exponential(%1.3lf)]\n", parent.getAdmixtureParm1() );
		else if (parent.getAdmixtureDist() == "uniform")
			stprint("      Prior on variance parameter of admixture proportions is Uniform(%1.3lf,%1.3lf)]\n", parent.getAdmixtureParm1(), parent.getAdmixtureParm2() );
		else if (parent.getAdmixtureDist() == "gamma")
			stprint("      Prior on variance parameter of admixture proportions is Gamma(%1.3lf,%1.3lf)]\n", parent.getAdmixtureParm1(), parent.getAdmixtureParm2() );
		else if (parent.getAdmixtureDist() == "fixed")
			stprint("      Variance parameter of admixture proportions is %1.3lf)\n", parent.getAdmixtureParm1() );
		}
	stprint("   Successfully set model\n");
#	if defined(MAC_GUI)
	stprint("\n");
#	endif
}

void NexusStructuramaBlock::doReadSamples(NexusFile &parent, istream &str, UserInterface *ui) {

	// read the file name
	string fileName = "", filePath = "";
	IoManager inMngr;
	map<string, string> pairs;
	readKeyValPairs(str, pairs);
	for (map<string, string>::iterator p=pairs.begin(); p != pairs.end(); p++)
		{
		string key = p->first;
		string val = p->second;
		int i = MyString::partialFind (parms, 0, key);
		if (i == -1)
			{
			Stmessage::error("Parameter " + key + " not understood");
			return;
			}
		else
			{
			string parm = parms[i];
			if (val == "")
				{
				Stmessage::error("No parameter given to parameter \"" + parm + "\"");
				return;
				}
			if (parm == "filename")
				{
				inMngr.parsePathFileNames(val);
				if ( inMngr.getFilePath() == "" )
					inMngr.setFilePath( inMngr.getCurDirectory() );
				fileName = inMngr.getFileName();
				filePath = inMngr.getFilePath();
				}
			else
				{
				Stmessage::warning("Parameter " + key + " not understood");
				return;
				}
			}
		}

	// check that there is a file named to be opened
	if( inMngr.getFileName() == "" ) 
		{
		Stmessage::error("No file name was specified");
		return;
		}
	if( inMngr.getFilePath() == "" ) 
		{
		Stmessage::error("Missing file path");
		return;
		}
		
	// open up the file manager and test for the presence of the path and file
	if ( inMngr.testDirectory() == false )
		{
		Stmessage::error("Cannot find directory \"" + inMngr.getFilePath() + "\"");
		return;
		}
	if ( inMngr.testFile() == false )
		{
		Stmessage::error("Cannot find file \"" + inMngr.getFileName() + "\"");
		return;
		}

	// read the file
	if ( parent.getMcmcSamplesPtr()->readFile( inMngr, parent ) == false )
		stprint("   Terminated reading of file \"%s\"\n", fileName.c_str() );
	else
		stprint("   Successfully read file \"%s\"\n", fileName.c_str() );
#	if defined(MAC_GUI)
	stprint("\n");
#	endif
}

void NexusStructuramaBlock::doSet(NexusFile &parent, istream &str, UserInterface *ui) {

	bool wu = parent.getWarnUser();
	
	map<string, string> pairs;
	readKeyValPairs(str, pairs);
	for (map<string, string>::iterator p=pairs.begin(); p != pairs.end(); p++)
		{
		string key = p->first;
		string val = p->second;
		int i = MyString::partialFind (parms, 0, key);
		if (i == -1)
			{
			Stmessage::error("Parameter " + key + " not understood");
			return;
			}
		else
			{
			string parm = parms[i];
			if (val == "")
				{
				Stmessage::error("No parameter given to \"" + parm + "\"");
				return;
				}
			if (parm == "warnuser")
				{

				if (val[0] == 'y')
					wu = true;
				else if (val[0] == 'n')
					wu = false;
				else
					{
					Stmessage::error("Parameter option \"" + val + "\" not understood");
					return;
					}
				}
			else
				{
				Stmessage::warning("Parameter " + key + " not understood");
				return;
				}
			}
		}

	parent.setWarnUser(wu);

	if (parent.getWarnUser() == true)
		stprint("   Prompting user enabled\n");
	else
		stprint("   Prompting user disabled\n");
#	if defined(MAC_GUI)
	stprint("\n");
#	endif
}

void NexusStructuramaBlock::doShowdata(NexusFile &parent, istream &str, UserInterface *ui) {

	if (parent.getSampleRead() == false)
		{
		Stmessage::warning("No observations have been read in");
		return;
		}

	string temp;
	str >> temp;
	if (temp != "")
		{
		Stmessage::warning("Too many parameters given to Showdata");
		return;
		}

	stprint("   Showing current data in memory\n");
	parent.samplePtr()->print();
#	if defined(MAC_GUI)
	stprint("\n");
#	endif
}

void NexusStructuramaBlock::doShowMeanPart(NexusFile &parent, istream &str, UserInterface *ui) {

	// get a pointer to the io manager
	IoManager *iop = parent.getMeanPartOutPtr();
	
	// get the current parameter values
	int bi = parent.getBurnInMeanPart();
	bool sf = parent.getSaveToFileMeanPart();
	
	// read the possible parameters for this command
	map<string, string> pairs;
	readKeyValPairs(str, pairs);
	for (map<string, string>::iterator p=pairs.begin(); p != pairs.end(); p++)
		{
		string key = p->first;
		string val = p->second;
		int i = MyString::partialFind (parms, 0, key);
		if (i == -1)
			{
			Stmessage::error("Parameter " + key + " not understood");
			return;
			}
		else
			{
			string parm = parms[i];
			if (parm == "burnin")
				{
				istringstream buf(val);
				int v;
				buf >> v;
				bi = v;
				}
			else if (parm == "savetofile")
				{
				vector<string> possibilities;
				possibilities.push_back("no");
				possibilities.push_back("yes");
				int pos = MyString::partialFind(possibilities, 0, val);
				if (pos == -1)
					{
					Stmessage::warning("Cannot interpret " + val);
					return;
					}
				string result = possibilities[pos];
				if (result == "yes")
					sf = true;
				else
					sf = false;
				}
			else if (parm == "outfile")
				{
				iop->parsePathFileNames(val);
				if ( iop->getFilePath() == "" )
					iop->setFilePath( iop->getCurDirectory() );
				}
			else
				{
				Stmessage::warning("Parameter " + key + " not understood");
				return;
				}
			}
		}
		
	// check that we have some samples
	if (parent.getMcmcSamplesPtr()->isEmpty() == true)
		{
		Stmessage::warning("No MCMC samples have been read into computer memory");
		return;
		}

	// open up the file manager and test for the presence of the path and file
	if ( iop->testDirectory() == false )
		{
		Stmessage::error("Cannot find directory \"" + iop->getFilePath() + "\"");
		return;
		}
		
	// set the parameters, in case they have changed
	parent.setBurnInMeanPart( bi );
	parent.setSaveToFileMeanPart( sf );
	
	// calculate the mean partition
	parent.getMcmcSamplesPtr()->calcMeanPart( parent );

#	if defined(MAC_GUI)
	stprint("\n");
#	endif
}

void NexusStructuramaBlock::doShowNumPops(NexusFile &parent, istream &str, UserInterface *ui) {

	// get a pointer to the io manager
	IoManager *iop = parent.getNumPopsOutPtr();
	
	// get the current parameter values
	int bi = parent.getBurnInNumPops();
	bool sf = parent.getSaveToFileNumPops();
	
	// read the possible parameters for this command
	map<string, string> pairs;
	readKeyValPairs(str, pairs);
	for (map<string, string>::iterator p=pairs.begin(); p != pairs.end(); p++)
		{
		string key = p->first;
		string val = p->second;
		int i = MyString::partialFind (parms, 0, key);
		if (i == -1)
			{
			Stmessage::error("Parameter " + key + " not understood");
			return;
			}
		else
			{
			string parm = parms[i];
			if (parm == "burnin")
				{
				istringstream buf(val);
				int v;
				buf >> v;
				bi = v;
				}
			else if (parm == "savetofile")
				{
				vector<string> possibilities;
				possibilities.push_back("no");
				possibilities.push_back("yes");
				int pos = MyString::partialFind(possibilities, 0, val);
				if (pos == -1)
					{
					Stmessage::warning("Cannot interpret " + val);
					return;
					}
				string result = possibilities[pos];
				if (result == "yes")
					sf = true;
				else
					sf = false;
				}
			else if (parm == "outfile")
				{
				iop->parsePathFileNames(val);
				if ( iop->getFilePath() == "" )
					iop->setFilePath( iop->getCurDirectory() );
				}
			else
				{
				Stmessage::warning("Parameter " + key + " not understood");
				return;
				}
			}
		}
		
	// check that we have some samples
	if (parent.getMcmcSamplesPtr()->isEmpty() == true)
		{
		Stmessage::warning("No MCMC samples have been read into computer memory");
		return;
		}

	// open up the file manager and test for the presence of the path and file
	if ( iop->testDirectory() == false )
		{
		Stmessage::error("Cannot find directory \"" + iop->getFilePath() + "\"");
		return;
		}
		
	// set the parameters, in case they have changed
	parent.setBurnInNumPops( bi );
	parent.setSaveToFileNumPops( sf );
	
	// calculate the mean partition
	parent.getMcmcSamplesPtr()->calcNumPops( parent );

#	if defined(MAC_GUI)
	stprint("\n");
#	endif
}

void NexusStructuramaBlock::doShowTogetherness(NexusFile &parent, istream &str, UserInterface *ui) {

	// get a pointer to the io manager
	IoManager *iop = parent.getTogethernessOutPtr();
	
	// get the current parameter values
	int bi = parent.getBurnInTogetherness();
	bool sf = parent.getSaveToFileTogetherness();
	
	// read the possible parameters for this command
	map<string, string> pairs;
	readKeyValPairs(str, pairs);
	for (map<string, string>::iterator p=pairs.begin(); p != pairs.end(); p++)
		{
		string key = p->first;
		string val = p->second;
		int i = MyString::partialFind (parms, 0, key);
		if (i == -1)
			{
			Stmessage::error("Parameter " + key + " not understood");
			return;
			}
		else
			{
			string parm = parms[i];
			if (parm == "burnin")
				{
				istringstream buf(val);
				int v;
				buf >> v;
				bi = v;
				}
			else if (parm == "savetofile")
				{
				vector<string> possibilities;
				possibilities.push_back("no");
				possibilities.push_back("yes");
				int pos = MyString::partialFind(possibilities, 0, val);
				if (pos == -1)
					{
					Stmessage::warning("Cannot interpret " + val);
					return;
					}
				string result = possibilities[pos];
				if (result == "yes")
					sf = true;
				else
					sf = false;
				}
			else if (parm == "outfile")
				{
				iop->parsePathFileNames(val);
				if ( iop->getFilePath() == "" )
					iop->setFilePath( iop->getCurDirectory() );
				}
			else
				{
				Stmessage::warning("Parameter " + key + " not understood");
				return;
				}
			}
		}
		
	// check that we have some samples
	if (parent.getMcmcSamplesPtr()->isEmpty() == true)
		{
		Stmessage::warning("No MCMC samples have been read into computer memory");
		return;
		}

	// open up the file manager and test for the presence of the path and file
	if ( iop->testDirectory() == false )
		{
		Stmessage::error("Cannot find directory \"" + iop->getFilePath() + "\"");
		return;
		}
		
	// set the parameters, in case they have changed
	parent.setBurnInTogetherness( bi );
	parent.setSaveToFileTogetherness( sf );
	
	// calculate the mean partition
	parent.getMcmcSamplesPtr()->calcTogetherness( parent );

#	if defined(MAC_GUI)
	stprint("\n");
#	endif
}

int NexusStructuramaBlock::interpretDistribution(string &qs, string &dist, double &val1, double &val2) {

	// parse out the tokens
	string s1 = "", s2 = "", s3 = "";
	int readNum = 0;
	for (int i=0; i<qs.size(); i++)
		{
		if ( qs[i] == '(' || qs[i] == ',' || qs[i] == ')' )
			{
			readNum++;
			continue;
			}
		if ( qs[i] == ' ' )
			continue;
		if (readNum == 0)
			s1 += qs[i];
		else if (readNum == 1)
			s2 += qs[i];
		else if (readNum == 2)
			s3 += qs[i];
		}
		
	// find the distribution
	vector<string> validDistributions;
	validDistributions.push_back( "fixed" );
	validDistributions.push_back( "uniform" );
	validDistributions.push_back( "gamma" );
	validDistributions.push_back( "exponential" );
	validDistributions.push_back( "fixed_expk" );
	MyString::lowerCase(s1);
	int pos = MyString::partialFind(validDistributions, 0, s1);
	if ( pos == -1 )
		{
		Stmessage::error("Unknown argument " + s1);
		return -1;
		}
	dist = validDistributions[pos];
	
	// get the parameters of the distribution
	if ( MyString::isNumber(s2) == false )
		{
		Stmessage::error("First argument of " + dist + " is not a number");
		return -1;
		}
	if ( MyString::isNumber(s3) == false )
		{
		Stmessage::error("First argument of " + dist + " is not a number");
		return -1;
		}
	val1 = -1.0;
	val2 = -1.0;
	int na = 0;
	if ( s2 != "" )
		{
		istringstream buf(s2);
		double v;
		buf >> v;
		val1 = v;
		na++;
		}
	if ( s3 != "" )
		{
		istringstream buf(s3);
		double v;
		buf >> v;
		val2 = v;
		na++;
		}
	return na;
}




