#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <stdio.h>
#include <string>
#include <sys/time.h>
#include <vector>
#include "DataNode.h"
#include "DataSequence.h"
#include "Sequence.h"
#include "SequenceList.h"
#include "SequenceReader.h"
#include "Primer.h"
#include "PrimerPair.h"
#include "PrimerPairList.h"
#include "PrimerCalculator.h"
#include "GlobalSettings.h"

using namespace std;
using namespace DeGenPrime;

int GlobalSettings::_ampLength = DEFAULT_AMPLICON_LENGTH;
int GlobalSettings::_beginningNucleotide = DEFAULT_BEGIN_NUCLEOTIDE;
int GlobalSettings::_endingNucleotide = DEFAULT_END_NUCLEOTIDE;
bool GlobalSettings::_measureByAmpliconSize = DEFAULT_MEASURE_BY_AMPLICON;
float GlobalSettings::_minTemp = DEFAULT_MIN_TEMP;
float GlobalSettings::_maxTemp = DEFAULT_MAX_TEMP;
float GlobalSettings::_primerConcentration = DEFAULT_PRIMER_CONC;
float GlobalSettings::_monovalentIonConcentration = DEFAULT_SALT_CONC;
int GlobalSettings::_maxPrimers = DEFAULT_MAX_PRIMERS;
float GlobalSettings::_thermodynamicTemperature = DEFAULT_THERMODYNAMIC_TEMPERATURE;

void ProcessTags(int argc, char *argv[]);
void PrintHelp();

int main(int argc, char *argv[])
{
	// Set up clock
	struct timeval begin, end;
	gettimeofday(&begin, 0);
	
	// Check if user wants help
	if(argc == 1 || (argc == 2 && (strcmp("--h", argv[1]) == 0 || strcmp("--help", argv[1]) == 0) ) )
	{
		PrintHelp();
	}

	// Process Tags
	if(argc != 2)
	{
		ProcessTags(argc, argv);
	}

	// Create Filename/path
	string filename = argv[argc - 1];
	string filepath = std::filesystem::current_path();

	// Open Output File Stream
	ofstream ofs;
	ofs.open(filepath + "/primers_" + filename.substr(0, filename.length() - 3) + "txt");

	// File Information
	ofs << "Full path of argument: " << filepath << endl;
	ofs << "Full path of file: " << filepath << "/" << filename << endl;

	// Open Input File
	ifstream ifs;
	ifs.open(filepath + "/" + filename);
	if(ifs.fail())
	{
		cout << "Failure to open file.\n";
		exit(BAD_INPUT_FILE);
	}
	ofs << "Input File opened." << endl;

	// Read Sequences
	SequenceReader read;
	SequenceList list = read.CreateList(ifs);

	// Test Sequences
	if(list.size() == 0)
	{
		cout << "There were no sequences in the input file.\n";
		ofs << "There were no sequences in the input file.\n";
		exit(BAD_INPUT_FILE);
	}
	else if(list.TestAlignment() == false)
	{
		cout << "File not aligned.  Running MAFFT, then re-running this program.\n";
		bool local = true;
		if(argc > 2)
		{
			for(int i = 1;i < argc;i++)
			{
				if(local == false)
				{
					break;
				}
				if(argv[i] == "--global" || argv[i] == "--g")
				{
					local = false;
				}
			}
		}
		string clustal = filename.substr(0, filename.find("."));
		clustal += local ? "_local.clust" : "_global.clust";
		string command1 = "mafft ";
		command1 += local ? "--localpair " : "--globalpair ";
		command1 += "--maxiterate 1000 --clustalout --quiet " + filename + " >" + clustal;
		string command2 = "";
		for(int i = 0;i < argc - 1;i++)
		{
			command2 += (string)argv[i] + " ";
		}
		command2 += clustal;
		const char *c1 = command1.c_str();
		const char *c2 = command2.c_str();
		system(c1);
		system(c2);
		exit(FILE_MISALIGNED);
	}
	
	// Sequence Information Before Filtering
	ofs << "Before Filter, the number of Sequences in the list is: " << list.size() << endl;
	ofs << list.PrintSequenceNames();
	list.FilterDashes();
	
	// Sequence Information After Filtering
	ofs << "After Filter, the number of Sequences in the list is: " << list.size() << endl;
	ofs << list.PrintSequenceNames() << endl;

	// Create Forward and Reverse DataSequences
	SequenceList reverse_list = list.InvRevList();
	DataSequence data = list.ProcessList();
	DataSequence rev = reverse_list.ProcessList();

	// Create Primer Calculators
	PrimerCalculator calc, rev_calc;
	if(GlobalSettings::GetMeasureByAmpliconSize()) // User wants minimum amplicon length
	{
		calc.InitializePrimers(data);
		rev_calc.InitializePrimers(rev);
	}
	else // User specified a range
	{
		GlobalSettings::SetMinimumAmplicon(GlobalSettings::GetEndingNucleotide() - GlobalSettings::GetBeginningNucleotide());
		calc.InitializeBoundedPrimers(data, GlobalSettings::GetBeginningNucleotide());
		int rev_lowerBound = data.RevIndex(GlobalSettings::GetEndingNucleotide());
		rev_calc.InitializeBoundedPrimers(rev, rev_lowerBound);
	}

	// Display number of possible primers, run filters and output filter percentages.
	cout << "Number of possible forward primers in specified size range, before Filters: [" << calc.size() << "]\n" << endl;
	ofs << "Number of possible forward primers in specified size range, before Filters: [" << calc.size() << "]\n" << endl;
	ofs << calc.FilterAll(data, list) << endl;
	cout << "After running Filters: [" << calc.size() << "]\n" << endl;

	cout << "Number of possible reverse primers in specified size range, before Filters: [" << rev_calc.size() << "]\n" << endl;
	ofs << "Number of possible reverse primers in specified size range, before Filters: [" << rev_calc.size() << "]\n" << endl;
	ofs << rev_calc.FilterAll(rev, reverse_list) << endl;
	cout << "After running Filters: [" << rev_calc.size() << "]\n" << endl;
	
	// Create PrimerPairList
	PrimerPairList pairlist;
	pairlist.CreateList(data, rev, calc.GetPrimers(), rev_calc.GetPrimers());
	cout << "Before filters, the number of forward-reverse primer pairs in this list is: " << pairlist.size() << endl;
	ofs << "Before filters, the number of forward-reverse primer pairs in this list is: " << pairlist.size() << endl;

	// Filter PrimerPairList
	ofs << pairlist.FilterAmpliconLength();
	ofs << pairlist.FilterTemperatureDifference();
	ofs << pairlist.FilterMessage("final", 0);
	cout << "After filters, the number of forward-reverse primer pairs in this list is: " << pairlist.size() << endl;

	// Sort PrimerPairList by least temperature difference.
	ofs << "\nSorting remaining " << pairlist.size();
	ofs << " primer pairs in the list by temperature difference." << endl;
	pairlist.Sort();

	// Loop through top desired primer pairs to filter them for annealing temperature
	// and output the final list of primer pairs.
	const int desiredpairs = GlobalSettings::GetMaximumReturnPrimers();
	ofs << "Running FilterAnnealingTemp on the top " << desiredpairs << " primer pairs." << endl;
	PrimerPairList top = pairlist.SubList(0,desiredpairs);
	if(top.size() == 0)
	{
		cout << "No primer pairs were found that meet your specifications." << endl;
		ofs << "No primer pairs were found that meet your specifications." << endl;
	}
	else
	{
		int nextIndex = desiredpairs;
		int filtercount = top.FilterAnnealingTemp(data, rev);
		int remaining = pairlist.size() - desiredpairs;
		int loopcount = 1;
		while(filtercount != 0 && remaining > 0)
		{
			int nextlength = (filtercount < remaining) ? filtercount : remaining;
			PrimerPairList next = pairlist.SubList(nextIndex,nextlength);
			filtercount = next.FilterAnnealingTemp(data, rev);
			loopcount++;
			nextIndex += nextlength;
			remaining -= nextlength;
			top.Append(next);
		}
		if(top.size() == 0)
		{
			cout << "No primer pairs were found that meet your specifications." << endl;
			ofs << "No primer pairs were found that meet your specifications." << endl;
		}
		else
		{
			cout << top.PrintAll(data, rev);
			ofs << top.PrintAll(data, rev);
		}
	}

	// Close input/output file streams.
	ifs.close();
	ofs.close();
	
	// Show closing messages and clock to user, then close the program.
	cout << "Output details saved to primers_" << filename.substr(0, filename.length() - 3) << "txt" << endl;
	gettimeofday(&end, 0);
	long seconds = end.tv_sec - begin.tv_sec;
	long microseconds = end.tv_usec - begin.tv_usec;
	double elapsed = seconds + microseconds*1e-6;
	printf("Program running time measured: %.3f seconds.\n", elapsed);
	cout << "Program complete." << endl;
	exit(PROGRAM_SUCCESS);
}

void ProcessTags(int argc, char *argv[])
{
	bool containsAmplicon = false;
	bool containsBegin = false;
	bool containsEnd = false;
	char *ptr;
	int value;
	for(int i = 1;i < argc - 1;i++)
	{
		if(strcmp("--h", argv[i]) == 0 || strcmp("--help", argv[i]) == 0)
		{
			PrintHelp();
		}

		// For most tags, they follow the format: '--tag:value'
		// We want to check what the tag is with strstr()
		// and store the given value to global settings with strchr()
		ptr = strchr(argv[i], ':') + 1;
		if(ptr != NULL)
		{
			value = atoi(ptr);
		}
		
		if(strstr(argv[i], "--amplicon:") != NULL)
		{
			GlobalSettings::SetMinimumAmplicon(value);
			containsAmplicon = true;
		}
		else if(strstr(argv[i], "--begin:") != NULL)
		{
			GlobalSettings::SetBeginningNucleotide(value);
			containsBegin = true;
		}
		else if(strstr(argv[i], "--end:") != NULL)
		{
			GlobalSettings::SetEndingNucleotide(value);
			containsEnd = true;
		}
		else if(strstr(argv[i], "--min_temp:") != NULL)
		{
			GlobalSettings::SetMinimumTemperature(value);
		}
		else if(strstr(argv[i], "--max_temp:") != NULL)
		{
			GlobalSettings::SetMaximumTemperature(value);
		}
		else if(strstr(argv[i], "--primer_conc:") != NULL)
		{
			GlobalSettings::SetPrimerConcentration(value);
		}
		else if(strstr(argv[i], "--salt_conc:") != NULL)
		{
			GlobalSettings::SetMonoIonConcentration(value);
		}
		else if(strstr(argv[i], "--max_primers:") != NULL)
		{
			GlobalSettings::SetMaximumReturnPrimers(value);
		}
		else
		{
			continue;
		}
	}

	// We need to check here if user has entered improper values for their tags.
	// close the program if any improper tags were entered.

	// User cannot enter a starting nucleotide without also entering an ending nucleotide
	if(containsBegin != containsEnd)
	{
		string first = containsBegin ? "--begin" : "--end";
		string second = containsBegin ? "--end" : "--begin";
		cout << "Syntax error.  DeGenPrime cannot have a " << first << " tag without a ";
		cout << second << " tag." << endl;
		exit(SETTINGS_FILE_NOT_FOUND);
	}

	// User cannot specify an amplicon length with beginning or ending nucleotides.
	if(containsAmplicon && (containsBegin || containsEnd))
	{
		cout << "Syntax error.  DeGenPrime cannot have an --amplicon tag ";
		cout << "with a --begin or --end tag." << endl;
		exit(SETTINGS_FILE_NOT_FOUND);
	}

	// User cannot specify an ending nucleotide less than the beginning nucleotide
	if(GlobalSettings::GetEndingNucleotide() <= GlobalSettings::GetBeginningNucleotide())
	{
		cout << "Syntax error.  The beginning nucleotide must be less than ";
		cout << "the ending nucleotide." << endl;
		exit(SETTINGS_FILE_NOT_FOUND);
	}

	// User cannot specify a minimum temperature greater than a maximum temperature
	if(GlobalSettings::GetMinimumTemperature() > GlobalSettings::GetMaximumTemperature())
	{
		cout << "Syntax error.  You cannot have a minimum temperature bigger ";
		cout << "than the maximum temperature." << endl;
		exit(SETTINGS_FILE_NOT_FOUND);
	}

	// Change global settings to measure specific sections if the users specified
	// a region to amplify.
	if(containsBegin)
	{
		GlobalSettings::SetMeasureByAmpliconSize(false);
	}
}

void PrintHelp()
{
	cout << "Syntax: ./DeGenPrime [--tags] <filename>\n";
	cout << "Valid tags include:\n";
	cout << "\t--amplicon:int, Set the minimum amplicon length.  This will not work with";
	cout << "--begin or --end tags.\n";
	cout << "\t--begin:int, Set the beginning nucleotide.  This must be used with --end ";
	cout << "and cannot be used with --amplicon.\n";
	cout << "\t--end:int, Set the ending nucleotide.  This must be used with --begin ";
	cout << "and cannot be used with --amplicon.\n";
	cout << "\t--global or --g, for lists of seqeuences that are misaligned, this tag specifies ";
	cout << "that the file should run mafft for global alignment.\n";
	cout << "\t--help or --h, prints this help menu.\n";
	cout << "\t--local or --l, for lists of sequences that are misaligned, this tag specifies ";
	cout << "that the file should run mafft for local alignment.\n";
	cout << "\t--min_temp:int, Sets the minimum primer melting temperature.  This has";
	cout << " a minimum value of " << MIN_PRIMER_TEMP << " (degrees celsius) and must be ";
	cout << "smaller than --max_temp.\n";
	cout << "\t--max_temp:int, Sets the maximum primer melting temperature.  This has";
	cout << " a maximum value of " << MAX_PRIMER_TEMP << " (degrees celsius) and must be ";
	cout << "larger than --min_temp.\n";
	cout << "\t--primer_conc:int, Sets the concentration of the PCR primer in nM.  This has ";
	cout << "a minimum value of " << MIN_PRIMER_CONC << " and this program will raise ";
	cout << "any value smaller to this value.\n";
	cout << "\t--salt_conc:int, Sets the concentration of monovalent ions in mM.  This has ";
	cout << "a minimum value of " << MIN_SALT_CONC << " and this program will raise ";
	cout << "any value smaller to this value.\n";
	cout << "\t--max_primers:int, Sets the maximum number of outputted primers.  This has ";
	cout << "a maximum value of " << MAX_PRIMER_RETURNS << " and this program will reduce ";
	cout << "any value larger to this value.\n";
	exit(PROGRAM_SUCCESS);
}