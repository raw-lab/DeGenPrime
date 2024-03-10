#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <string>
#include <vector>
#include "config.h"
#include "datanode.h"
#include "datasequence.h"
#include "format.h"
#include "sequence.h"
#include "sequencelist.h"
#include "sequencereader.h"
#include "primer.h"
#include "primerpair.h"
#include "primerpairlist.h"
#include "primercalculator.h"
#include "globalsettings.h"

using namespace std;
using namespace DeGenPrime;

int GlobalSettings::_ampLength = DEFAULT_AMPLICON_LENGTH;
int GlobalSettings::_beginningNucleotide = DEFAULT_BEGIN_NUCLEOTIDE;
float GlobalSettings::_deltag = DEFAULT_DELTA_G;
int GlobalSettings::_endingNucleotide = DEFAULT_END_NUCLEOTIDE;
bool GlobalSettings::_measureByAmpliconSize = DEFAULT_MEASURE_BY_AMPLICON;
bool GlobalSettings::_proteinSequence = DEFAULT_PROTEIN_SEQUENCE;
bool GlobalSettings::_beginflag = DEFAULT_BEGIN_FLAG;
bool GlobalSettings::_endflag = DEFAULT_END_FLAG;
float GlobalSettings::_minTemp = DEFAULT_MIN_TEMP;
float GlobalSettings::_maxTemp = DEFAULT_MAX_TEMP;
int GlobalSettings::_maxLen = DEFAULT_MAX_PRIMER_LENGTH;
int GlobalSettings::_minLen = DEFAULT_MIN_PRIMER_LENGTH;
float GlobalSettings::_primerConcentration = DEFAULT_PRIMER_CONC;
float GlobalSettings::_monovalentIonConcentration = DEFAULT_SALT_CONC;
int GlobalSettings::_maxPrimers = DEFAULT_MAX_PRIMERS;
float GlobalSettings::_thermodynamicTemperature = DEFAULT_THERMODYNAMIC_TEMPERATURE;
bool GlobalSettings::_nonDegenerate = true;
bool GlobalSettings::_testRun = DEFAULT_RUN_TEST;
bool GlobalSettings::_invRevRun = false;
bool GlobalSettings::_SearchFwd = false;
bool GlobalSettings::_SearchRev = false;
bool GlobalSettings::_sortbytemp = true;
bool GlobalSettings::_userTemp = true;
string GlobalSettings::_testStr = "";
string GlobalSettings::_invRevValue = "";
string GlobalSettings::_searchFwdArg = "";
string GlobalSettings::_searchRevArg = "";
string GlobalSettings::_inputfile = "";
string GlobalSettings::_outputfile = "";


void ProcessTags(int argc, char *argv[]);
void PrintHelp();
string ConservedRegions(std::vector<Primer> primers);
string TestValue(DataSequence data, bool details);
string InvRev(DataSequence data);
string Banner(string message);

int main(int argc, char *argv[])
{
	
	// Check if user wants help or wants to test a k-mer
	if(argc == 1 || (argc == 2 && (strcmp("--h", argv[1]) == 0 || 
		strcmp("--help", argv[1]) == 0) ) )
	{
		PrintHelp();
	}
	/*
	else if(argc == 2 && (strstr(argv[1], "--test:") != NULL ||
		strstr(argv[1], "--invrev:") != NULL))
	{
		argc++;
	}*/

	// Process Tags
	ProcessTags(argc, argv);
	/*
	if(argc != 2)
	{
		ProcessTags(argc, argv);
	}*/

	// Create Filename/path/output strings
	//std::string filename = argv[argc - 1];
	std::string filename = GlobalSettings::GetInputFile();
	std::size_t found = filename.find_last_of(".");
	std::string primer_output = "";
	//std::string csv_output = "";
	std::string detail_output = "";
	std::string line_output = "";

	// Open Input File
	ifstream ifs;
	ifs.open(filename);
	if(ifs.fail())
	{
		cout << "Failure to open file.\n";
		exit(BAD_INPUT_FILE);
	}

	// Read Sequences
	SequenceReader read;
	SequenceList list = read.CreateList(ifs);
	// Test Sequences
	if(list.size() == 0)
	{
		cout << "There were no sequences in the input file.\n";
		exit(BAD_INPUT_FILE);
	}
	else if(GlobalSettings::GetProteinSequence())
	{
		ofstream os;
		//os.open(filename.substr(0, found) + "_protein.fasta");
		os.open(GlobalSettings::GetOutputFile());
		os << list.DecodeProteins() << endl;
		cout << "Decoded the proteins in the file.  Output saved to: ";
		//cout << filename.substr(0, found) + "_protein.faa" << endl;
		cout << GlobalSettings::GetOutputFile() << endl;
		exit(PROGRAM_SUCCESS);
	}
	else if(list.TestAlignment() == false)
	{
		cout << "Error. DeGenPrime cannot run on an unaligned list" << endl;
		cout << "Try: mafft --localpair --maxiterate 1000 --clustalout --quiet ";
		cout << GlobalSettings::GetInputFile() << " > ";
		string str = GlobalSettings::GetInputFile().substr(0, found);
		cout << str << ".clust" << endl;
		exit(PROGRAM_SUCCESS);
		/*
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
				if(strstr(argv[i], "--g") != NULL || strstr(argv[i], "--global") != NULL)
				{
					local = false;
				}
			}
		}
		string clustal = filename.substr(0,found);
		clustal += local ? "_l.clust" : "_g.clust";
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

		// Attempt to align the sequence
		try
		{
			system(c1);
		}
		catch(int exception)
		{
			cout << "There was an error trying to run MAFFT on this file." << endl;
			exit(FILE_MISALIGNED);
		}
		system(c2);
		exit(PROGRAM_SUCCESS);
		*/
	}
	else if(list.size() == 1)
	{
		GlobalSettings::SetNonDegenerate(false);
	}
	
	// Open Output File Stream
	ofstream ofs;
	//ofs.open(filename.substr(0, found) + ".dgp");
	if(GlobalSettings::GetOutputFile() == "")
	{
		string out = GlobalSettings::GetInputFile().substr(0, found) + ".csv";
		GlobalSettings::SetOutputFile(out);
	}
	ofs.open(GlobalSettings::GetOutputFile());
	if(ofs.fail())
	{
		cout << "Error opening output file.";
		exit(BAD_INPUT_FILE);
	}
	
	// Sequence Information Before Filtering
	detail_output += Banner(" Sequence Information ");
	detail_output += Format("Sequence Count: " + to_string(list.size()), STR_FORMAT, Alignment::Center);
	detail_output += "\n\n" + list.PrintSequenceNames() + "\n";

	// Make sure last Nucleotide isn't out-of-bounds
	int last = list.GetSequenceList()[0].size() < GlobalSettings::GetEndingNucleotide() ? 
		list.GetSequenceList()[0].size() : GlobalSettings::GetEndingNucleotide();
	GlobalSettings::SetEndingNucleotide(last);

	// Create Forward and Reverse DataSequences
	SequenceList reverse_list = list.InvRevList();
	DataSequence data = list.ProcessList();
	DataSequence rev = reverse_list.ProcessList();

	// Check DataSequence for conserved regions
	//detail_output += Banner(" Conserved Regions ");
	int conserved_count = 0;
	int begin_index = 0;
	bool conserved_start = false;
	bool conserved_region = false;
	std::vector<Primer> conserved_fwd_primers;
	for(int i = 0;i < data.size() - GlobalSettings::GetMinimumPrimerLength();i++)
	{
		char c = data.GetDataSequence()[i].GetCode();
		bool isConserved = (c == 'C' || c == 'G' || c == 'A' || c == 'T');
		if(isConserved)
		{
			conserved_count++;
			if(conserved_start == false)
			{
				conserved_start = true;
				begin_index = i;
			}
			if(conserved_count >= GlobalSettings::GetMinimumPrimerLength())
			{
				conserved_region = true;
			}
		}
		else 
		{
			if(conserved_start == true)
			{
				conserved_start = false;
				if(conserved_region)
				{
					Primer p(begin_index, conserved_count);
					conserved_fwd_primers.push_back(p);
				}
			}
			conserved_count = 0;
			conserved_region = false;
		}
		// Check the last region
		if(conserved_region && i == (data.size() - GlobalSettings::GetMinimumPrimerLength() - 1))
		{
			Primer p(begin_index, conserved_count);
			conserved_fwd_primers.push_back(p);
		}
	}
	/*
	cout << "Seg Fault after this point" << endl;
	cout << "Conserved_fwd_primers size = " << conserved_fwd_primers.size() << endl;
	line_output = "There ";
	line_output += (conserved_fwd_primers.size() == 1) ? "is " : "are ";
	line_output += to_string(conserved_fwd_primers.size());
	line_output += " conserved ";
	line_output += (conserved_fwd_primers.size() == 1) ? "region." : "regions.";
	detail_output += Format(line_output, STR_FORMAT, Alignment::Center) + "\n\n";
	*/
	if(conserved_fwd_primers.size() == 0 && GlobalSettings::GetNonDegenerate())
	{
		//line_output = "Insufficient conserved regions in sequences to find primers.";
		//detail_output += Format(line_output, STR_FORMAT, Alignment::Center) + "\n";
		std::vector<int> empty_int;
		//detail_output += data.Consensus(empty_int, empty_int, false);
		//cout << detail_output << endl;
		//ofs << detail_output << endl;
		ofs.close();
		cout << "Insufficient conserved regions for primers." << endl;
		exit(PROGRAM_SUCCESS);
	}
	else if(GlobalSettings::GetNonDegenerate())
	{
		//detail_output += ConservedRegions(conserved_fwd_primers) + "\n";
		int cand_pair_regions = 0;
		int amp;
		if(GlobalSettings::GetMeasureByAmpliconSize())
		{
			amp = GlobalSettings::GetMinimumAmplicon();
		}
		else
		{
			int first = (0 >= GlobalSettings::GetBeginningNucleotide()) ?
				0 : GlobalSettings::GetBeginningNucleotide();
			int last2 = (data.size() <= GlobalSettings::GetEndingNucleotide()) ?
				data.size() - 1 : GlobalSettings::GetEndingNucleotide();
			amp = last2 - first;
		}
		if(conserved_fwd_primers.size() == 1)
		{
			if(conserved_fwd_primers[0].Length() >= amp)
			{
				int region_len = ceil(((float)conserved_fwd_primers[0].Length())/((float)amp));
				cand_pair_regions += region_len;
			}
			else
			{
				cout << "The conserved region is smaller than minimum amplicon or size bounds." << endl;
				ofs << "The conserved region is smaller than minimum amplicon or size bounds." << endl;
				//std::vector<int> empty_int;
				//cout << data.Consensus(empty_int, empty_int, false) << endl;
				//ofs << data.Consensus(empty_int, empty_int, false) << endl;
				ofs.close();
				exit(PROGRAM_SUCCESS);
			}
		}
		else
		{
			/*
			// Is conserved region big enough to make multiple primer regions
			for(int i = 0;i < conserved_fwd_primers.size();i++)
			{
				if(conserved_fwd_primers[i].Length() >= amp)
				{
					int region_len = ceil(((float)conserved_fwd_primers[0].Length())/((float)amp));
					cand_pair_regions += region_len;
				}
			}*/

			// Is there enough space between the conserved regions to make primers.
			for(int i = 0;i < conserved_fwd_primers.size() - 1;i++)
			{
				for(int j = conserved_fwd_primers.size() - 1;j > i;j--)
				{
					int primer_amp = ((conserved_fwd_primers[j].Index() + 
						conserved_fwd_primers[j].Length()) - 
						conserved_fwd_primers[i].Index());
					if(primer_amp >= amp)
					{
						cand_pair_regions++;
					}
				}
			}
		}
		if(cand_pair_regions == 0)
		{
			//line_output = "Insufficient conservation to find primers within bounds.";
			//detail_output += Format(line_output, STR_FORMAT, Alignment::Left) + "\n";
			//std::vector<int> empty_int;
			//detail_output += data.Consensus(empty_int, empty_int, false);
			//cout << detail_output;
			//ofs << detail_output;
			ofs.close();
			exit(PROGRAM_SUCCESS);
		}
		/*
		else
		{
			line_output = "There ";
			line_output += (cand_pair_regions == 1) ? "is " : "are ";
			line_output += to_string(cand_pair_regions);
			line_output += " candidate pair ";
			line_output += (cand_pair_regions == 1) ? "mapping." : "mappings.";
			detail_output += Format(line_output, STR_FORMAT, Alignment::Center) + "\n\n";
		}*/
	}

	// Consensus Sequence Information
	std::vector<int> Indeces, Ranges;
	std::vector<Primer> conserved_rev_primers;
	for(int i = 0;i < conserved_fwd_primers.size();i++)
	{
		Indeces.push_back(conserved_fwd_primers[i].Index());
		Ranges.push_back(conserved_fwd_primers[i].Length());
		Primer r(data.RevIndex(conserved_fwd_primers[i].Index() + 
			conserved_fwd_primers[i].Length() - 1), conserved_fwd_primers[i].Length());
		conserved_rev_primers.push_back(r);
	}
	//detail_output += Banner(" Consensus Sequence ");
	//detail_output += data.Consensus(Indeces, Ranges, true) + "\n";
	
	// Create Primer Calculators
	PrimerCalculator calc, rev_calc;
	if(GlobalSettings::GetNonDegenerate()) // User wants conserved regions
	{
		calc.InitializeFromRegion(conserved_fwd_primers, data);
		rev_calc.InitializeFromRegion(conserved_rev_primers, rev);
	}
	else if(GlobalSettings::GetMeasureByAmpliconSize()) // User wants minimum amplicon length
	{
		calc.InitializePrimers(data);
		rev_calc.InitializePrimers(rev);
	}
	else // User specified a range
	{
		calc.InitializeBoundedPrimers(data, GlobalSettings::GetBeginningNucleotide());
		int rev_lowerBound = data.RevIndex(GlobalSettings::GetEndingNucleotide());
		rev_calc.InitializeBoundedPrimers(rev, rev_lowerBound);

		int range = GlobalSettings::GetEndingNucleotide() - GlobalSettings::GetBeginningNucleotide();
		int amp = GlobalSettings::GetMinimumAmplicon();
		int min = (range >= amp) ? range : amp;
		if(data.size() < min)
		{
			GlobalSettings::SetMinimumAmplicon(data.size());
		}
		else
		{
			GlobalSettings::SetMinimumAmplicon(min);
		}
	}

	// Abort program if either calculator has no primers
	if(calc.size() == 0 || rev_calc.size() == 0)
	{
		//line_output = "There were insufficient primers found for this data.";
		//detail_output += Format(line_output, STR_FORMAT, Alignment::Left) + "\n";
		//cout << detail_output;
		//ofs << detail_output;
		//ofs << line_output;
		ofs.close();
		cout << "Insufficient primers found for this data." << endl;
		exit(NO_PRIMERS_FOUND);
	}

	// Display number of possible primers, run filters and output filter percentages.
	/*
	if(GlobalSettings::GetNonDegenerate() == false)
	{
		detail_output += Banner(" Forward Primers ");
		detail_output += calc.FilterAll(data) + "\n";
		detail_output += Banner(" Reverse Primers ");
		detail_output += rev_calc.FilterAll(rev) + "\n";
	}*/

	if(GlobalSettings::GetSearchFwd() || GlobalSettings::GetSearchRev())
	{
		int index;
		string d_output = "";
		string l_output = "";
		d_output += Banner(" Search Mode ");
		if(GlobalSettings::GetSearchFwd())
		{
			if(calc.size() == 0)
			{
				index = -1;
				l_output = "No primers found in the calculator.";
			}
			else
			{
				index = calc.IndexOf(data, GlobalSettings::GetSearchFwdArg());
				l_output = "Forward Primer: \'" + GlobalSettings::GetSearchFwdArg();
			}
			if(index != -1)
			{
				l_output += "\' found at index: " + to_string(index);
				d_output += Format(l_output, STR_FORMAT, Alignment::Center) + "\n";
			}
			else
			{
				l_output += "\' not found.";
				d_output += Format(l_output, STR_FORMAT, Alignment::Center) + "\n";
				l_output = "\'" + GlobalSettings::GetSearchFwdArg() + "\' details:";
				d_output += Format(l_output, STR_FORMAT, Alignment::Left) + "\n";
				DataSequence d(GlobalSettings::GetSearchFwdArg());
				d_output += TestValue(d, true) + "\n";
			}
		}
		if(GlobalSettings::GetSearchRev())
		{
			if(rev_calc.size() == 0)
			{
				index = -1;
			}
			else
			{
				index = rev_calc.IndexOf(rev, GlobalSettings::GetSearchRevArg());
				l_output = "Reverse Primer: \'" + GlobalSettings::GetSearchRevArg();
			}
			if(index != -1)
			{
				l_output += "\' found at index: " + to_string(index);
				d_output += Format(l_output, STR_FORMAT, Alignment::Left) + "\n";
			}
			else
			{
				l_output += "\' not found.";
				d_output += Format(l_output, STR_FORMAT, Alignment::Left) + "\n";
				l_output = "\'" + GlobalSettings::GetSearchRevArg() + "\' details:";
				d_output += Format(l_output, STR_FORMAT, Alignment::Left) + "\n";
				DataSequence d(GlobalSettings::GetSearchRevArg());
				d_output += TestValue(d, true) + "\n";
			}
		}
		//ofs << detail_output;
		ifs.close();
		ofs.close();
		cout << d_output << endl;
		exit(TEST_MODE);
	}

	// Get Partitions of PrimerPairList
	PrimerPairList maxlist, pairlist, top;
	if(calc.size() != 0 && rev_calc.size() != 0)
	{
		//detail_output += Banner(" Primer Pair Filtering ");
		//detail_output += maxlist.CreateList(data, rev, calc.GetPrimers(), rev_calc.GetPrimers());
		const int part = pairlist.PartitionCount(calc.size(), rev_calc.size());
		//const int part = maxlist.PartitionCount();
		const int len_part = sqrt(part);
		bool done = false;
		bool move_horiz = false;
		int count = 1;
		int fwd_len = calc.size() / len_part;
		int rev_len = rev_calc.size() / len_part;
		int x_mult = 0;
		int y_mult = 0;
		int x_start = 0;
		int y_start = 0;
		int x_end = x_start + fwd_len;
		int y_end = y_start + rev_len;

		// Declare Filtering variables
		const int desiredpairs = GlobalSettings::GetMaximumReturnPrimers();
		int limit = pow(desiredpairs, 3);
		if(limit > part)limit = part;
		/*if(part != 1)
		{
			line_output = "Dividing pair lists into " + to_string(part) + " partitions.";
			detail_output += Format(line_output, STR_FORMAT, Alignment::Center) + "\n\n";
			if(limit < part)
			{
				line_output = "Checking top " + to_string(limit) + " partitions.";
				detail_output += Format(line_output, STR_FORMAT, Alignment::Center) + "\n\n";
			}
		}*/
		int nextIndex, filtercount, remaining, nextlength;
		int goodprimers = 0;

		// Run data partition in loop
		while(count <= part && count <= limit)
		{
			if(done)
			{
				break;
			}
			bool isSquare = ceil((double)sqrt(count)) == floor((double)sqrt(count));
			bool isOneLessThanSquare = ceil((double)sqrt(count + 1))
				== floor((double)sqrt(count + 1));

			// Create Primer Pair List
			pairlist.CreateFromRange(data, rev, calc.GetPrimers(),
				rev_calc.GetPrimers(), x_start, x_end, y_start, y_end);
			//line_output = "Partition #" + to_string(count);
			//line_output += " Pairs: " + to_string(pairlist.size());
			//detail_output += Format(line_output, STR_FORMAT, Alignment::Left) + "\n";

			// Set next data block partition
			if(count == 1)
			{
				move_horiz = true;
				x_mult = 1;
			}
			else if(isSquare)
			{
				y_mult = 0;
				x_mult = sqrt((double)count);
				move_horiz = true;
			}
			else if(isOneLessThanSquare)
			{
				x_mult++;
				move_horiz = true;
			}
			else if(move_horiz)
			{
				int place_holder = y_mult;
				y_mult = x_mult;
				x_mult = place_holder;
				move_horiz = false;
			}
			else
			{
				int place_holder = y_mult;
				y_mult = x_mult + 1;
				x_mult = y_mult;
				move_horiz = true;
			}

			// Adjust parameters
			x_start = x_mult * fwd_len;
			y_start = y_mult * rev_len;
			x_end = x_start + fwd_len > calc.size() ? calc.size() : x_start + fwd_len;
			y_end = y_start + rev_len > rev_calc.size() ? rev_calc.size() : y_start + rev_len;
			count++;
			
			
			// Filter PrimerPairList
			/*
			detail_output += pairlist.FilterMessage("start", 0) + "\n";
			detail_output += pairlist.FilterAmpliconLength() + "\n";
			detail_output += pairlist.FilterTemperatureDifference() + "\n";
			detail_output += pairlist.FilterMessage("final", 0) + "\n";
			*/

			// Make sure we still have primers to work with
			if(pairlist.size() == 0)
			{
				continue;
			}

			// Sort PrimerPairList by least temperature difference.
			pairlist.Sort();

			// Loop through top desired primer pairs to filter them for annealing temperature
			// and output the final list of primer pairs.
			remaining = pairlist.size();
			nextIndex = 0;
			filtercount = 0;
			nextlength = desiredpairs - goodprimers;
			if(nextlength > remaining)
			{
				nextlength = remaining;
			}
			do
			{
				PrimerPairList subPairList = pairlist.SubList(nextIndex, nextlength);
				remaining -= nextlength;
				nextIndex += nextlength;
				top.Append(subPairList);
				filtercount = top.FilterAnnealingTemp(data, rev, goodprimers);
				filtercount += top.FilterUnique();
				goodprimers = MAX_PRIMER_RETURNS - filtercount;
				nextlength = (filtercount < remaining) ? filtercount : remaining;
			} while(filtercount != 0 && remaining > 0);

			
			if(top.size() < desiredpairs)
			{
				nextIndex = 0;
				nextlength = desiredpairs - top.size();
			}
			else
			{
				done = true;
			}
		}
	}
	else // At least one of the fwd or rev lists was empty
	{
		//line_output = "At least one of the forward or reverse primer lists was empty.";
		//detail_output += Format(line_output, STR_FORMAT, Alignment::Left) + "\n";
		//ofs << detail_output;
		ofs.close();
		cout << "At least one of the forward or reverse primer lists was empty." << endl;
		exit(PROGRAM_SUCCESS);
	}

	// Print output
	//primer_output += Banner(" Results ");
	if(top.size() == 0)
	{
		//line_output = "No primer pairs were found for these specifications.";
		//primer_output += Format(line_output, STR_FORMAT, Alignment::Center) + "\n";
		cout << "No primer pairs were found for these specifications." << endl;
	}
	else
	{
		//line_output = top.PrintAll(data, rev);
		//primer_output += line_output;
		//ofs_csv.open(filename.substr(0, found) + ".csv");
		//ofs_csv << top.CreateCSV(data, rev);
		//ofs_csv.close();
		ofs << top.CreateCSV(data, rev);
	}
	//ofs << primer_output << endl;
	//ofs << detail_output << endl;

	// Close input/output file streams.
	ifs.close();
	ofs.close();
	
	// Show closing messages then close the program.
	//cout << "Output details saved to " << filename.substr(0, found) << ".dgp" << endl;
	cout << "Output details saved to " << GlobalSettings::GetOutputFile() << endl;
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
	//for(int i = 1;i < argc - 1;i++)
	for(int i = 1;i < argc;i++)
	{
		if(strcmp("--h", argv[i]) == 0 || strcmp("--help", argv[i]) == 0)
		{
			PrintHelp();
		}
		else if(strstr(argv[i], "--protein") != NULL)
		{
			GlobalSettings::SetProteinSequence(true);
			continue;
		}
		else if(strstr(argv[i], "--degenerate") != NULL)
		{
			GlobalSettings::SetNonDegenerate(false);
			continue;
		}
		else if(strstr(argv[i], "--global") != NULL ||
			strstr(argv[i], "--local") != NULL ||
			strstr(argv[i], "--g") != NULL ||
			strstr(argv[i], "--l") != NULL)
		{
			continue;
		}

		// For most tags, they follow the format: '--tag:value'
		// We want to check what the tag is with strstr()
		// and store the given value to global settings with strchr()
		ptr = strchr(argv[i], ':') + 1;
		if(ptr != NULL && (strstr(argv[i], "--test:") == NULL || 
			(strstr(argv[i], "--search_fwd:") == NULL) ||
			(strstr(argv[i], "--search_rev:") == NULL)))
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
		else if(strstr(argv[i], "--delta_g:") != NULL)
		{
			GlobalSettings::SetDeltaG(value);
		}
		else if(strstr(argv[i], "--end:") != NULL)
		{
			GlobalSettings::SetEndingNucleotide(value);
			containsEnd = true;
		}
		else if(strstr(argv[i], "--max_primer_len:") != NULL ||
			strstr(argv[i], "--max_primer_length:") != NULL)
		{
			GlobalSettings::SetMaximumPrimerLength(value);
		}
		else if(strstr(argv[i], "--min_primer_len:") != NULL ||
			strstr(argv[i], "--min_primer_length:") != NULL)
		{
			GlobalSettings::SetMinimumPrimerLength(value);
		}
		else if(strstr(argv[i], "--min_temp:") != NULL)
		{
			GlobalSettings::SetMinimumTemperature(value);
			GlobalSettings::SetUserTemp(true);
		}
		else if(strstr(argv[i], "--max_temp:") != NULL)
		{
			GlobalSettings::SetMaximumTemperature(value);
			GlobalSettings::SetUserTemp(true);
		}
		else if(strstr(argv[i], "--primer_conc:") != NULL)
		{
			GlobalSettings::SetPrimerConcentration(value);
		}
		else if(strstr(argv[i], "--salt_conc:") != NULL)
		{
			GlobalSettings::SetMonoIonConcentration(value);
		}
		else if(strstr(argv[i], "--search_fwd:") != NULL)
		{
			string str = ptr;
			GlobalSettings::SetSearchFwd(true);
			GlobalSettings::SetSearchFwdArg(str);
		}
		else if(strstr(argv[i], "--search_rev:") != NULL)
		{
			string str = ptr;
			GlobalSettings::SetSearchRev(true);
			GlobalSettings::SetSearchRevArg(str);
		}
		else if(strstr(argv[i], "--max_primers:") != NULL)
		{
			GlobalSettings::SetMaximumReturnPrimers(value);
		}
		else if(strstr(argv[i], "--test:") != NULL)
		{
			string str = ptr;
			GlobalSettings::SetTestValue(str);
			GlobalSettings::SetRunTest(true);
		}
		else if(strstr(argv[i], "--invrev:") != NULL)
		{
			string str = ptr;
			GlobalSettings::SetInvRevValue(str);
			GlobalSettings::SetRunInvRev(true);
		}
		else if(strstr(argv[i], "--input_file:") != NULL)
		{
			string str = ptr;
			GlobalSettings::SetInputFile(str);
		}
		else if(strstr(argv[i], "--output_file:") != NULL)
		{
			string str = ptr;
			GlobalSettings::SetOutputFile(str);
		}
		else
		{
			cout << "Warning: Unrecognized tag \'" << argv[i] << "\'" << endl;
			continue;
		}
	}

	// We need to check here if user has entered improper values for their tags.
	// close the program if any improper tags were entered.

	// The tags '--test' and '--search' are incompatible.
	if(GlobalSettings::GetRunTest() && (GlobalSettings::GetSearchFwd() ||
		GlobalSettings::GetSearchRev()))
	{
		cout << "ERROR: '--search' tags are incompatiable with ";
		cout << "the '--test' tag." << endl;
		exit(SETTINGS_FILE_NOT_FOUND);
	}

	// Make sure the user hasn't specified a min primer size > max primer size
	if(GlobalSettings::GetMaximumPrimerLength() < GlobalSettings::GetMinimumPrimerLength())
	{
		cout << "ERROR: cannot set maximum primer length < ";
		cout << "minimum primer length." << endl;
		exit(SETTINGS_FILE_NOT_FOUND);
	}

	if(GlobalSettings::GetMaximumPrimerLength() < MIN_PRIMER_LENGTH)
	{
		cout << "ERROR: cannot set a maximum primer length ";
		cout << "less than " << MIN_PRIMER_LENGTH << "." << endl;
		exit(SETTINGS_FILE_NOT_FOUND);
	}

	if(GlobalSettings::GetMinimumPrimerLength() > MAX_PRIMER_LENGTH)
	{
		cout << "ERROR: cannot set a minimum primer length ";
		cout << "greater than " << MAX_PRIMER_LENGTH << "." << endl;
		exit(SETTINGS_FILE_NOT_FOUND);
	}

	if(GlobalSettings::GetRunTest())
	{
		DataSequence data;
		for(char c : GlobalSettings::GetTestValue())
		{
			char c2;
			if(c == 'I' || c == 'i')c2 = 'G';
			else c2 = c;
			DataNode node(c2,c2,1.0);
			data.PushBack(node);
		}
		string message = TestValue(data, true);
		cout << message;
		exit(TEST_MODE);
	}
	else if(GlobalSettings::GetRunInvRev())
	{
		DataSequence data;
		for(char c : GlobalSettings::GetInvRevValue())
		{
			char c2;
			if(c == 'I' || c == 'i')c2 = 'G';
			else c2 = c;
			DataNode node(c2, c2, 1.0);
			data.PushBack(node);
		}
		string message = InvRev(data);
		cout << message << endl;
		exit(PROGRAM_SUCCESS);
	}

	// User cannot enter a starting nucleotide without also entering an ending nucleotide
	if(containsBegin != containsEnd)
	{
		if(containsBegin)
		{
			GlobalSettings::SetEndFlag(true);
		}
		else
		{
			GlobalSettings::SetBeginFlag(true);
		}
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
	if(containsBegin || containsEnd)
	{
		GlobalSettings::SetMeasureByAmpliconSize(false);
	}
}

void PrintHelp()
{
	cout << "degenprime version " << DeGenPrime_VERSION << "\n";
	cout << "Syntax: ./degenprime [--tags] <filename>\n";
	cout << "Valid tags include:\n";
	cout << "\t--amplicon:int, Set the minimum amplicon length.  This will not work with";
	cout << "--begin or --end tags.\n";
	cout << "\t--begin:int, Set the beginning nucleotide.  This must be used with --end ";
	cout << "and cannot be used with --amplicon.\n";
	cout << "\t--end:int, Set the ending nucleotide.  This must be used with --begin ";
	cout << "and cannot be used with --amplicon.\n";
	cout << "\t--global or --g, for lists of sequences that are misaligned, this tag specifies ";
	cout << "that the file should run MAFFT for global alignment.\n";
	cout << "\t--help or --h, prints this help menu.\n";
	cout << "\t--local or --l, for lists of sequences that are misaligned, this tag specifies ";
	cout << "that the file should run MAFFT for local alignment.\n";
	cout << "\t--min_temp:int, Sets the minimum primer melting temperature.  This has";
	cout << " a minimum value of " << MIN_PRIMER_TEMP << " (degrees Celsius) and must be ";
	cout << "smaller than --max_temp.\n";
	cout << "\t--max_primers:int, Sets the maximum number of output primers.  This has ";
	cout << "a maximum value of " << MAX_PRIMER_RETURNS << " and this program will reduce ";
	cout << "any value larger to this value.\n";
	cout << "\t--max_temp:int, Sets the maximum primer melting temperature.  This has";
	cout << " a maximum value of " << MAX_PRIMER_TEMP << " (degrees Celsius) and must be ";
	cout << "larger than --min_temp.\n";
	cout << "\t--primer_conc:int, Sets the concentration of the PCR primer in nM.  This has ";
	cout << "a minimum value of " << MIN_PRIMER_CONC << " and this program will raise ";
	cout << "any value smaller to this value.\n";
	cout << "\t--protein, Tells the program that the input sequence is a protein sequence and ";
	cout << "the program should unwrap the protein sequence into its base nucleotides.  This ";
	cout << "will produce degenerate nucleotide codes whenever there is any ambiguity.\n";
	cout << "\t--salt_conc:int, Sets the concentration of monovalent ions in mM.  This has ";
	cout << "a minimum value of " << MIN_SALT_CONC << " and this program will raise ";
	cout << "any value smaller to this value.\n";
	cout << "\t--test:string, The string represents a single primer.  Runs the primer through all filters.";
	cout <<  "Returns the thermodynamic values of this primer as well as any filters this primer would not pass and";
	cout << "its calculated penalty.  This tag is incompatible with --search tags.";
	cout << "Any primer smaller or larger than the size limits will show primer outside size range.";


	exit(PROGRAM_SUCCESS);
}

string TestValue(DataSequence data, bool details)
{
	string message = "";
	string line = "";
	if(details)message += data.Print() + "\n";
	PrimerCalculator prime, clone;
	prime.InitializeTestPrimer(data);
	clone.InitializeTestPrimer(data);
	bool flag = false;
	if(data.isEmpty())
	{
		message += "Primer is located in a non conserved region.\n";
	}
	string valid_chars = "AaCcGgTtWwSsRrYyKkMmBbDdHhVvNn-";
	for(int i = 0;i < data.size();i++)
	{
		if(flag)
		{
			break;
		}
		char c = data.GetDataSequence()[i].GetCode();
		if(valid_chars.find(c) == string::npos)
		{
			flag = true;
		}
	}
	if(flag)
	{
		message = "Error. Primer contained one or more invalid characters.\n";
		return message;
	}
	else if(data.size() < MIN_PRIMER_LENGTH || data.size() > MAX_PRIMER_LENGTH)
	{
		message += "Primer not within allowed size for filtering.\n";
	}
	else
	{
		line = "Penalty: ";
		line += Format(data.Penalty(), 2);
		message += line + "\n";
		message += "Filtered by:\n";
		string trash = prime.FilterDeletions(data);
		if(prime.size() < clone.size())
		{
			prime = clone;
			message += "\tFilterDeletions\n";
			flag = true;
		}
		trash = prime.FilterDegeneracy(data);
		if(prime.size() < clone.size())
		{
			prime = clone;
			message += "\tFilterDegeneracy (";
			int count = 0;
			for(int i = 0;i < data.size();i++)
			{
				char c = data.GetDataSequence()[i].GetCode();
				if(c == 'A' || c == 'C' || c == 'G'
					|| c == 'T' || c == '-')
				{
					continue;
				}
				else
				{
					count++;
				}
			}
			message += to_string(count) + " out of " + to_string(data.size());
			message += ")\n";
			flag = true;
		}
		trash = prime.FilterGCContent(data);
		if(prime.size() < clone.size())
		{
			prime = clone;
			message += "\tFilterGCContent (";
			float ratio = data.GCRatio();
			bool space = false;
			if(ratio < MIN_GC_TOTAL_RATIO || ratio > MAX_GC_TOTAL_RATIO)
			{
				line = "GC: ";
				line += Format((float)100.0 * data.GCRatio(), 2) + "%";
				message += line;
				space = true;
			}
			DataSequence ending = data.SubSeq(data.size() - 5, 5);
			ratio = ending.GCRatio();
			if(ratio > MAX_GC_EXTREMA_RATIO)
			{
				line = "";
				if(space)line += " ";
				line += "End GC: ";
				line += Format((float)100.0 * ratio, 2) + "%";
				message += line;
			}
			message += ")\n";
			flag = true;
		}
		trash = prime.FilterRepeats(data);
		if(prime.size() < clone.size())
		{
			prime = clone;
			message += "\tFilterRepeats\n";
			flag = true;
		}
		trash = prime.FilterComplementaryEnds(data);
		if(prime.size() < clone.size())
		{
			prime = clone;
			message += "\tFilterComplementaryEnds\n";
			flag = true;
		}
		trash = prime.FilterTemperature(data, 0.0);
		if(prime.size() < clone.size())
		{
			prime = clone;
			message += "\tFilterTemperature\n";
			flag = true;
		}
		if(flag == false)
		{
			message += "\t(None)\n";
		}
	}
	return message;
}

string InvRev(DataSequence data)
{
	DataSequence rev = data.RevSeq();
	DataSequence revinv = rev.InvSeq();
	string ret = revinv.Codes();
	return ret;
}

string Banner(string message)
{
	string ret = "";
	string line = "";
	int len = message.length() + 2;
	for(int i = 0;i < len;i++)
	{
		line += "%";
	}
	string line2 = Format(line, STR_FORMAT, Alignment::Center);
	string line3 = Format("%" + message + "%", STR_FORMAT, Alignment::Center);
	ret += line2 + "\n" + line3 + "\n" + line2 + "\n\n";
	return ret;
}

string ConservedRegions(std::vector<Primer> primers)
{
	string ret = "";
	int largest = 0;
	for(Primer p : primers)
	{
		int final_ind = p.Index() + p.Length() - 1;
		if(final_ind > largest)
		{
			largest = final_ind;
		}
	}
	int digs = digits(largest);
	// ' (####-####) '
	const int reg_size = (digs * 2) + 5;
	const int max_regs_per_line = STR_FORMAT / reg_size;
	const int square_threshold = (max_regs_per_line) * (max_regs_per_line - 1);
	int root = 1;
	while(pow(root, 2) < primers.size())
	{
		root++;
	}
	int regs_per_line = (pow(root, 2) > square_threshold) ? max_regs_per_line : root;
	int form_size = STR_FORMAT / regs_per_line;
	//int reg_count = 0;
	string region = "";
	for(int i = 0;i < primers.size();i++)
	{
		if(i != 0 && i % regs_per_line == 0)
		{
			ret += "\n\n";
		}
		int start = primers[i].Index();
		int end = primers[i].Index() + primers[i].Length() - 1;
		region = "(" + Format(start, digs) + "-";
		region += Format(end, digs) + ")";
		ret += Format(region, form_size, Alignment::Center);
	}
	ret += "\n";
	return ret;
}
