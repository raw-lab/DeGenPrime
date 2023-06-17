#include <cstdlib>
#include <cstring>
#include <cmath>
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
bool GlobalSettings::_proteinSequence = DEFAULT_PROTEIN_SEQUENCE;
bool GlobalSettings::_beginflag = DEFAULT_BEGIN_FLAG;
bool GlobalSettings::_endflag = DEFAULT_END_FLAG;
float GlobalSettings::_minTemp = DEFAULT_MIN_TEMP;
float GlobalSettings::_maxTemp = DEFAULT_MAX_TEMP;
float GlobalSettings::_primerConcentration = DEFAULT_PRIMER_CONC;
float GlobalSettings::_monovalentIonConcentration = DEFAULT_SALT_CONC;
int GlobalSettings::_maxPrimers = DEFAULT_MAX_PRIMERS;
float GlobalSettings::_thermodynamicTemperature = DEFAULT_THERMODYNAMIC_TEMPERATURE;
bool GlobalSettings::_testRun = DEFAULT_RUN_TEST;
bool GlobalSettings::_SearchFwd = false;
bool GlobalSettings::_SearchRev = false;
bool GlobalSettings::_DoSearchFile = false;
bool GlobalSettings::_sortbytemp = true;
string GlobalSettings::_testStr = "";
string GlobalSettings::_searchFwdArg = "";
string GlobalSettings::_searchRevArg = "";
string GlobalSettings::_searchFile = "";


void ProcessTags(int argc, char *argv[]);
void PrintHelp();
string TestValue(DataSequence data, bool details);
string Analysis(ifstream& ifs);

int main(int argc, char *argv[])
{
	// Set up clock
	struct timeval begin, end;
	gettimeofday(&begin, 0);
	
	// Check if user wants help or wants to test a k-mer
	if(argc == 1 || (argc == 2 && (strcmp("--h", argv[1]) == 0 || 
		strcmp("--help", argv[1]) == 0) ) )
	{
		PrintHelp();
	}
	else if(argc == 2 && strstr(argv[1], "--test:") != NULL)
	{
		argc++;
	}

	// Process Tags
	if(argc != 2)
	{
		ProcessTags(argc, argv);
	}

	// Create Filename/path
	string filename = argv[argc - 1];
	string filepath = std::filesystem::current_path();
	std::size_t found = filename.find_last_of(".");

	// Open Input File
	ifstream ifs;
	ifs.open(filepath + "/" + filename);
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
		std::size_t found = filename.find_last_of(".");
		os.open(filepath + "/" + filename.substr(0, found) + "_protein.fasta");
		cout << "Decoded proteins: " << list.DecodeProteins() << endl;
		os << "" << list.DecodeProteins() << endl;
		cout << "Decoded the proteins in the file.  Output saved to: ";
		cout << filepath + "/" + filename.substr(0, found) + "_protein.fasta" << endl;
		exit(PROGRAM_SUCCESS);
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
				if(strstr(argv[i], "--g") != NULL || strstr(argv[i], "--global") != NULL)
				{
					local = false;
				}
			}
		}
		std::size_t found = filename.find_last_of(".");
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
		exit(FILE_MISALIGNED);
	}

	// Open Output File Stream
	ofstream ofs;
	ofs.open(filepath + "/" + filename.substr(0, found) + ".dgp");

	// File Information
	ofs << "Full path of argument: " << filepath << endl;
	ofs << "Full path of file: " << filepath << "/" << filename << endl;
	ofs << "Input File opened." << endl;
	
	// Sequence Information Before Filtering
	ofs << "Before Sequence Filters, the number of Sequences are: " << list.size() << endl;
	ofs << list.PrintSequenceNames();
	int last = list.GetSequenceList()[0].size() < GlobalSettings::GetEndingNucleotide() ? 
		list.GetSequenceList()[0].size() : GlobalSettings::GetEndingNucleotide();
	GlobalSettings::SetEndingNucleotide(last);
	SequenceList filter_list = list.FilterDashes();
	
	// Sequence Information After Filtering
	if(filter_list.size() != 0)
	{	
		ofs << "\nAfter Sequence Filters, these " << list.size() " sequences were filtered."  << endl;
		ofs << filter_list.PrintSequenceNames() << endl;
	}

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
		calc.InitializeBoundedPrimers(data, GlobalSettings::GetBeginningNucleotide());
		int rev_lowerBound = data.RevIndex(GlobalSettings::GetEndingNucleotide());
		cout << "value of rev_Lowerbound: " << rev_lowerBound << endl;
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

	// Display number of possible primers, run filters and output filter percentages.
	cout << "Number of possible forward primers in specified size range, before Filters: [" << calc.size() << "]\n" << endl;
	ofs << "Number of possible forward primers in specified size range, before Filters: [" << calc.size() << "]\n" << endl;
	ofs << calc.FilterAll(data) << endl;
	cout << "After running Filters: [" << calc.size() << "]\n" << endl;

	cout << "Number of possible reverse primers in specified size range, before Filters: [" << rev_calc.size() << "]\n" << endl;
	ofs << "Number of possible reverse primers in specified size range, before Filters: [" << rev_calc.size() << "]\n" << endl;
	ofs << rev_calc.FilterAll(rev) << endl;
	cout << "After running Filters: [" << rev_calc.size() << "]\n" << endl;
	
	// Check PrimerPairs Can Exist and Not Searching
	if(calc.size() == 0 || rev_calc.size() == 0)
	{
		cout << "No primer pairs were found that meet your specifications." << endl;
		ofs << "No primer pairs were found that meet your specifications." << endl;
		if(GlobalSettings::GetDoSearchFile() == false &&
			GlobalSettings::GetSearchFwd() == false &&
			GlobalSettings::GetSearchRev() == false)
		{
			exit(PROGRAM_SUCCESS);
		}
	}
	cout << "Total combinations of primer pairs:" << (calc.size() * rev_calc.size()) << endl;

	// Assign Quality to Primers and Sort
	for(Primer prime : calc.GetPrimers())
	{
		DataSequence sub = data.SubSeq(prime.Index(), prime.Length());
		prime.SetQuality(sub.Quality());
	}
	for(Primer rev_prime : rev_calc.GetPrimers())
	{
		DataSequence rev_sub = rev.SubSeq(rev_prime.Index(), rev_prime.Length());
		rev_prime.SetQuality(rev_sub.Quality());
	}
	calc.Sort();
	rev_calc.Sort();

	if(GlobalSettings::GetSearchFwd() || GlobalSettings::GetSearchRev())
	{
		int index;
		if(GlobalSettings::GetSearchFwd())
		{
			index = calc.IndexOf(data, GlobalSettings::GetSearchFwdArg());
			cout << "Forward Primer: \'" << GlobalSettings::GetSearchFwdArg();
			if(index != -1)
			{
				cout << "\' found at index: " << index << endl;
			}
			else
			{
				cout << "\' not found." << endl;
				cout << "\'" << GlobalSettings::GetSearchFwdArg() << "\' details:" << endl;
				DataSequence d(GlobalSettings::GetSearchFwdArg());
				cout <<  TestValue(d, true) << endl;
			}
		}
		if(GlobalSettings::GetSearchRev())
		{
			index = rev_calc.IndexOf(rev, GlobalSettings::GetSearchRevArg());
			cout << "Reverse Primer: \'" << GlobalSettings::GetSearchRevArg();
			if(index != -1)
			{
				cout << "\' found at index: " << index << endl;
			}
			else
			{
				cout << "\' not found." << endl;
				cout << "\'" << GlobalSettings::GetSearchRevArg() << "\' details:" << endl;
				DataSequence d(GlobalSettings::GetSearchRevArg());
				cout <<  TestValue(d, true) << endl;
			}
		}

		ifs.close();
		ofs.close();
		exit(PROGRAM_SUCCESS);
	}

	// Get Partitions of PrimerPairList
	PrimerPairList pairlist, top;
	if(calc.size() != 0 && rev_calc.size() != 0)
	{
		const int part = pairlist.PartitionCount(calc.size(), rev_calc.size());
		if(part != 1)
		{
			cout << "Splitting into " << part;
			cout << " partitions with approximately 1600 primerpairs each.\n" << endl;
		}
		const int len_part = sqrt(part);
		bool start = true;
		bool done = false;
		bool move_horiz = false;
		int count = 1;
		int fwd_len = ceil((float)calc.size()/(float)len_part);
		int rev_len = ceil((float)rev_calc.size()/(float)len_part);
		int x_mult = 0;
		int y_mult = 0;
		int x_start = 0;
		int y_start = 0;
		int x_end = x_start + fwd_len;
		int y_end = y_start + rev_len;

		// Declare Filtering variables
		const int desiredpairs = GlobalSettings::GetMaximumReturnPrimers();
		int nextIndex, filtercount, remaining, goodprimers, nextlength;

		// Run data partition in loop
		while(count <= part)
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
			cout << "Partition #" << to_string(count);
			cout << " Pairlist size: " << to_string(pairlist.size()) << endl;
			ofs << "Partition #" << to_string(count);
			ofs << " Pairlist size: " << to_string(pairlist.size()) << endl;

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
			string filter_msg = pairlist.FilterAmpliconLength();
			ofs << filter_msg;
			cout << filter_msg;
			filter_msg = pairlist.FilterTemperatureDifference();
			ofs << filter_msg;
			cout << filter_msg;
			filter_msg = pairlist.FilterMessage("final", 0);
			ofs << filter_msg;
			cout << filter_msg;

			// Make sure we still have primers to work with
			if(pairlist.size() == 0)
			{
				ofs << "Moving to next partition." << endl;
				cout << "Moving to next partition." << endl;
				continue;
			}

			// Sort PrimerPairList by least temperature difference.
			ofs << "\nSorting remaining " << pairlist.size();
			ofs << " primer pairs in the list by temperature difference." << endl;
			cout << "\nSorting remaining " << pairlist.size();
			cout << " primer pairs in the list by temperature difference." << endl;
			pairlist.Sort();

			// Loop through top desired primer pairs to filter them for annealing temperature
			// and output the final list of primer pairs.

			ofs << "Running FilterAnnealingTemp on the top " << desiredpairs << " primer pairs." << endl;
			if(start)
			{
				top = pairlist.SubList(0,desiredpairs);
			}
			if(top.size() == 0)
			{
				continue;
			}
			else
			{
				remaining = pairlist.size() - desiredpairs;
				if(start)
				{
					nextIndex = desiredpairs;
					filtercount = top.FilterAnnealingTemp(data, rev, 0);
					goodprimers = MAX_PRIMER_RETURNS - filtercount;
					start = false;
				}
				while(filtercount != 0 && remaining > 0)
				{
					nextlength = (filtercount < remaining) ? filtercount : remaining;
					PrimerPairList next = pairlist.SubList(nextIndex,nextlength);
					filtercount = next.FilterAnnealingTemp(data, rev, goodprimers);
					goodprimers = MAX_PRIMER_RETURNS - filtercount;
					nextIndex += nextlength;
					remaining -= nextlength;
					top.Append(next);
				}
			}
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

	// Check if user wanted to run a search for primers then find
	// primers and exit the program.
	if(GlobalSettings::GetDoSearchFile())
	{
		string bould_file = filepath + "/" + GlobalSettings::GetSearchFile();
		ifstream b_file;
		b_file.open(bould_file);
		if(b_file.fail())
		{
			cout << "Error.  Could not open target search file \'" << bould_file << "\'" << endl;
			exit(BAD_INPUT_FILE);
		}
		cout << "Searching /'" << bould_file << "/' for sequences and primers.\n";
		ofs << "Searching /'" << bould_file << "/' for sequences and primers.\n";

		// The target input file is in boulder_io format.
		// The sequence names will begin with 'SEQUENCE_ID='
		// The forward primers will begin with 'PRIMER_<Direction>_<ID#>='
		//		followed by the primer nucleotide sequence
		
		string line = "";
		while(getline(b_file, line))
		{
			int p_l_count = 0;
			int p_r_count = 0;
			bool r_detail_flag = false;
			bool l_detail_flag = false;
			string message = "";
			string banner = "";
			if(line.find("SEQUENCE_ID=") == string::npos)
			{
				continue;
			}
			else
			{
				line.erase(0,12);
				message = "% Testing primers from sequence \'" + line + "\' %";
				for(int i = 0;i < message.length();i++)
				{
					banner += "%";
				}
				cout << banner << endl << message << endl << banner << endl << endl;
				ofs << banner << endl << message << endl << banner << endl << endl;

				string p_l = "PRIMER_LEFT_" + to_string(p_l_count) + "_SEQUENCE=";
				string p_r = "PRIMER_RIGHT_" + to_string(p_r_count) + "_SEQUENCE=";
				string p_l_data;
				string p_r_data;
				while(getline(b_file, line) && (((p_l_count < 5) && (p_r_count < 5)) ||
					r_detail_flag || l_detail_flag))
				{
					
					int index;
					if(line.find(p_l) == string::npos && line.find(p_r) == string::npos &&
						line.find(p_l_data) == string::npos && line.find(p_r_data) == string::npos)
					{
						continue;
					}
					if(line.find(p_l) != string::npos)
					{
						cout << line << endl;
						ofs << line << endl;
						line.erase(0,p_l.length());
						cout << "\'" << line << "\'";
						ofs << "\'" << line << "\'";
						index = calc.IndexOf(data, line);
						if(index != -1)
						{
							cout << " forward primer found at index: " << index << endl << endl;
							ofs << " forward primer found at index: " << index << endl << endl;
							p_l = "PRIMER_LEFT_" + to_string(++p_l_count) + "_SEQUENCE=";
						}
						else
						{
							cout << " forward primer not found.\n" << endl;
							ofs << " forward primer not found.\n" << endl;
							l_detail_flag = true;
							p_l_data = "PRIMER_LEFT_" + to_string(p_r_count) + "=";
							p_l = line;
						}
					}
					else if(line.find(p_r) != string::npos)
					{
						cout << line << endl;
						ofs << line << endl;
						line.erase(0,p_r.length());
						cout << "\'" << line << "\'";
						ofs << "\'" << line << "\'";
						index = rev_calc.IndexOf(rev, line);
						if(index != -1)
						{
							cout << " reverse primer found at index: " << index << endl << endl;
							ofs << " reverse primer found at index: " << index << endl << endl;
							p_r = "PRIMER_RIGHT_" + to_string(++p_r_count) + "_SEQUENCE=";
						}
						else
						{
							cout << " reverse primer not found.\n" << endl;
							ofs << " reverse primer not found.\n" << endl;
							r_detail_flag = true;
							p_r_data = "PRIMER_RIGHT_" + to_string(p_r_count) + "=";
							p_r = line;
						}
					}
					else if((line.find(p_l_data) != string::npos) && l_detail_flag)
					{
						cout << "Data of: " << line << endl;
						ofs << "Data of: " << line << endl;
						line.erase(0,p_l_data.length());
						std::size_t comma = line.find_last_of(",");
						index = stoi(line.substr(0,comma));
						line.erase(0,comma + 1);
						int len = stoi(line);
						DataSequence p_d = data.SubSeq(index, len);
						DataSequence l(p_l);
						cout << "Sequence Codes: " << p_d.Codes() << endl;
						ofs  << "Sequence Codes: " << p_d.Codes() << endl;
						cout << "   Sequence MC: " << p_d.MC() << endl;
						ofs  << "   Sequence MC: " << p_d.MC() << endl;
						cout << "     Primer MC: " << l.MC() << endl;
						ofs  << "     Primer MC: " << l.MC() << endl;
						cout << "                ";
						ofs  << "                ";
						int mismatch_count = 0;
						int deletions_count = 0;
						for(int i = 0;i < l.size();i++)
						{
							if(l.GetDataSequence()[i].GetMostCommon() == 
								p_d.GetDataSequence()[i].GetMostCommon())
							{
								cout << " ";
								ofs << " ";
							}
							else
							{
								cout << "*";
								ofs << "*";
								mismatch_count++;
							}
						}
						cout << endl;
						ofs << endl;
						if(mismatch_count != 0)
						{
							string match = "Number of mismatches: " + to_string(mismatch_count) + "\n";
							cout << match;
							ofs << match;
						}
						string test_val = TestValue(p_d, false);
						if(test_val.find("(None)") != string::npos)
						{
							cout << "The replacement primer \'" << p_d.MC() << "\' is good.\n";
							ofs << "The replacement primer \'" << p_d.MC() << "\' is good.\n";
						}
						cout << test_val << endl;
						ofs << test_val << endl;
						p_l = "PRIMER_LEFT_" + to_string(++p_l_count) + "_SEQUENCE=";
						l_detail_flag = false;
					}
					else if((line.find(p_r_data) != string::npos) && r_detail_flag)
					{
						cout << "Data of: " << line << endl;
						ofs << "Data of: " << line << endl;
						line.erase(0,p_r_data.length());
						std::size_t comma = line.find_last_of(",");
						index = stoi(line.substr(0,comma));
						index = rev.RevIndex(index);
						line.erase(0,comma + 1);
						int len = stoi(line);
						DataSequence p_d = rev.SubSeq(index, len);
						DataSequence r(p_r);
						cout << "Sequence Codes: " << p_d.Codes() << endl;
						ofs  << "Sequence Codes: " << p_d.Codes() << endl;
						cout << "   Sequence MC: " << p_d.MC() << endl;
						ofs  << "   Sequence MC: " << p_d.MC() << endl;
						cout << "     Primer MC: " << r.MC() << endl;
						ofs  << "     Primer MC: " << r.MC() << endl;
						cout << "                ";
						ofs  << "                ";
						int mismatch_count = 0;
						for(int i = 0;i < r.size();i++)
						{
							if(r.GetDataSequence()[i].GetMostCommon() == 
								p_d.GetDataSequence()[i].GetMostCommon())
							{
								cout << " ";
								ofs << " ";
							}
							else
							{
								cout << "*";
								ofs << "*";
								mismatch_count++;
							}
						}
						cout << endl;
						ofs << endl;
						if(mismatch_count != 0)
						{
							string match = "Number of mismatches: " + to_string(mismatch_count) + "\n";
							cout << match;
							ofs << match;
						}
						string test_val = TestValue(p_d, false);
						if(test_val.find("(None)") != string::npos)
						{
							cout << "The replacement primer \'" << p_d.MC() << "\' is good.\n";
							ofs << "The replacement primer \'" << p_d.MC() << "\' is good.\n";
						}
						cout << test_val << endl;
						ofs << test_val << endl;
						p_r = "PRIMER_RIGHT_" + to_string(++p_r_count) + "_SEQUENCE=";
						r_detail_flag = false;
					}
				} // end of Reading sequence primers
			} // End of reading sequence
		} // end of reading file
		b_file.close();

		ifstream input;
		input.open(filepath + "/" + filename.substr(0, found) + ".dgp");

		ofstream summary;
		summary.open(filepath + "/" + filename.substr(0, found) + "_summary.dgp");
		string analysis = Analysis(input);
		summary << analysis << endl;
		ofs << analysis << endl;
		
		input.close();
		summary.close();
	}

	// Close input/output file streams.
	ifs.close();
	ofs.close();
	
	// Show closing messages and clock to user, then close the program.
	cout << "Output details saved to primers_" << filename.substr(0, found) << ".txt" << endl;
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
		else if(strstr(argv[i], "--protein") != NULL)
		{
			GlobalSettings::SetProteinSequence(true);
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
		else if(strstr(argv[i], "--search_file:") != NULL)
		{
			string str = ptr;
			GlobalSettings::SetSearchFile(str);
			GlobalSettings::SetDoSearchFile(true);
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
		else
		{
			continue;
		}
	}

	// We need to check here if user has entered improper values for their tags.
	// close the program if any improper tags were entered.

	// The tags '--test' and '--search' are incompatible.
	if(GlobalSettings::GetRunTest() && (GlobalSettings::GetSearchFwd() ||
		GlobalSettings::GetRunTest() || GlobalSettings::GetDoSearchFile()))
	{
		cout << "ERROR: '--search' tags are incompatiable with ";
		cout << "the '--test' tag." << endl;
		exit(SETTINGS_FILE_NOT_FOUND);
	}

	// The tag '--search_file' is not usable with '--search_fwd' or '--search_rev'.
	if(GlobalSettings::GetDoSearchFile() && (GlobalSettings::GetSearchFwd() ||
		GlobalSettings::GetSearchRev()))
	{
		cout << "ERROR: '--search_file' tag is not compatible with ";
		cout << "'--search_fwd' or '--search_rev' tags." << endl;
		exit(SETTINGS_FILE_NOT_FOUND);
	}

	if(GlobalSettings::GetRunTest())
	{
		DataSequence data;
		for(char c : GlobalSettings::GetTestValue())
		{
			DataNode node(c,c,1.0);
			data.PushBack(node);
		}
		string message = TestValue(data, true);
		cout << message;
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
	cout << "Syntax: ./DeGenPrime [--tags] <filename>\n";
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
	cout << "\t--max_primers:int, Sets the maximum number of output primers.  This has ";
	cout << "a maximum value of " << MAX_PRIMER_RETURNS << " and this program will reduce ";
	cout << "any value larger to this value.\n";
	exit(PROGRAM_SUCCESS);
}

string TestValue(DataSequence data, bool details)
{
	string message = "";
	if(details)message += data.Print() + "\n";
	PrimerCalculator prime, clone;
	prime.InitializeTestPrimer(data);
	clone.InitializeTestPrimer(data);
	bool flag = false;
	if(data.isEmpty())
	{
		message += "Primer is empty.\n";
	}
	else if(data.size() < MIN_PRIMER_LENGTH || data.size() > MAX_PRIMER_LENGTH)
	{
		message += "Primer not within allowed size for filtering.\n";
	}
	else
	{
		message += "Filtered by:\n";
		string trash = prime.FilterDeletions(data);
		if(prime.size() < clone.size())
		{
			prime = clone;
			message += "\tFilterDeletions";
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
			if(ratio < MIN_GC_TOTAL_RATIO)
			{
				message += " primer GC ratio: " + to_string(data.GCRatio());
				message += " less than " + to_string(MIN_GC_TOTAL_RATIO);
			}
			else if(ratio > MAX_GC_TOTAL_RATIO)
			{
				message += "primer GC ratio: " + to_string(data.GCRatio());
				message += " more than " + to_string(MAX_GC_TOTAL_RATIO) + " ";
			}
			DataSequence ending = data.SubSeq(data.size() - 6, 5);
			ratio = ending.GCRatio();
			if(ratio > MAX_GC_EXTREMA_RATIO)
			{
				message += "primer end GC ratio: " + to_string(ratio);
				message += " more than " + to_string(MAX_GC_EXTREMA_RATIO) + " ";
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
		cout << "Past FilterRepeats" << endl;
		trash = prime.FilterComplementaryEnds(data);
		if(prime.size() < clone.size())
		{
			prime = clone;
			message += "\tFilterComplementaryEnds\n";
			flag = true;
		}
		trash = prime.FilterHairpins(data);
		if(prime.size() < clone.size())
		{
			prime = clone;
			message += "\tFilterHairpins\n";
			flag = true;
		}
		trash = prime.FilterDimers(data);
		if(prime.size() < clone.size())
		{
			prime = clone;
			message += "\tFilterDimers (";
			DataSequence ending = data.SubSeq(data.size() - 6, 5);
			message += "ending delta g: " + to_string(ending.Gibbs());
			message += ")\n";
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

string Analysis(ifstream& ifs)
{
	int seq_total = 0;
	int p_total = 0;
	int replace = 0;
	int good, degen, del, gc, rep, comp_end, hair, dimer, temp, mismatch;
	good = degen = del = gc = rep = comp_end = hair = dimer = temp = mismatch = -2;
	float ratio;

	string ret = "";
	string line = "";
	while(getline(ifs, line))
	{
		if(line.find("Testing primers from sequence") != string::npos)seq_total++;
		else if(line.find("_SEQUENCE=") != string::npos)p_total++;
		else if(line.find("primer found") != string::npos)good++;
		else if(line.find("replacement primer") != string::npos)replace++;
		else if(line.find("Number of mismatches:") != string::npos)mismatch++;
		else if(line.find("FilterDegeneracy") != string::npos)degen++;
		else if(line.find("FilterDeletions") != string::npos)del++;
		else if(line.find("FilterGCContent") != string::npos)gc++;
		else if(line.find("FilterRepeats") != string::npos)rep++;
		else if(line.find("FilterComplementaryEnds") != string::npos)comp_end++;
		else if(line.find("FilterHairpins") != string::npos)hair++;
		else if(line.find("FilterDimers") != string::npos)dimer++;
		else if(line.find("FilterTemperature") != string::npos)temp++;
	}
	good += 2;

	string banner = "";
	string message = "%    Analysis summary of suggested primers    %";
	for(int i = 0;i < message.length();i++)
	{
		banner += "%";
	}
	ret += banner + "\n" + message + "\n" + banner + "\n\n";
	ret += "\t\tGeneral\n";
	ret += "                Total Sequences: " + to_string(seq_total) + "\n";
	ret += "                  Total Primers: " + to_string(p_total) + "\n\n";
	ret += "                Passing Primers: " + to_string(good) + "\n";
	ret += " Acceptable Replacement Primers: " + to_string(replace) + "\n";
	ret += "                Failing Primers: " + to_string(p_total - (good + replace)) + "\n";
	ratio = ((float)(good + replace))/((float)p_total);
	ret += "                  Passing Ratio: " + to_string(ratio * 100.0) + "%\n\n";
	ret += "\t\tFailure Breakdown\n";
	ratio = ((float)mismatch)/((float)p_total) * 100.0;
	ret += "\t  Mismatches: " + to_string(mismatch) + " (";
	ret += to_string(ratio) + "%)\n";
	ratio = ((float)degen)/((float)p_total) * 100.0;
	ret += "\t  Degeneracy: " + to_string(degen) + " (";
	ret += to_string(ratio) + "%)\n";
	ratio = ((float)del)/((float)p_total) * 100.0;
	ret += "\t   Deletions: " + to_string(del) + " (";
	ret += to_string(ratio) + "%)\n";
	ratio = ((float)gc)/((float)p_total) * 100.0;
	ret += "\t  GC Content: " + to_string(gc) + " (";
	ret += to_string(ratio) + "%)\n";
	ratio = ((float)rep)/((float)p_total) * 100.0;
	ret += "\t Repetitions: " + to_string(rep) + " (";
	ret += to_string(ratio) + "%)\n";
	ratio = ((float)comp_end)/((float)p_total) * 100.0;
	ret += "\t  Comp. Ends: " + to_string(comp_end) + " (";
	ret += to_string(ratio) + "%)\n";
	ratio = ((float)hair)/((float)p_total) * 100.0;
	ret += "\t    Hairpins: " + to_string(hair) + " (";
	ret += to_string(ratio) + "%)\n";
	ratio = ((float)dimer)/((float)p_total) * 100.0;
	ret += "\t      Dimers: " + to_string(dimer) + " (";
	ret += to_string(ratio) + "%)\n";
	ratio = ((float)temp)/((float)p_total) * 100.0;
	ret += "\t Temperature: " + to_string(temp) + " (";
	ret += to_string(ratio) + "%)\n\n";

	return ret;
}