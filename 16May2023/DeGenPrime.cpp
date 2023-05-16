#include <cstdlib>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
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
// #include "GlobalSettings.h"

using namespace std;
using namespace DeGenPrime;

int main(int argc, char *argv[])
{
	
	if(argc == 1)
	{
		cout << "Syntax: ./DeGenPrime <filename>\n";
		exit(PROGRAM_SUCCESS);
	}
	
	string filename = argv[argc - 1];
	string filepath = std::filesystem::current_path();
	
	// File Information
	cout << "Full path of argument: " << filepath << endl;
	cout << "Full path of file: " << filepath << "/" << filename << endl;
	

	ifstream is;
	is.open(filepath + "/SETTINGS.csv");
	if(is.fail())
	{
		cout << "Failure to open settings file.\n";
		exit(SETTINGS_FILE_NOT_FOUND);
	}
	cout << "Opened Settings File." << endl;
	is.close();

	ifstream ifs;
	ifs.open(filepath + "/" + filename);
	if(ifs.fail())
	{
		cout << "Failure to open file.\n";
		exit(BAD_INPUT_FILE);
	}
	cout << "Input File opened." << endl;

	SequenceReader read;
	SequenceList list = read.CreateList(ifs);

	if(list.size() == 0)
	{
		cout << "There were no sequences in the input file.\n";
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
		command1 += "--maxiterate 1000 --clustalout " + filename + " >" + clustal;
		string command2 = "./DeGenPrime " + clustal;
		const char *c1 = command1.c_str();
		const char *c2 = command2.c_str();
		system(c1);
		system(c2);
		exit(FILE_MISALIGNED);
	}
	
	// Sequence Filter Information
	// cout << "Before Filter, the number of Sequences in the list is: " << list.size() << endl;
	list.PrintSequenceNames();
	list.FilterDashes();
	
	// cout << "After Filter, the number of Sequences in the list is: " << list.size() << endl;
	// list.PrintSequenceNames();
	// cout << endl;
	

	SequenceList reverse_list = list.InvRevList();

	DataSequence data = list.ProcessList();
	DataSequence rev = reverse_list.ProcessList();
	
	/*
	// Test Annealing
	int anneal_index = 700;
	int product_length = 200;
	DataSequence anneal = data.SubSeq(anneal_index,20);
	DataSequence product = data.SubSeq(anneal_index, product_length);
	anneal.Print(anneal_index,37.0,50,50);
	cout << "Basic Anneal for amplicon length " << product_length << ": " << anneal.BasicAnneal(product, 50, 50) << endl;
	*/
	// TEST ENTROPY, ENTHALPY, GIBBS

	/*
	//1st TEST (Sequence)
	Sequence test_seq;
	test_seq.SetName("Test Sequence");
	test_seq.PushBack("CAGCGCCACATACATCAT");

	SequenceList test_list;
	test_list.PushBack(test_seq);

	DataSequence testData = test_list.ProcessList();
	testData.Print(0, 37.0, 50, .250);
	*/
	/*
	// 2nd TEST (3' Ends)
	Sequence tester_seq;
	tester_seq.SetName("Test Sequence");
	tester_seq.PushBack("ATCGC");

	SequenceList tester_list;
	tester_list.PushBack(tester_seq);

	DataSequence testerData = tester_list.ProcessList();
	testerData.Print(0, 25.0, 50, 50, 50);
	*/

	// Primer Calculator Section
	PrimerCalculator calc;
	calc.InitializePrimers(data, 0);
	int filterchange = 0;

	cout << "Number of possible forward primers in specified size range, before Filters: [" << calc.size() << "]\n" << endl;
	PrimerCalculator test_calc;
	test_calc.SetPrimers(calc.GetPrimers());
	
	calc.FilterDegeneracy(data);
	cout << "After FilterDegeneracy: Filtered [";
	cout << test_calc.size() - calc.size() - filterchange << "] Primers. (" << 100.0 * (float)(test_calc.size() - calc.size())/(float)(test_calc.size());
	cout << "% of total filtered)" << endl;
	filterchange = test_calc.size() - calc.size();

	calc.FilterDeletions(data, list);
	cout << "After FilterDeletions: Filtered [";
	cout << test_calc.size() - calc.size() - filterchange << "] Primers. (" << 100.0 * (float)(test_calc.size() - calc.size())/(float)(test_calc.size());
	cout << "% of total filtered)" << endl;
	filterchange = test_calc.size() - calc.size();

	calc.FilterGCContent(data);
	cout << "After FilterGCContent: Filtered [";
	cout << test_calc.size() - calc.size() - filterchange << "] Primers. (" << 100.0 * (float)(test_calc.size() - calc.size())/(float)(test_calc.size());
	cout << "% of total filtered)" << endl;
	filterchange = test_calc.size() - calc.size();

	calc.FilterComplementaryEnds(data);
	cout << "After FilterComplementaryEnds: Filtered [";
	cout << test_calc.size() - calc.size() - filterchange << "] Primers. (" << 100.0 * (float)(test_calc.size() - calc.size())/(float)(test_calc.size());
	cout << "% of total filtered)" << endl;
	filterchange = test_calc.size() - calc.size();

	calc.FilterHairpins(data);
	cout << "After FilterHairpins: Filtered [";
	cout << test_calc.size() - calc.size() - filterchange << "] Primers. (" << 100.0 * (float)(test_calc.size() - calc.size())/(float)(test_calc.size());
	cout << "% of total filtered)" << endl;
	filterchange = test_calc.size() - calc.size();

	
	calc.FilterDimers(data);
	cout << "After FilterDimers: Filtered [";
	cout << test_calc.size() - calc.size() - filterchange << "] Primers. (" << 100.0 * (float)(test_calc.size() - calc.size())/(float)(test_calc.size());
	cout << "% of total filtered)" << endl;
	filterchange = test_calc.size() - calc.size();
	
	
	calc.FilterRepeats(data);
	cout << "After FilterRepeats: Filtered [";
	cout << test_calc.size() - calc.size() - filterchange << "] Primers. (" << 100.0 * (float)(test_calc.size() - calc.size())/(float)(test_calc.size());
	cout << "% of total filtered)" << endl;
	filterchange = test_calc.size() - calc.size();

	int offset = 0;
	do
	{
		float tenthOfDifference = (float)(MAX_PRIMER_TEMP - MIN_PRIMER_TEMP)/ 10.0;
		calc.FilterTemperature(data, offset * tenthOfDifference);
		cout << "After FilterTemperature[" << offset << "]: Filtered [";
		cout << test_calc.size() - calc.size() - filterchange << "] Primers. (" << 100.0 * (float)(test_calc.size() - calc.size())/(float)(test_calc.size());
		cout << "% of total filtered)" << endl;
		filterchange = test_calc.size() - calc.size();
		offset += 1;
	} while(calc.size() > 400 && offset < 10);
	
	calc.PrintSize();
	cout << endl;


	PrimerCalculator rev_calc;
	rev_calc.InitializePrimers(rev, 0);
	filterchange = 0;
	
	cout << "Number of possible reverse primers in specified size range, before Filters: [" << rev_calc.size() << "]\n" << endl;
	test_calc.SetPrimers(rev_calc.GetPrimers());
	
	rev_calc.FilterDegeneracy(rev);
	cout << "After FilterDegeneracy: Filtered [";
	cout << test_calc.size() - rev_calc.size() - filterchange << "] Primers. (" << 100.0 * (float)(test_calc.size() - rev_calc.size())/(float)(test_calc.size());
	cout << "% of total filtered)" << endl;
	filterchange = test_calc.size() - rev_calc.size();

	rev_calc.FilterDeletions(rev, reverse_list);
	cout << "After FilterDeletions: Filtered [";
	cout << test_calc.size() - rev_calc.size() - filterchange << "] Primers. (" << 100.0 * (float)(test_calc.size() - rev_calc.size())/(float)(test_calc.size());
	cout << "% of total filtered)" << endl;
	filterchange = test_calc.size() - rev_calc.size();

	rev_calc.FilterGCContent(rev);
	cout << "After FilterGCContent: Filtered [";
	cout << test_calc.size() - rev_calc.size() - filterchange << "] Primers. (" << 100.0 * (float)(test_calc.size() - rev_calc.size())/(float)(test_calc.size());
	cout << "% of total filtered)" << endl;
	filterchange = test_calc.size() - rev_calc.size();

	rev_calc.FilterComplementaryEnds(rev);
	cout << "After FilterComplementaryEnds: Filtered [";
	cout << test_calc.size() - rev_calc.size() - filterchange << "] Primers. (" << 100.0 * (float)(test_calc.size() - rev_calc.size())/(float)(test_calc.size());
	cout << "% of total filtered)" << endl;
	filterchange = test_calc.size() - rev_calc.size();

	rev_calc.FilterHairpins(rev);
	cout << "After FilterHairpins: Filtered [";
	cout << test_calc.size() - rev_calc.size() - filterchange << "] Primers. (" << 100.0 * (float)(test_calc.size() - rev_calc.size())/(float)(test_calc.size());
	cout << "% of total filtered)" << endl;
	filterchange = test_calc.size() - rev_calc.size();

	
	rev_calc.FilterDimers(rev);
	cout << "After FilterDimers: Filtered [";
	cout << test_calc.size() - rev_calc.size() - filterchange << "] Primers. (" << 100.0 * (float)(test_calc.size() - rev_calc.size())/(float)(test_calc.size());
	cout << "% of total filtered)" << endl;
	filterchange = test_calc.size() - rev_calc.size();
	

	rev_calc.FilterRepeats(rev);
	cout << "After FilterRepeats: Filtered [";
	cout << test_calc.size() - rev_calc.size() - filterchange << "] Primers. (" << 100.0 * (float)(test_calc.size() - rev_calc.size())/(float)(test_calc.size());
	cout << "% of total filtered)" << endl;
	filterchange = test_calc.size() - rev_calc.size();

	offset = 0;
	do
	{
		float tenthOfDifference = (float)(MAX_PRIMER_TEMP - MIN_PRIMER_TEMP) / 10.0;
		rev_calc.FilterTemperature(rev, offset * tenthOfDifference);
		cout << "After FilterTemperature[" << offset << "]: Filtered [";
		cout << test_calc.size() - rev_calc.size() - filterchange << "] Primers. (";
		cout << 100.0 * (float)(test_calc.size() - rev_calc.size())/(float)(test_calc.size());
		cout << "% of total filtered)" << endl;
		filterchange = test_calc.size() - rev_calc.size();
		offset++;
	} while(rev_calc.size() > 400 && offset < 10);

	rev_calc.PrintSize();
	cout << endl;
	
	
	PrimerPairList pairlist;
	pairlist.CreateList(data, rev, calc.GetPrimers(), rev_calc.GetPrimers());
	cout << "Before filters, ";
	pairlist.PrintSize();
	filterchange = 0;

	PrimerPairList testPairList(data, rev, pairlist.GetPairs());

	pairlist.FilterAmpliconLength();
	cout << "After FilterAmpliconLength: Filtered [";
	cout << testPairList.size() - pairlist.size() - filterchange << "] Primer Pairs. (";
	cout << 100.0 * (float)(testPairList.size() - pairlist.size())/(float)(testPairList.size());
	cout << "% of total filtered)" << endl;
	filterchange = testPairList.size() - pairlist.size();
	
	pairlist.FilterTemperatureDifference();
	cout << "After FilterTemperatureDifference: Filtered [";
	cout << testPairList.size() - pairlist.size() - filterchange << "] Primer Pairs. (";
	cout << 100.0 * (float)(testPairList.size() - pairlist.size())/(float)(testPairList.size());
	cout << "% of total filtered)" << endl;
	filterchange = testPairList.size() - pairlist.size();

	cout << "After filters, ";
	pairlist.PrintSize();

	cout << "\nSorting remaining " << pairlist.size();
	cout << " primer pairs in the list by temperature difference." << endl;
	pairlist.Sort();
	cout << "Finished Sorting!" << endl;

	int desiredpairs = 5;
	cout << "Running FilterAnnealingTemp on the top " << desiredpairs << " primer pairs." << endl;
	PrimerPairList top = pairlist.SubList(0,desiredpairs);
	int nextIndex = desiredpairs;
	int filtercount = top.FilterAnnealingTemp(data, rev);
	while(filtercount != 0)
	{
		PrimerPairList next = pairlist.SubList(nextIndex,filtercount);
		filtercount = next.FilterAnnealingTemp(data, rev);
		top.Append(next);
	}
	
	top.PrintAll(data, rev);
	

	ifs.close();
	return 0;
}
