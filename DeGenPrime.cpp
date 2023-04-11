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

using namespace std;
using namespace DeGenPrime;

int main(int argc, char *argv[])
{
	
	if(argc != 2)
	{
		cout << "Syntax: ./DeGenPrime <filename>\n";
		return 0;
	}

	
	string filename = argv[1];
	string argument2 = std::filesystem::current_path();
	/*
	// File Information
	cout << "Full path of argument: " << argument2 << endl;
	cout << "Full path of file: " << argument2 << "/" << filename << endl;
	*/

	ifstream ifs;
	ifs.open(argument2 + "/" + filename);
	if(ifs.fail())
	{
		cout << "Failure to open file.\n";
		exit(1);
	}
	// cout << "File opened." << endl;

	SequenceReader read;
	SequenceList list = read.CreateList(ifs);

	/*
	// Sequence Filter Information
	cout << "Before Filter, the number of Sequences in the list is: " << list.size() << endl;
	list.PrintSequenceNames();
	list.FilterDashes();
	
	cout << "After Filter, the number of Sequences in the list is: " << list.size() << endl;
	list.PrintSequenceNames();
	cout << endl;
	*/

	SequenceList reverse_list = list.InvRevList();

	DataSequence data = list.ProcessList();
	DataSequence rev = reverse_list.ProcessList();
	
	/*
	// Test Annealing
	int anneal_index = 700;
	DataSequence anneal = data.SubSeq(anneal_index,20);
	anneal.Print(anneal_index,25.0,50,50,50);
	cout << "Annealing Temp of [" << anneal_index << "] is: " << anneal.BasicAnneal(data, anneal_index) << endl;
	*/
	// TEST ENTROPY, ENTHALPY, GIBBS

	
	//1st TEST (Sequence)
	Sequence test_seq;
	test_seq.SetName("Test Sequence");
	test_seq.PushBack("CGTTGA");

	SequenceList test_list;
	test_list.PushBack(test_seq);

	DataSequence testData = test_list.ProcessList();
	testData.Print(0, 37.0, 50, 50, 50);
	
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
	/*
	// Primer Calculator Section
	PrimerCalculator calc;
	calc.InitializePrimers(data, 2000);
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
	*/
	/*
	PrimerCalculator rev_calc;
	rev_calc.InitializePrimers(rev, 2000);
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
	*/
	/*
	PrimerPairList pairlist;
	pairlist.CreateList(data, rev, calc.GetPrimers(), rev_calc.GetPrimers());
	cout << "Before filters, ";
	pairlist.PrintSize();
	filterchange = 0;

	PrimerPairList testPairList(data, rev, pairlist.GetPairs());

	pairlist.FilterAmpliconLength();
	cout << "After FilterAmpliconLength: Filtered [";
	cout << "testPairList.size() - pairlist.size() - filterchange << "] Primer Pairs. (";
	cout << 100.0 * (float)(testPairList.size() - pairlist.size())/(float)(testPairList.size());
	cout << "% of total filtered)" << endl;
	filterchange = testPairList.size() p pairlist.size();
	
	pairlist.FilterTemperatureDifference();
	cout << "After FilterTemperatureDifference: Filtered [";
	cout << testPairList.size() - pairlist.size() - filterchange << "] Primer Pairs. (";
	cout << 100.0 * (float)(testPairList.size() - pairlist.size())/(float)(testPairList.size());
	cout << "% of total filtered)" << endl;
	filterchange = testPairList.size() - pairlist.size();

	pairlist.FilterCrossDimers();
	cout << "After FilterCrossDimers: Filtered [";
	cout << testPairList.size() - pairlist.size() - filterchange << "] Primer Pairs. (";
	cout << 100.0 * (float)(testPairList.size() - pairlist.size())/(float)(testPairList.size());
	cout << "% of total filtered)" << endl;
	filterchange = testPairList.size() - pairlist.size();

	cout << "After filters, ";
	pairlist.PrintSize();
	*/

	ifs.close();
	return 0;
}