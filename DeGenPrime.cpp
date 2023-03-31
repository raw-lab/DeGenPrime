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

	cout << "Full path of argument: " << argument2 << endl;
	cout << "Full path of file: " << argument2 << "/" << filename << endl;

	ifstream ifs;
	ifs.open(argument2 + "/" + filename);
	if(ifs.fail())
	{
		cout << "Failure to open file.\n";
		exit(1);
	}
	cout << "File opened." << endl;

	SequenceReader read;
	SequenceList list = read.CreateList(ifs);
	cout << "Before Filter, the number of Sequences in the list is: " << list.size() << endl;
	list.PrintSequenceNames();
	list.FilterDashes();
	
	cout << "After Filter, the number of Sequences in the list is: " << list.size() << endl;
	list.PrintSequenceNames();

	SequenceList reverse_list = list.InvRevList();

	DataSequence data = list.ProcessList();
	DataSequence rev = reverse_list.ProcessList();

	/*
	std::vector<char> test_chars;
	for(char c : "GATTACAGATTACAGATTACA")
	{
		test_chars.push_back(c);
	}
	test_chars.pop_back();

	Sequence test_seq;
	test_seq.SetName("Test Sequence");
	test_seq.SetList(test_chars);

	SequenceList test_list;
	test_list.PushBack(test_seq);

	DataSequence testData = test_list.ProcessList();
	DataSequence testInvRev = testData.InvSeq().RevSeq();
	testData.Print(0);
	testInvRev.Print(0);
	*/
	/*
	cout << "Printing all 'N' indeces." << endl;
	for(int i = 0;i < 180;i++)
	{
		if(data.GetDataSequence()[i].GetCode() == 'N')
		{
			cout << "N found at index: [" << i << "]" << endl;
			data.GetDataSequence()[i].Print();
			cout << endl;
		}
	}
	*/
	
	PrimerCalculator calc;
	calc.InitializePrimers(data, 0);

	cout << "Number of possible primers in specified size range, before Filters: [" << calc.size() << "]\n" << endl;
	PrimerCalculator test_calc;
	test_calc.SetPrimers(calc.GetPrimers());
	
	calc.FilterDegeneracy(data);
	cout << "After FilterDegeneracy: Filtered [";
	cout << test_calc.size() - calc.size() << "] Primers. (" << 100.0 * (float)(test_calc.size() - calc.size())/(float)(test_calc.size());
	cout << "%)" << endl;

	calc.FilterDeletions(data, list);
	cout << "After FilterDeletions: Filtered [";
	cout << test_calc.size() - calc.size() << "] Primers. (" << 100.0 * (float)(test_calc.size() - calc.size())/(float)(test_calc.size());
	cout << "%)" << endl;

	calc.FilterGCContent(data);
	cout << "After FilterGCContent: Filtered [";
	cout << test_calc.size() - calc.size() << "] Primers. (" << 100.0 * (float)(test_calc.size() - calc.size())/(float)(test_calc.size());
	cout << "%)" << endl;
	
	//calc.FilterGibbs(data, 37.0, 1.0);
	//cout << "After FilterGibbs: Filtered [";
	//cout << test_calc.size() - calc.size() << "] Primers. (" << 100.0 * (float)(test_calc.size() - calc.size())/(float)(test_calc.size());
	//cout << "%)" << endl;

	calc.FilterRepeats(data);
	cout << "After FilterRepeats: Filtered [";
	cout << test_calc.size() - calc.size() << "] Primers. (" << 100.0 * (float)(test_calc.size() - calc.size())/(float)(test_calc.size());
	cout << "%)" << endl;

	calc.FilterComplementaryEnds(data);
	cout << "After FilterComplementaryEnds: Filtered [";
	cout << test_calc.size() - calc.size() << "] Primers. (" << 100.0 * (float)(test_calc.size() - calc.size())/(float)(test_calc.size());
	cout << "%)" << endl;	

	calc.FilterTemperature(data);
	cout << "After FilterTemperature: Filtered [";
	cout << test_calc.size() - calc.size() << "] Primers. (" << 100.0 * (float)(test_calc.size() - calc.size())/(float)(test_calc.size());
	cout << "%)" << endl;
	

	calc.PrintSize();
	//calc.PrintAll();
	
	/*
	DataNode node1('A', 'A', 1.0);
	DataNode node2('C', 'C', 1.0);
	DataNode node3('G', 'G', 1.0);
	DataNode node4('T', 'T', 1.0);

	cout << "Enthalpy of " << node1.GetMostCommon() << " and " << node1.GetMostCommon() << ": " << node1.Enthalpy(node1) << "\t";
	cout << "Enthalpy of " << node1.GetMostCommon() << " and " << node2.GetMostCommon() << ": " << node1.Enthalpy(node2) << "\t";
	cout << "Enthalpy of " << node1.GetMostCommon() << " and " << node3.GetMostCommon() << ": " << node1.Enthalpy(node3) << "\t";
	cout << "Enthalpy of " << node1.GetMostCommon() << " and " << node4.GetMostCommon() << ": " << node1.Enthalpy(node4) << endl;

	cout << "Entropy of " << node1.GetMostCommon() << " and " << node1.GetMostCommon() << ": " << node1.Entropy(node1) << "\t";
	cout << "Entropy of " << node1.GetMostCommon() << " and " << node2.GetMostCommon() << ": " << node1.Entropy(node2) << "\t";
	cout << "Entropy of " << node1.GetMostCommon() << " and " << node3.GetMostCommon() << ": " << node1.Entropy(node3) << "\t";
	cout << "Entropy of " << node1.GetMostCommon() << " and " << node4.GetMostCommon() << ": " << node1.Entropy(node4) << endl;
	
	cout << "Gibbs of " << node1.GetMostCommon() << " and " << node1.GetMostCommon() << ": " << node1.Gibbs(node1, 37.0) << "\t";
	cout << "Gibbs of " << node1.GetMostCommon() << " and " << node2.GetMostCommon() << ": " << node1.Gibbs(node2, 37.0) << "\t";
	cout << "Gibbs of " << node1.GetMostCommon() << " and " << node3.GetMostCommon() << ": " << node1.Gibbs(node3, 37.0) << "\t";
	cout << "Gibbs of " << node1.GetMostCommon() << " and " << node4.GetMostCommon() << ": " << node1.Gibbs(node4, 37.0) << endl;
	*/
	/*
	for(Primer p : calc.GetPrimers())
	{
		DataSequence temp = data.SubSeq(p.Index(), p.Length());
		temp.Print(p.Index());
	}
	*/
	/*
	test_calc.FilterDegeneracy(data);
	cout << "FilterDegeneracy: Filtered [";
	cout << calc.size() - test_calc.size() << "] Primers. (" << 100.0 *(float)(calc.size() - test_calc.size())/(float)(calc.size());
	cout << "%)" << endl;
	// calc.PrintSize();
	
	test_calc.SetPrimers(calc.GetPrimers());
	test_calc.FilterDeletions(data, list);
	cout << "FilterDeletions: Filtered [";
	cout << calc.size() - test_calc.size() << "] Primers. (" << 100.0 * (float)(calc.size() - test_calc.size())/(float)(calc.size());
	cout << "%)" << endl;
	// calc.PrintSize();

	// test_calc.SetPrimers(calc.GetPrimers());
	test_calc.FilterGCContent(data);
	cout << "FilterGCContent: Filtered [";
	cout << calc.size() - test_calc.size() << "] Primers. (" << 100.0 * (float)(calc.size() - test_calc.size())/(float)(calc.size());
	cout << "%)" << endl;
	//calc.PrintSize();
	
	test_calc.SetPrimers(calc.GetPrimers());
	test_calc.FilterRepeats(data);
	cout << "FilterRepeats: Filtered [";
	cout << calc.size() - test_calc.size() << "] Primers. (" << 100.0 * (float)(calc.size() - test_calc.size())/(float)(calc.size());
	cout << "%)" << endl;
	// calc.PrintSize();
	
	test_calc.SetPrimers(calc.GetPrimers());
	test_calc.FilterComplementaryEnds(data);
	cout << "FilterComplementaryEnds: Filtered [";
	cout << calc.size() - test_calc.size() << "] Primers. (" << 100.0 * (float)(calc.size() - test_calc.size())/(float)(calc.size());
	cout << "%)" << endl;
	// calc.PrintSize();

	test_calc.SetPrimers(calc.GetPrimers());
	test_calc.FilterTemperature(data);
	cout << "FilterTemperature: Filtered [";
	cout << calc.size() - test_calc.size() << "] Primers. (" << 100.0 * (float)(calc.size() - test_calc.size())/(float)(calc.size());
	cout << "%)" << endl;
	*/
	/*
	calc.PrintSize();
	calc.PrintAll();
	for(Primer p : calc.GetPrimers())
	{
		DataSequence temp = data.SubSeq(p.Index(), p.Length());
		temp.Print(p.Index());
		cout << "Primer Code: " << endl;
		temp = temp.InvSeq();
		temp.Print(p.Index());
	}

	PrimerCalculator rev_calc;
	rev_calc.InitializePrimers(rev);
	
	cout << "Rev, before Filters: " << endl;
	rev_calc.PrintSize();
	
	rev_calc.FilterDegeneracy(rev);
	cout << "Rev, after FilterDegeneracy: " << endl;
	rev_calc.PrintSize();

	rev_calc.FilterDeletions(rev, reverse_list);
	cout << "Rev, after FilterDeletions: " << endl;
	rev_calc.PrintSize();

	rev_calc.FilterGCContent(rev);
	cout << "Rev, after FilterGCContent: " << endl;
	rev_calc.PrintSize();
	
	// rev_calc.FilterRepeats(rev);
	// cout << "Rev, after FilterRepeats: " << endl;
	// rev_calc.PrintSize();
	
	rev_calc.FilterComplementaryEnds(rev);
	cout << "Rev, after FilterComplementaryEnds: " << endl;
	rev_calc.PrintSize();

	rev_calc.FilterTemperature(rev);
	cout << "Rev, after FilterTemperature: " << endl;
	rev_calc.PrintSize();
	*/
	// rev_calc.PrintAll();

	ifs.close();
	return 0;
}