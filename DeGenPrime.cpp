#include <cstdlib>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <iterator>
#include <vector>
#include "DataNode.h"
#include "DataSequence.h"
#include "Sequence.h"
#include "SequenceList.h"
#include "SequenceReader.h"

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

	cout << "After building, the size of the list is: " << list.size() << endl;
	list.PrintSequenceNames();
	
	DataSequence data = list.ProcessList();
	data.PrintNonNSequences();

	std::vector<int>indeces = data.IndecesOfNonNSubsequences(20);
	cout << "These Are Indeces where a primer of length 20 would have no 'N's\n";
	for(int i = 0;i < indeces.size();i++)
	{
		cout << "Index [" << indeces[i] << "]\t";
	}
	cout << "\nTotal number of indeces: " << indeces.size() << endl;

	for(int i = 0;i < indeces.size();i++)
	{
		DataSequence sub;
		sub.SetList(data.SubSeq(indeces[i], 20));
		sub.Print(indeces[i]);
	}

	ifs.close();
	return 0;
}