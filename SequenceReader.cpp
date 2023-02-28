#include <fstream>
#include <vector>
#include <string>

#include "SequenceReader.h"
#include "Sequence.h"
#include "SequenceList.h"

using namespace std;

namespace DeGenPrime
{
	
	SequenceReader::SequenceReader() { }

	SequenceList SequenceReader::CreateList(ifstream &ifs)
	{
		FindEndOfHeader(ifs);
		SequenceList list;
		string line;
		while(getline(ifs, line))
		{
			if(ParseLine(line))
			{
				Sequence seq = CreateSequence(line);
				list.PushBack(seq);
			}
		}
		list.PopBack(); // The last line needs to be removed.
		return list;
	}

	Sequence SequenceReader::CreateSequence(string line)
	{
		Sequence seq;
		string name = line.substr(0, 16);
		string codes = line.erase(0, 16);
		seq.SetName(name);
		std::vector<char> list;
		for(int i=0; i<codes.length(); i++)
		{
			list.push_back(codes[i]);
		}
		seq.SetList(list);
		return seq;
	}

	bool SequenceReader::ParseLine(string line)
	{
		bool _ret = true;
		char first = (char)line[0];
		if(first == '\n' || first == ' ' || first == '\t')
		{
			_ret = false;
		}
		return _ret;
	}

	void SequenceReader::FindEndOfHeader(ifstream &ifs)
	{
		int _count = 0;
		char next;
		ifs.get(next);
		while(next != '\n')
		{
			_count++;
			ifs.get(next);
		}
		while(next == '\n')
		{
			_count++;
			ifs.get(next);
		}
		_start_position = _count;
		ReturnToBeginning(ifs);
		return;
	}

	void SequenceReader::ReturnToBeginning(ifstream &ifs)
	{
		ifs.seekg(_start_position, ios::beg);
	}

} // End of DeGenPrime