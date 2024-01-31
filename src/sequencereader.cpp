#include <cctype>
#include <fstream>
#include <vector>
#include <string>

#include "sequencereader.h"
#include "sequence.h"
#include "sequencelist.h"
#include "global.h"

using namespace std;

namespace DeGenPrime
{
	
	SequenceReader::SequenceReader() { }

	FileType SequenceReader::IdentifyFileType(ifstream &ifs)
	{
		FileType ret;
		int forward_arrow = (int)'>';
		if(ifs.peek() == forward_arrow)
		{
			ret = fasta;
		}
		else
		{
			ret = clust;
		}
		return ret;
	}

	SequenceList SequenceReader::CreateList(ifstream &ifs)
	{
		SequenceList list;
		Sequence seq;
		FileType type = IdentifyFileType(ifs);
		string line;
		int counter = 1;
		char next;
		switch(type)
		{
			case clust:
				FindEndOfHeader(ifs);
				while(getline(ifs, line))
				{
					if(ParseLine(line))
					{
						seq = CreateSequence(line);
						list.PushBack(seq);
					}
				}
				if(list.size() > 1)
				{
					list.PopBack(); // The last sequence needs to be removed.
				}
				break;
			case fasta:
				while(getline(ifs,line))
				{
					counter++;
					if(line[0] == '>') // New Sequence
					{
						list.PushBack(seq); // The first sequence will be empty
						line = line.erase(0,1); // erase the '>'
						seq.SetName(line);
					}
					else // Adding to existing sequence
					{
						for(int i = 0;i < line.length();i++)
						{
							seq.PushBack(line[i]);
						}
					}
				}
				// Push the last collected sequence and delete the empty first sequence
				list.PushBack(seq);
				list.Erase(0);
				break;
			default:
				break;
		}
		return list;
	}

	Sequence SequenceReader::CreateSequence(string line)
	{
		Sequence seq;
		string line2 = "";
		for(char c : line)
		{
			if(c != '\r')
			{
				line2 += c;
			}
		}
		string name = line2.substr(0, 16);
		string codes = line2.erase(0, 16);
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
		return ( std::isblank((char)line[0]) == 1 ? false : true);
	}

	void SequenceReader::FindEndOfHeader(ifstream &ifs)
	{
		int _count = 0;
		char next;
		ifs.get(next);
		while(next != '\n' && next != '\r')
		{
			_count++;
			ifs.get(next);
		}
		while(next == '\n' || next == '\r')
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