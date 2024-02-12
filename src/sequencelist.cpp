#include <algorithm>
#include <cctype>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include "datanode.h"
#include "datasequence.h"
#include "sequencelist.h"
#include "sequence.h"
#include "globalsettings.h"
#include "format.h"
#include "global.h"

using namespace std;

namespace DeGenPrime
{
	SequenceList::SequenceList() { }
	
	SequenceList SequenceList::InvRevList()
	{
		SequenceList rev_inv_list;
		for(int i = 0;i < _list.size();i++)
		{
			Sequence seq = _list[i];
			seq.Invert();
			seq.Reverse();
			rev_inv_list.PushBack(seq);
		}
		return rev_inv_list;
	}

	DataSequence SequenceList::ProcessList()
	{
		DataSequence data_seq;

		int SizeOfChars = _list[0].GetCodes().size();
		// If there is only one sequence on the list, use constructors
		if(size() == 1)
		{
			std::vector<char> char_list = _list[0].GetCodes();
			for(int i = 0; i < SizeOfChars;i++)
			{
				char c = char_list[i];
				if(c == '\n' || c == ' '  || c == '\r')
				{
					continue;
				}
				DataNode node(char_list[i], char_list[i], 1.0);
				data_seq.PushBack(node);
			}
		}
		else
		{
			// Loop through the length of entire sequence of chars
			for(int i = 0;i < SizeOfChars;i++)
			{
				std::vector<char> char_list = CharsAt(i);
				DataNode node(char_list);
				data_seq.PushBack(node);
			}
		}
		return data_seq;
	}
	
	void SequenceList::SetList(vector<Sequence> catalog) 
	{
		_list.clear(); 
		for(Sequence seq : catalog)
		{
			_list.push_back(seq);
		}
	}
	void SequenceList::Clear()
	{
		_list.clear();
	}
	void SequenceList::Erase(int index)
	{
		_list.erase(_list.begin() + index);
	}
	void SequenceList::PushBack(Sequence seq)
	{
		std::vector<Sequence>::iterator iter;
		for(iter=_list.begin();iter != _list.end();iter++)
		{
			if(iter->GetName() == seq.GetName())
			{
				for(char c : seq.GetCodes())
				{
					iter->PushBack(c);
				}
				return;
			}
		}
		_list.push_back(seq);
		return;
	}
	void SequenceList::PopBack() { _list.pop_back(); }
	void SequenceList::Sort()
	{
		sort(_list.begin(), _list.end());
	}
	SequenceList SequenceList::FilterDashes()
	{
		SequenceList ret;
		int begin = GlobalSettings::GetBeginningNucleotide();
		int ending = GlobalSettings::GetEndingNucleotide();
		for(int i = 0;i < _list.size();i++)
		{
			int dash_count = 0;
			for(int j = begin;j < ending;j++)
			{
				if(_list[i].GetCodes()[j] == '-')
				{
					dash_count++;
				}
			}
			float tolerance = (float)dash_count / (float)_list[i].GetCodes().size();
			if(tolerance >= MAX_DASH_HORIZONTAL_TOLERANCE)
			{
				ret.PushBack(_list[i]);
				Erase(i);
			}
		}
		return ret;
	}

	void SequenceList::RemoveDashes()
	{
		std::vector<Sequence> seq_temp;
		std::vector<char> chars;
		for(int i = 0;i < _list.size();i++)
		{
			Sequence seq(_list[i].GetName());
			for(int j = 0;j < _list[i].GetCodes().size();j++)
			{
				char c = _list[i].GetCodes()[j];
				if(c != '-')
				{
					chars.push_back(c);
				}
			}
			seq.SetList(chars);
			chars.clear();
			seq_temp.push_back(seq);
		}
		SetList(seq_temp);
	}

	string SequenceList::PrintSequenceNames() const
	{
		string ret = "";
		string line;
		for(int i = 0;i < _list.size();i++)
		{
			line = "";
			string seq_num = Format(i + 1,4);
			line += "Sequence(" + seq_num + "):[";
			line += Format(_list[i].GetName(), 16, Alignment::Right);
			line += "] Size:[" + to_string(_list[i].size()); // All should be same size
			line += "]";
			ret += Format(line, STR_FORMAT, Alignment::Center);
			ret += "\n";
		}
		return ret;
	}

	std::string SequenceList::CreateFasta() const
	{
		string fasta = "";
		for(Sequence seq : _list)
		{
			fasta += seq.Fasta();
		}
		return fasta;
	}

	std::string SequenceList::DecodeProteins() const
	{
		string ret = "";
		int count = 0;
		for(Sequence seq : _list)
		{
			ret += ">" + seq.GetName() + "\n";
			for(char c : seq.GetCodes())
			{
				ret += Codon(c);
				count += 3;
				if(count == 69)
				{
					ret += "\n";
					count = 0;
				}
			}
			ret += "\n";
		}
		return ret;
	}

	std::string SequenceList::Section(int index, int length) const
	{
		string ret = "";
		for(int i = 0;i < _list.size();i++)
		{
			Sequence seq = GetSequenceList()[i];
			ret += seq.GetName();
			for(int j = index;j < index + length;j++)
			{
				ret += string(1, _list[i].GetCodes()[j]);
			}
			ret += "\n";
		}
		return ret;
	}

	bool SequenceList::TestAlignment() const
	{
		bool aligned = true;
		std::vector<int> sizes;
		if(_list.size() == 1)
		{
			aligned = true;
		}
		else
		{
			for(Sequence seq : _list)
			{
				sizes.push_back(seq.size());
			}
			for(int i = 1;i < _list.size();i++)
			{
				if(aligned == false)
				{
					break;
				}
				if(sizes[0] != sizes[i])
				{
					aligned = false;
				}
			}
		}
		return aligned;
	}

	std::vector<char> SequenceList::CharsAt(int index) const
	{
		std::vector<char> char_list;
		for(int i = 0;i < _list.size();i++)
		{
			char_list.push_back(_list[i].GetCodes()[index]);
		}
		return char_list;
	}

	std::vector<Sequence> SequenceList::GetSequenceList() const { return _list; }
	int SequenceList::size() const { return _list.size(); }
	int SequenceList::IndexOf(std::string name)
	{
		int index = -1;
		for(int i = 0;i < _list.size();i++)
		{
			string seq_name = _list[i].GetName();
			seq_name = seq_name.substr(0,name.length());
			if(name.compare(seq_name) == 0)
			{
				index = i;
				break;
			}
		}
		return index;
	}

	std::string SequenceList::Codon(char c) const
	{
		std::map<char, std::string> str = 
		{
			{'A', "GCT"}, {'R', "CGT"}, {'N', "AAT"}, {'D', "GAT"},
    		{'C', "TGT"}, {'Q', "CAA"}, {'E', "GAA"}, {'G', "GGT"},
    		{'H', "CAT"}, {'I', "ATT"}, {'L', "CTT"}, {'K', "AAA"},
    		{'M', "ATG"}, {'F', "TTT"}, {'P', "CCT"}, {'S', "TCT"},
    		{'T', "ACT"}, {'W', "TGG"}, {'Y', "TAT"}, {'V', "GTT"},
    		{'*', "TAA"} // Stop codon
		};

		return str[c];
		/*
		switch(toupper(c))
		{
			case 'A': // Alanine
				ret += "GCN";
				break;
			case 'R': // Arginine
				ret += "MGN";
				break;
			case 'D': // Aspartic acid
				ret += "GAY";
				break;
			case 'N': // Asparagine
				ret += "AAY";
				break;
			case 'C': // Cysteine
				ret += "TGY";
				break;
			case 'E': // Glutamic acid
				ret += "GAR";
				break;
			case 'Q': // Glutamine
				ret += "CAR";
				break;
			case 'G': // Glycine
				ret += "GGN";
				break;
			case 'H': // Histidine
				ret += "CAY";
				break;
			case 'I': // Isoleucine
				ret += "ATH";
				break;
			case 'L': // Leucine
				ret += "YTN";
				break;
			case 'K': // Lysine
				ret += "AAR";
				break;
			case 'M': // Methionine, START
				ret += "ATG";
				break;
			case 'F': // Phenylalanine
				ret += "TTY";
				break;
			case 'P': // Proline
				ret += "CCN";
				break;
			case 'S': // Serine
				ret += "WSN";
				break;
			case 'T': // Threonine
				ret += "ACN";
				break;
			case 'W': // Tryptophan
				ret += "TGG";
				break;
			case 'Y': // Tyrosine
				ret += "TAY";
				break;
			case 'V': // Valine
				ret += "GTN";
				break;
			case '*': // Stop Codon
				ret += "TRR";
				break;
			default:
				break;
		}
		*/
	}
} // End of DeGenPrime