#include <cctype>
#include <iostream>
#include <string>
#include <vector>
#include "DataNode.h"
#include "DataSequence.h"
#include "SequenceList.h"
#include "Sequence.h"
#include "GlobalSettings.h"
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

		// Loop through the length of entire sequence of chars
		int SizeOfChars = _list[0].GetCodes().size();
		for(int i = 0;i < SizeOfChars;i++)
		{
			std::vector<char> char_list = CharsAt(i);
			DataNode node(char_list);
			data_seq.PushBack(node);
		}
		return data_seq;
	}
	
	void SequenceList::SetList(vector<Sequence> catalog) {_list = catalog; }
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
	SequenceList SequenceList::FilterDashes()
	{
		SequenceList ret;
		int begin = GlobalSettings::GetBeginningNucleotide();
		int ending = GlobalSettings::GetEndingNucleotide();
		if(_list.size() < 2)return _list;
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
		for(Sequence seq : _list)
		{
			for(int i = seq.size() - 1;i >= 0;i--)
			{
				if(seq.GetCodes()[i] == '-')
				{
					seq.Erase(i);
				}
			}
		}
	}
	string SequenceList::PrintSequenceNames() const
	{
		string ret = "";
		for(int i = 0;i < _list.size();i++)
		{
			ret += "Sequence(" + to_string(i);
			ret += "):[";
			ret += _list[i].GetName() + "] Size:[";
			ret += to_string(_list[i].size());
			ret += "]\n";
		}
		return ret;
	}

	string SequenceList::DecodeProteins() const
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
		}
		return ret;
	}

	void SequenceList::PrintSection(int index, int length) const
	{
		cout << "Sequence Data from " << index << " to ";
		cout << index + length - 1 << endl;
		for(int i = 0;i < _list.size();i++)
		{
			Sequence seq = GetSequenceList()[i];
			string name = seq.GetName();
			cout << name;
			for(int j = index;j < index + length;j++)
			{
				cout << _list[i].GetCodes()[j];
			}
		}
		cout << endl;
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

	std::string SequenceList::Codon(char c) const
	{
		string ret = "";
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
		return ret;
	}
} // End of DeGenPrime