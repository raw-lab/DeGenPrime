#include <iostream>
#include <vector>
#include "DataNode.h"
#include "DataSequence.h"
#include "SequenceList.h"
#include "Sequence.h"
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
	void SequenceList::FilterDashes()
	{
		for(int i = _list.size() - 1;i > -1;i--)
		{
			int dash_count = 0;
			for(int j = 0;j < _list[i].GetCodes().size();j++)
			{
				if(_list[i].GetCodes()[j] == '-')
				{
					dash_count++;
				}
			}
			float tolerance = (float)dash_count / (float)_list[i].GetCodes().size();
			if(tolerance >= MAX_DASH_HORIZONTAL_TOLERANCE)
			{
				Erase(i);
			}
		}
		return;
	}
	void SequenceList::PrintSequenceNames() const
	{
		for(int i = 0;i < _list.size();i++)
		{
			cout << "name=" << _list[i].GetName() << "*";
			cout << " size =" << _list[i].size() << endl;
		}
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
} // End of DeGenPrime