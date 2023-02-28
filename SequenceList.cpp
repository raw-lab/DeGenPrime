#include <iostream>
#include <vector>
#include "DataNode.h"
#include "DataSequence.h"
#include "SequenceList.h"
#include "Sequence.h"

using namespace std;

namespace DeGenPrime
{
	SequenceList::SequenceList() { }

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
	void SequenceList::PrintSequenceNames()
	{
		std::vector<Sequence>::iterator iter;
		for(iter = _list.begin();iter != _list.end();iter++)
		{
			cout << "name=" << iter->GetName() << "*";
			cout << " size=" << iter->size() << endl;
		}
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