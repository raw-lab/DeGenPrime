// DataSequence.cpp
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include "DataNode.h"
#include "DataSequence.h"

using namespace std;

namespace DeGenPrime
{
	DataSequence::DataSequence() { }

	void DataSequence::SetList(std::vector<DataNode> catalog) {_list = catalog; }
	void DataSequence::PushBack(DataNode node) { _list.push_back(node); }
	void DataSequence::PopBack() { _list.pop_back(); }

	void DataSequence::Print(int id)
	{
		double accum_ratio = 0.0;
		double accum_weight = 0.0;
		cout << "DataSequence ID [" << id;
		cout << "]\tLength: [" << _list.size() << "]" << endl;
		cout << "Sequence Codes: ";
		for(int i = 0;i < _list.size();i++)
		{
			cout << _list[i].GetCode();
			/* Used to calculate standard deviation
			accum_ratio += (_list[i].Ratio() - AverageRatio()) *
						(_list[i].Ratio() - AverageRatio());
			accum_weight += (_list[i].WeightedRatio() - AverageWeighted()) *
						(_list[i].WeightedRatio() - AverageWeighted());
			*/
		}
		cout << "\nSequence Most Commons: ";
		for(int i = 0;i < _list.size();i++)
		{
			cout << _list[i].GetMostCommon();
		}
		cout << "\nSequence Average Match Ratio: " << AverageRatio() << endl;
		// accum_percent = sqrt(accum_ratio);
		// accum_weight = sqrt(accum_weight);
		// cout << "\t(std: " << accum_percent << "%)" << endl;
		cout << "Sequence Average Weighted Ratio: " << AverageWeighted();
		// cout << " %(std: " << accum_weight << "%)\n" << endl;
		cout << endl << endl;
		return;
	}
	void DataSequence::PrintNonNSequences()
	{
		int counter = 0;
		float ratio = 0.0;
		float weighted = 0.0;
		for(int i = 0;i < _list.size();i++)
		{
			if(_list[i].GetCode() == 'N')
			{
				continue;
			}
			else
			{
				cout << "[Index: " << i << "] ";
				_list[i].Print();
				cout << endl;
				counter++;
				ratio += _list[i].Ratio();
				weighted += _list[i].WeightedRatio();
			}
		}
		ratio /= counter;
		weighted /= counter;
		cout << "Total Non-N Sequences: [" << counter << "]" << endl;
		cout << "Average Match Score: [" << ratio << endl;
		cout << "Average Weighted Score: [" << weighted << endl;
	}

	std::vector<DataNode> DataSequence::SubSeq(int startIndex, int length)
	{
		std::vector<DataNode> sub_seq;
		for(int i = startIndex;i < startIndex+length;i++)
		{
				sub_seq.push_back(_list[i]);
		}
		return sub_seq;

	}

	std::vector<int> DataSequence::IndecesOfNonNSubsequences(int length)
	{
		std::vector<int> indeces;
		int ending_index = _list.size() - length;
		for(int i = 0;i < ending_index; ) // Third argument intentionally blank.
		{
			bool was_no_n = true;
			int next_index = i + 1;
			// Base nested for loop on i and length
			for(int j = i;j < i + length;j++)
			{
				if(_list[j].GetCode() == 'N')
				{
					next_index = j + 1;
					was_no_n = false;
				}
			}
			if(was_no_n)
			{
				indeces.push_back(i);
			}
			i = next_index;
		}
		return indeces;
	}

	std::vector<DataNode> DataSequence::GetDataSequence() const { return _list; }
	int DataSequence::size() const { return _list.size(); }
	float DataSequence::AverageRatio() const
	{
		int counter = 0;
		float accum = 0.0;
		for(int i = 0;i < _list.size();i++)
		{
			accum += _list[i].Ratio();
			counter++;
		}
		accum /= counter;
		return accum;
	}
	float DataSequence::AverageWeighted() const
	{
		int counter = 0;
		float accum = 0.0;
		for(int i = 0;i < _list.size();i++)
		{
			accum += _list[i].WeightedRatio();
			counter++;
		}
		accum /= counter;
		return accum;
	}
} // End of DeGenPrime