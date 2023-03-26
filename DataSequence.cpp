// DataSequence.cpp
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include "DataNode.h"
#include "DataSequence.h"
#include "global.h"

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
		}
		cout << endl << "Temperature: " << Temperature() << endl;
		cout << "Enthalpy: " << Enthalpy() << endl;
		cout << "Entropy (1M Na+): " << Entropy(1.0) << endl;
		cout << "Gibbs (25 C, 1M Na+): " << Gibbs(25.0, 1.0) << endl;
		return;
	}

	DataSequence DataSequence::SubSeq(int startIndex, int length)
	{
		DataSequence sub_seq;
		for(int i = startIndex;i < startIndex+length;i++)
		{
				sub_seq.PushBack(_list[i]);
		}
		return sub_seq;

	}

	DataSequence DataSequence::InvSeq()
	{
		DataSequence inv_seq;
		for(DataNode node : _list)
		{
			DataNode inv_node = node.InvNode();
			inv_seq.PushBack(inv_node);
		}
		return inv_seq;
	}

	DataSequence DataSequence::RevSeq()
	{
		DataSequence rev_seq;
		for(int i = _list.size() - 1; i > -1;i--)
		{
			rev_seq.PushBack(_list[i]);
		}
		return rev_seq;
	}

	float DataSequence::Enthalpy() const
	{
		float accumulated_enthalpy = 0.0;
		for(int i = 0;i < _list.size() - 1;) // Third argument intentionally blank.
		{
			// Skip deletions in first nucleotide
			if(_list[i].GetMostCommon() == '-')
			{
				continue;
			}
			int j = i + 1;
			
			// Skip deletions in second nucleotide
			while(_list[j].GetMostCommon() == '-' && j < _list.size())
			{
				j++;
			}
			
			// Make sure not out of bounds
			if(j < _list.size())
			{
				accumulated_enthalpy += _list[i].Enthalpy(_list[j]);
			}

			// Increment the loop counter
			i = j;
		}
		return accumulated_enthalpy;
	}

	float DataSequence::Entropy(float salt_concentration) const
	{
		float accumulated_entropy = 0.0;
		for(int i = 0;i < _list.size() - 1;) // Third argument intentionally blank.
		{
			// Skip deletions in first nucleotide
			if(_list[i].GetMostCommon() == '-')
			{
				continue;
			}
			int j = i + 1;
			
			// Skip deletions in second nucleotide
			while(_list[j].GetMostCommon() == '-' && j < _list.size())
			{
				j++;
			}
			
			// Make sure not out of bounds
			if(j < _list.size())
			{
				accumulated_entropy += _list[i].Entropy(_list[j]);
			}

			// Increment the loop counter
			i = j;
		}
		accumulated_entropy += 0.368 * (float)(_list.size()) * log(salt_concentration * 50.0/1000.0);
		// accumulated_entropy += (0.368 * (float)(_list.size()) * log(salt_concentration)); (is this wrong?)
		return accumulated_entropy;
	}

	float DataSequence::Gibbs(float temperature, float salt_concentration) const
	{
		float gibbs = 0.0;
		DataSequence p;
		
		// If a nearest neighbor is a deletion, we need to skip over that particular 
		// DataNode and get a subsequence that does not contain any deletions.
		for(int i = 0; i < _list.size();i++)
		{
			if(_list[i].GetMostCommon() != '-')
			{
				p.PushBack(_list[i]);
			}
		}
		
		// Gibbs Energy = Enthalpy - (Temperature * Salt Corrected Entropy)
		gibbs = p.Enthalpy();
		gibbs -= ((temperature + 273.15) * (p.Entropy(salt_concentration) / 1000.0));
		return gibbs;
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

	int DataSequence::IndexOf(DataSequence data) const
	{
		int index = -1;
		bool outerflag = false;
		for(int i = 0;i < _list.size() - data.size();i++)
		{
			bool innerflag = false;
			for(int j = 0;j < data.size();j++)
			{
				if(innerflag)
				{
					break;
				}
				else if(_list[i+j].GetMostCommon() != data.GetDataSequence()[j].GetMostCommon())
				{
					innerflag = true;
					break;
				}
			}
			if(!innerflag)
			{
				outerflag = true;
				index = i;
			}
		}
		return index;
	}

	float DataSequence::Temperature() const
	{
		float temperature = 0.0;
		if(_list.size() > MAX_PRIMER_LENGTH || _list.size() < MIN_PRIMER_LENGTH)
		{
			temperature = -1.0;
		}
		else
		{
			int gc = 0;
			int at = 0;
			for(int i = 0;i < _list.size();i++)
			{
				switch(_list[i].GetMostCommon())
				{
					case 'A':
					case 'T':
						at++;
						break;
					case 'C':
					case 'G':
						gc++;
						break;
					default:
						break;
				}
			}
			temperature = 64.9 + (41.0 * (gc - 16.4) / (float)(gc + at) ); // Source, northwestern univ.
		}
		return temperature;
	}

	bool DataSequence::checkMatch(DataSequence data) const
	{
		if(data.size() != _list.size())
		{
			return false;
		}
		for(int i = 0;i < data.size();i++)
		{
			if(data.GetDataSequence()[i].GetMostCommon() != _list[i].GetMostCommon())
			{
				return false;
			}
		}
		return true;
	}
} // End of DeGenPrime