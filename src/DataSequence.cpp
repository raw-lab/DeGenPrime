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

	void DataSequence::Print(int id, float temperature, float salt_conc, float primer_conc)
	{
		double accum_ratio = 0.0;
		double accum_weight = 0.0;
		cout << "DataSequence ID [" << id;
		cout << "]\tLength: [" << _list.size() << "]" << endl;
		cout << "Sequence Codes: ";
		PrintCodes();
		cout << "\nBasic Melting Temperature: " << BasicTemperature() << endl;
		cout << "Nearest Neighbor Melting Temperature: " << NNMeltingTemperature(salt_conc, primer_conc) << endl;
		return;
	}

	void DataSequence::PrintCodes()
	{
		for(int i = 0;i < _list.size();i++)
		{
			cout << _list[i].GetCode();
		}
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
		DataSequence p;

		//	Get a subsequence that contains no deletions
		for(int i = 0;i < _list.size();i++)
		{
			if(_list[i].GetMostCommon() != '-')
			{
				p.PushBack(_list[i]);
			}
		}

		// Add up NN enthalpies of the subsequence
		for(int i = 0; i < p.size() - 1;i++)
		{
			accumulated_enthalpy += p.GetDataSequence()[i].Enthalpy(p.GetDataSequence()[i + 1]);
		}

		// Add the initialization enthalpy of the first and last sequences.
		char first = p.GetDataSequence()[0].GetMostCommon();
		char last = p.GetDataSequence()[p.size() - 1].GetMostCommon();

		accumulated_enthalpy += (first == 'C' || first == 'G') ? 0.1 : 2.3;
		accumulated_enthalpy += (last == 'C' || last == 'G') ? 0.1 : 2.3;
		
		return accumulated_enthalpy;
	}

	float DataSequence::Entropy() const
	{
		float accumulated_entropy = 0.0;
		DataSequence p;
		
		//	Get a subsequence that contains no deletions
		for(int i = 0;i < _list.size();i++)
		{
			if(_list[i].GetMostCommon() != '-')
			{
				p.PushBack(_list[i]);
			}
		}

		// Add up NN entropies of the subsequence
		for(int i = 0; i < p.size() - 1;i++)
		{
			accumulated_entropy += p.GetDataSequence()[i].Entropy(p.GetDataSequence()[i + 1]);
		}

		// Add the initialization enthalpy of the first and last sequences.
		char first = p.GetDataSequence()[0].GetMostCommon();
		char last = p.GetDataSequence()[p.size() - 1].GetMostCommon();

		accumulated_entropy += (first == 'C' || first == 'G') ? -2.8 : 4.1;
		accumulated_entropy += (last == 'C' || last == 'G') ? -2.8 : 4.1;

		return accumulated_entropy;
	}

	float DataSequence::Gibbs(float temperature) const
	{
		float accumulated_gibbs = 0.0;
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
		
		for(int i = 0;i < p.size() - 1;i++)
		{
			accumulated_gibbs += p.GetDataSequence()[i].Gibbs(p.GetDataSequence()[i + 1], temperature);
		}

		char first = p.GetDataSequence()[0].GetMostCommon();
		char last = p.GetDataSequence()[p.size() - 1].GetMostCommon();

		float accum = (first == 'C' || first == 'G') ? .98 : 1.03;
		accum = (last == 'C' || last == 'G') ? .98 : 1.03;

		accumulated_gibbs += (first == 'C' || first == 'G') ? 0.98 : 1.03;
		accumulated_gibbs += (last == 'C' || last == 'G') ? 0.98 : 1.03;

		return accumulated_gibbs;
	}

	std::vector<DataNode> DataSequence::GetDataSequence() const { return _list; }
	int DataSequence::size() const { return _list.size(); }

	int DataSequence::RevIndex(int index) const
	{
		return (_list.size() - index - 1);
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

	float DataSequence::BasicTemperature() const
	{
		float temperature = 0.0;
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
		if(_list.size() < 14)
		{
			temperature = (2.0 * at) + (gc * 4.0);
		}
		else
		{
			temperature = 64.9 + (41.0 * (gc - 16.4) / (float)(gc + at) ); // Source, northwestern univ.
		}
		return temperature;
	}

	float DataSequence::RlnK(float primer_conc) const
	{
		float r = 1.9872;
		r *= log(pow(10,9)/primer_conc);
		return r;
	}

	float DataSequence::MonoIonMod(float salt_conc) const
	{
		return (16.6 * (log(salt_conc * pow(10,-3))/log(10.0)));
	}

	float DataSequence::NNMeltingTemperature(float salt_conc, float primer_conc) const
	{
		float numerator = -1.0 * Enthalpy();
		float denominator = -1.0 * Entropy();

		float rlnk = RlnK(primer_conc/4.0);
		float salt_mod = 0.368 * size() * log(salt_conc*pow(10,-3));

		denominator += rlnk + salt_mod;
		denominator /= 1000.0;

		float offset = MonoIonMod(salt_conc);
		float kelvin = (numerator/denominator) + offset;
		return (kelvin - 273.15);
	}

	float DataSequence::BasicAnneal(DataSequence product, float salt_conc, float primer_conc)
	{
		float temperature = 0.3 * NNMeltingTemperature(salt_conc, primer_conc);
		temperature += (0.7 * product.NNMeltingTemperature(salt_conc, primer_conc));
		temperature -= 14.9;
		return temperature;
	}

	int DataSequence::CountMatches(DataSequence data) const
	{
		int matches = 0;
		if(data.size() != _list.size())
		{
			return -1;
		}
		for(int i = 0; i < data.size();i++)
		{
			if(data.GetDataSequence()[i].GetMostCommon() == _list[i].GetMostCommon())
			{
				matches++;
			}
		}
		return matches;
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