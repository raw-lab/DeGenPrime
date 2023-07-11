// DataSequence.cpp
#include <bits/stdc++.h>
#include <cctype>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include "DataNode.h"
#include "DataSequence.h"
#include "GlobalSettings.h"
#include "Primer.h"
#include "global.h"

using namespace std;

namespace DeGenPrime
{
	DataSequence::DataSequence() { }
	DataSequence::DataSequence(std::string str) 
	{
		for(char c : str)
		{
			DataNode node(c,c,1.0);
			_list.push_back(node);
		}
	}

	void DataSequence::SetList(std::vector<DataNode> catalog) {_list = catalog; }
	void DataSequence::PushBack(DataNode node) { _list.push_back(node); }
	void DataSequence::PopBack() { _list.pop_back(); }
	void DataSequence::Erase(int index)
	{
		_list.erase(_list.begin() + index);
	} 

	std::string DataSequence::Print()
	{
		string ret = "Sequence MC: ";
		ret += MC();
		ret += "\nSequence Codes: ";
		ret += Codes();
		ret += "\nLength: " + to_string(_list.size()) + " bps\n";
		ret += "Enthalpy: " + to_string(Enthalpy()) + " kcal\n";
		ret += "Entropy: " + to_string(Entropy()) + " cal\n";
		ret += "Gibbs: " + to_string(Gibbs()) + " kcal\n";
		ret += "GC Content: " + to_string(GCRatio()) + "\n";
		ret += "Basic Melting Temperature: " + to_string(BasicTemperature()) + " deg C\n";
		ret += "Nearest Neighbor Melting Temperature: " + to_string(NNMeltingTemperature()) + " deg C\n";
		ret += "\n";
		return ret;
	}

	std::string DataSequence::Codes() const
	{
		string ret = "";
		for(int i = 0;i < _list.size();i++)
		{
			ret += _list[i].GetCode();
		}
		return ret;
	}

	std::string DataSequence::MC()
	{
		string ret = "";
		for(int i = 0;i < _list.size();i++)
		{
			ret += _list[i].GetMostCommon();
		}
		return ret;
	}

	std::string DataSequence::SectionCodes(int index, int length)
	{
		string ret = "Con (";
		int temp = index;
		int mod = 4;
		while(temp > 0)
		{
			temp /= 10;
			mod--;
		}
		for(int i = 0;i < mod;i++)
		{
			ret += "0";
		}
		ret += to_string(index) + "-";
		temp = index + length;
		mod = 4;
		while(temp > 0)
		{
			temp /= 10;
			mod--;
		}
		for(int i = 0;i < mod;i++)
		{
			ret += "0";
		}
		ret += to_string(index + length - 1) + ") ";
		string codes = Codes();
		codes = codes.substr(index, length);
		ret += codes + "\n";
		return ret;
	}

	std::string DataSequence::SeqInfo()
	{
		string ret = "";
		for(int i = 0;i < _list.size();i++)
		{
			ret += "Index [" + to_string(i);
			ret += "] " + _list[i].NodeInfo();
		}
		return ret;
	}

	string DataSequence::Consensus(std::vector<int> fwd_ind, std::vector<int> fwd_len, bool cons)
	{
		string codes = Codes();
		string ret = "";
		for(int i = 0;i < 76;i++)
		{
			ret += "-";
		}
		int len = ret.length();
		string ret2 = ret;
		ret += "\n\tSize of dashes: " + to_string(len);
		ret += "\n" + ret2 + "\n";
		transform(codes.begin(), codes.end(), codes.begin(), ::tolower);
		if(cons)
		{
			ret += "\t(conserved regions capitalized.)\n"; 
			for(int i = 0;i < fwd_ind.size();i++)
			{
				int len = fwd_len[i];
				string sub = codes.substr(fwd_ind[i], fwd_len[i]);
				transform(sub.begin(), sub.end(), sub.begin(), ::toupper);
				codes.replace(fwd_ind[i], fwd_len[i], sub);
			}
		}
		int final_index = codes.length() - 1;
		int lines = codes.length() / 60;
		int index = 0;
		// Model
		// ' (#####-#####): '
		for(int i = 0;i < lines;i++)
		{
			index = i * 60;
			int mod = 5;
			int temp = index;
			if(temp == 0)mod--;
			while(temp > 0)
			{
				temp /= 10;
				mod--;
			}
			ret += " (";
			for(int j = 0;j < mod;j++)
			{
				ret += "0";
			}
			ret += to_string(index) + "-";
			mod = 5;
			temp = index + 59;
			while(temp > 0)
			{
				temp /= 10;
				mod--;
			}
			for(int j = 0;j < mod;j++)
			{
				ret += "0";
			}
			ret += to_string(index + 59) + "): " + codes.substr(i*60,60) + "\n";
		}
		index += 60;
		codes.erase(0,index);
		int mod = 5;
		int temp = index;
		while(temp > 0)
		{
			temp /= 10;
			mod--;
		}
		ret += " (";
		for(int j = 0;j < mod;j++)
		{
			ret += "0";
		}
		ret += to_string(index) + "-";
		mod = 5;
		temp = final_index;
		while(temp > 0)
		{
			temp /= 10;
			mod--;
		}
		for(int j = 0;j < mod;j++)
		{
			ret += "0";
		}
		ret += to_string(final_index) + "): " + codes + "\n";
		return ret;
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

	float DataSequence::Gibbs() const
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
			accumulated_gibbs += p.GetDataSequence()[i].Gibbs(p.GetDataSequence()[i + 1]);
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
	int DataSequence::ActualSize() const
	{
		int s = 0;
		for(int i = 0;i < _list.size();i++)
		{
			if(_list[i].GetMostCommon() != '-')s++;
		}
		return s;
	}

	int DataSequence::RevIndex(int index) const
	{
		int ret = _list.size();
		ret -= index;
		ret--;
		return (ret);
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
			// Source, https://www.pnas.org/doi/10.1073/pnas.95.4.1460
			temperature = 64.9 + (41.0 * (gc - 16.4) / (float)(gc + at) ); 
		}
		return temperature;
	}

	float DataSequence::RlnK() const
	{
		float r = 1.987;
		r *= log((pow(10,-9) * (GlobalSettings::GetPrimerConcentration())/4.0));
		return r;
	}

	float DataSequence::MonoIonMod() const
	{
		return (16.6 * (log(GlobalSettings::GetMonoIonConcentration() * pow(10,-3))/log(10.0)));
	}

	float DataSequence::NNMeltingTemperature() const
	{
		float numerator = Enthalpy();
		float denominator = Entropy();
		float rlnk = RlnK();
		float salt_mod = 0.368 * (size() - 1) * log(GlobalSettings::GetMonoIonConcentration()*pow(10,-3));
		denominator += salt_mod;
		denominator += rlnk;
		denominator /= 1000.0;
		float kelvin = numerator/denominator;
		return (kelvin - 273.15);
	}

	float DataSequence::GCRatio() const
	{
		int gc = 0;
		for(int i = 0;i < size();i++)
		{
			if(_list[i].GetMostCommon() == 'C' ||
				_list[i].GetMostCommon() == 'G')
			{
				gc++;
			}
		}
		float ratio = (float)gc;
		ratio /= (float)size();
		return ratio;
	}

	float DataSequence::BasicAnneal(DataSequence product)
	{
		
		float temperature = 0.3 * NNMeltingTemperature();
		// cout << "Tm * 0.3 = " << temperature << endl;
		temperature += (0.7 * product.ProductMelt());
		// cout << "Add 0.7 * product melt = " << (0.7 * product.ProductMelt()) << endl; 
		temperature -= 14.9;
		// cout << "Subtract 14.9, final result = " << temperature << endl;
		return temperature;
	}

	float DataSequence::ProductMelt() const
	{
		float temperature = 81.5;
		// float temperature = 0.0;
		temperature += 16.6 * log(GlobalSettings::GetMonoIonConcentration()* pow(10,-3))/log(10);
		temperature += 41.0 * GCRatio();
		// Different sources use different values for the following constant.
		// All values range between 500 and 750.
		float b = 675.0;
		b /= ActualSize();
		temperature -= b;
		return temperature;
	}

	float DataSequence::AverageRatio() const
	{
		float average = 0.0;
		for(DataNode node : _list)
		{
			average += node.Ratio();
		}
		average /= _list.size();
		return average;
	}

	float DataSequence::Penalty() const
	{
		float penalty = 0.0; // The lower this number the higher the quality.

		// Primers outside of size range should not be considered
		if(size() < MIN_PRIMER_LENGTH || size() > MAX_PRIMER_LENGTH)
		{
			return 10000.0;
		}

		// Internal repetition
		// First we want to measure nucleotide repetition in the sequence
		int most_repeats = 0;
		for(char c : {'A', 'C', 'G', 'T'})
		{
			int char_match_count = 0;
			for(int i = 0;i < _list.size() - 1;i++)
			{
				if(_list[i].GetMostCommon() == c &&
					_list[i+1].GetMostCommon() == c)
				{
					char_match_count++;
				}
			}
			if(char_match_count > most_repeats)
			{
				most_repeats = char_match_count;
			}
		}
		switch(most_repeats)
		{
			case 0:
			case 1:
			case 2:
				break;
			case 3:
				penalty += 0.5;
				break;
			case 4:
				penalty += 1.0;
				break;
			case 5:
				penalty += 2.0;
				break;
			default:
				penalty += (0.75 * most_repeats);
				break;
		}

		// GC Clamp
		char c = _list[size() - 1].GetMostCommon();
		if(c != 'G' && c != 'C')
		{
			penalty += 1.0;
		}

		// Dimers
		string last = "";
		for(int i = 0;i < 5;i++)
		{
			int index = _list.size() - 6 + i;
			if(_list[index].GetMostCommon() == '-')
			{
				continue;
			}
			last += string(1, _list[index].GetMostCommon());
		}
		DataSequence ending(last);
		float gibbs = ending.Gibbs();
		penalty += (gibbs < -3.0) ? ((-1.0 * gibbs) - 3.0) : 0;

		// Hairpins
		for(int j = 0;j + 4 < _list.size();j++)
		{
			char first = _list[j].GetMostCommon();
			char second = _list[j+4].GetMostCommon();
			if(first == '-' || second == '-')
			{
				continue;
			}
			switch(second)
			{
				case 'A':
					second = 'T';
					break;
				case 'C':
					second = 'G';
					break;
				case 'G':
					second = 'C';
					break;
				case 'T':
					second = 'A';
					break;
				default:
					continue;
					break;
			}
			if(first != second)
			{
				continue;
			}
			first = _list[j+1].GetMostCommon();
			if(first != 'G')
			{
				continue;
			}
			second = _list[j+3].GetMostCommon();
			if(second == 'A')
			{
				penalty += 0.5;
			}
		}
		for(int j = 0;j + 5 < _list.size();j++)
		{
			char first = _list[j].GetMostCommon();
			char second = _list[j+5].GetMostCommon();
			if(first == '-' || second == '-')
			{
				continue;
			}
			switch(second)
			{
				case 'A':
					second = 'T';
					break;
				case 'C':
					second = 'G';
					break;
				case 'G':
					second = 'C';
					break;
				case 'T':
					second = 'A';
					break;
				default:
					continue;
					break;
			}
			if(first != second)
			{
				continue;
			}
			first = _list[j+1].GetMostCommon();
			if(first != 'G')
			{
				continue;
			}
			second = _list[j+4].GetMostCommon();
			if(second == 'A')
			{
				penalty += 1.0;
			}
		}

		// GC Content
		if(GCRatio() < MIN_GC_TOTAL_RATIO)
		{
			penalty += 25.0 * (MIN_GC_TOTAL_RATIO - GCRatio());
		}
		else if(GCRatio() > MAX_GC_TOTAL_RATIO)
		{
			penalty += 25.0 * (GCRatio() - MAX_GC_TOTAL_RATIO);
		}
		float gc = 0.0;
		for(int i = _list.size() - 6;i < _list.size();i++)
		{
			if(_list[i].GetMostCommon() == 'C' || _list[i].GetMostCommon() == 'G')gc++;
		}
		gc /= 5.0;
		if(gc > MAX_GC_EXTREMA_RATIO)
		{
			penalty += 27.0 * (gc - MAX_GC_EXTREMA_RATIO);
		}

		// Temperature
		if(NNMeltingTemperature() < MIN_PRIMER_TEMP)
		{
			penalty += 10.0 * (MIN_PRIMER_TEMP - NNMeltingTemperature());
		}
		else if(NNMeltingTemperature() > MAX_PRIMER_TEMP)
		{
			penalty += 10.0 * (NNMeltingTemperature() - MAX_PRIMER_TEMP);
		}

		// Degeneracy and Deletions
		string nondegen = "ACGT";
		string degen2 = "WSRYKM";
		string degen3 = "BDHV";
		int del_count = 0;
		for(int i = 0;i < _list.size();i++)
		{
			string s = string(1, _list[i].GetCode());
			if(nondegen.find(s) != string::npos)
			{
				continue;
			}
			else if(s == "-")
			{
				del_count++;
				continue;
			}
			else if(degen2.find(s) != string::npos)
			{
				penalty += 2.0;
				continue;
			}
			else if(degen3.find(s) != string::npos)
			{
				penalty += 6.0;
				continue;
			}
			else
			{
				penalty += 12.0;
			}
		}
		penalty += (8.0 * del_count);
		return penalty;
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

	bool DataSequence::isEmpty() const
	{
		bool ret = true;
		for(DataNode node : _list)
		{
			if(ret == false)
			{
				break;
			}
			if(node.GetMostCommon() != '-')
			{
				ret = false;
			}
		}
		return ret;
	}
} // End of DeGenPrime