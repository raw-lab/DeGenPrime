// DataSequence.cpp
#include <algorithm>
#include <cctype>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include "datanode.h"
#include "datasequence.h"
#include "format.h"
#include "globalsettings.h"
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
		string ret = Format("Basic Information", STR_FORMAT, Alignment::Center) + "\n";
		string line = "Sequence Codes: " + Codes();
		ret += Format(line, 45, Alignment::Right);
		line = "       Length: " + to_string(_list.size()) + " bps";
		ret += Format(line, STR_FORMAT - 40, Alignment::Left) + "\n";
		line = "Sequence MC: " + MC();
		ret += Format(line, 45, Alignment::Right);
		line = "   GC Content: ";
		line += Format((float)100.0 * GCRatio(), 2) + "%";
		ret += Format(line, STR_FORMAT - 45, Alignment::Left) + "\n";
		ret += Format("Thermodynamic Information", STR_FORMAT, Alignment::Center) + "\n";
		line = "Enthalpy: ";
		line += Format(Enthalpy(), 2) + " kcal";
		ret += Format(line, 25, Alignment::Center);
		line = "Entropy: ";
		line += Format(Entropy(), 2) + " cal";
		ret += Format(line, 25, Alignment::Center);
		line = "Gibbs: ";
		line += Format(Gibbs(), 2) + " kcal";
		ret += Format(line, STR_FORMAT - 50, Alignment::Center) + "\n";
		ret += Format("Melting Temperature", STR_FORMAT, Alignment::Center) + "\n";
		line = "Basic: ";
		line += Format(BasicTemperature(), 2) + " deg C";
		ret += Format(line, 30, Alignment::Center);
		line = "Nearest Neighbor: ";
		line += Format(NNMeltingTemperature(), 2) + " deg C";
		ret += Format(line, STR_FORMAT - 30, Alignment::Center) + "\n";
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

	string DataSequence::Consensus(std::vector<int> fwd_ind, std::vector<int> fwd_len, bool cons)
	{
		string ret = "";
		string codes = Codes();
		string start_num, end_num, range;
		transform(codes.begin(), codes.end(), codes.begin(), ::tolower);
		if(cons)
		{
			ret += Format("Conserved regions capitalized.", STR_FORMAT, Alignment::Center) + "\n";
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
		for(int i = 0;i < lines;i++)
		{
			index = i * 60;
			start_num = Format(index,5);
			end_num = Format(index + 59, 5);
			range = "(" + start_num + "-" + end_num + "): ";
			ret += Format(range, 16, Alignment::Right);
			ret += codes.substr(i*60,60) + "\n";
		}
		index += 60;
		codes.erase(0,index);
		start_num = Format(index, 5);
		end_num = Format(final_index, 5);
		range = "(" + start_num + "-" + end_num + "): ";
		ret += Format(range, 16, Alignment::Right) + codes + "\n";
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

		if(p.size() < 2)return 0.0; // Avoids Seg faults for mini-seqs.

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

		if(p.size() < 2)return 0.0; // Avoids seg faults for mini-seqs.

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

		if(p.size() < 2)return 0.0; // Avoids seg faults for mini-seqs.
		
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
		temperature += (0.7 * product.ProductMelt());
		temperature -= 14.9;
		return temperature;
	}

	float DataSequence::ProductMelt() const
	{
		float temperature = 81.5;
		temperature += 16.6 * log(GlobalSettings::GetMonoIonConcentration()* pow(10,-3))/log(10);
		temperature += 41.0 * GCRatio();
		// Different sources use different values for the following constant.
		// All values range between 500 and 750.
		float b = 675.0;
		b /= ActualSize();
		temperature -= b;
		return temperature;
	}

	float DataSequence::Penalty() const
	{
		float penalty = 0.0; // The lower this number the higher the quality.
		// Primers outside of size range should not be considered
		int penalty_mod = 1;
		if(size() < GlobalSettings::GetMinimumPrimerLength() ||
			size() > GlobalSettings::GetMaximumPrimerLength())
		{
			return 1000.0;
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
		float temp = NNMeltingTemperature();
		if(temp < GlobalSettings::GetMinimumTemperature())
		{
			if(GlobalSettings::GetUserTemp())
			{
				penalty += 1000;
			}
			else
			{
				penalty += 10.0 * (GlobalSettings::GetMinimumTemperature() - temp);
			}
		}
		else if(temp > GlobalSettings::GetMaximumTemperature())
		{
			if(GlobalSettings::GetUserTemp())
			{
				penalty += 1000;
			}
			else
			{
				penalty += 10.0 * (temp - GlobalSettings::GetMaximumTemperature());
			}
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