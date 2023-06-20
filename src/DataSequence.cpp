// DataSequence.cpp
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include "DataNode.h"
#include "DataSequence.h"
#include "GlobalSettings.h"
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

	std::string DataSequence::Codes()
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

	void DataSequence::Consensus(std::string filename, std::ofstream& ofs)
	{
		string ret = filename + "\n\n";
		string codes = Codes();
		int final_index = codes.length() - 1;
		int lines = codes.length() / 60;
		int remainder = codes.length() % 60;
		int index = 0;
		// Model
		// 'Con (####-####) '
		string start = "";
		for(int i = 0;i < lines;i++)
		{
			index = i * 60;
			ofs << "Con (" << std::setfill('0') << std::setw(4) << index;
			ofs << "-" << std::setfill('0') << std::setw(4) << index + 59;
			ofs << ") " << codes.substr(index,60) << endl;
		}
		index += 60;
		codes.erase(0,index);
		ofs << "Con (" << std::setfill('0') << std::setw(4) << index;
		ofs << "-" << std::setfill('0') << std::setw(4) << final_index;
		ofs << ") " << codes << endl;
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
		// temperature += (0.7 * product.NNMeltingTemperature()); -- Old, slow, incorrect
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
		float b = 600.0;
		b /= size();
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

	float DataSequence::Quality() const
	{
		float quality = 0.0;
		// The lower this number the higher the quality.
		
		// Primers outside of size range should not be considered
		if(size() < MIN_PRIMER_LENGTH || size() > MAX_PRIMER_LENGTH)
		{
			return 0.0;
		}

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
			case 3:
				quality += 1.0;
				break;
			case 4:
				quality += 2.0;
				break;
			case 5:
				quality += 3.0;
				break;
			default:
				break;
		}

		// Check if there is a GC clamp
		int gc = 0;
		for(int i = 1;i <=5;i++)
		{
			if(_list[size() - i].GetMostCommon() == 'G' ||
				_list[size() - i].GetMostCommon() == 'C')
			{
				gc++;
			}
		}
		switch(gc)
		{
			case 0:
				quality += 5.0;
				break;
			case 1:
				quality += 1.0;
				break;
			default:
				break;
		}
		return quality;
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