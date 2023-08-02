// PrimerPairList.cpp

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "primerpairlist.h"
#include "primerpair.h"
#include "primer.h"
#include "datasequence.h"
#include "globalsettings.h"
#include "format.h"
#include "global.h"

using namespace std;

namespace DeGenPrime
{
	PrimerPairList::PrimerPairList() 
	{
		_OriginalSize = 0;
	}
	PrimerPairList::PrimerPairList(DataSequence fwd_seq, DataSequence rev_seq, std::vector<PrimerPair> pair_list)
	{
		_pairs = pair_list;
		_OriginalSize = size();
	}

	PrimerPairList PrimerPairList::SubList(int startIndex, int length)
	{
		PrimerPairList sub_list;
		int amt = (_pairs.size() <= length) ? _pairs.size() : length; 
		for(int i = startIndex;i < startIndex+amt;i++)
		{
				sub_list.PushBack(_pairs[i]);
		}
		return sub_list;

	}

	void PrimerPairList::Append(PrimerPairList list)
	{
		for(int i = 0;i < list.size();i++)
		{
			_pairs.push_back(list.GetPairs()[i]);
		}
	}

	void PrimerPairList::CreateFromRange(DataSequence fwd_seq, DataSequence rev_seq,
		std::vector<Primer> fwd_list, std::vector<Primer> rev_list, 
		int fwd_begin, int fwd_end, int rev_begin, int rev_end)
	{
		int x_begin = fwd_begin >= 0 ? fwd_begin : 0;
		int x_end = fwd_end > fwd_seq.size() ? fwd_seq.size() : fwd_end;
		int y_begin = rev_begin >= 0 ? rev_begin : 0;
		int y_end = rev_end > rev_seq.size() ? rev_seq.size() : rev_end;
		for(int i = x_begin;i < x_end;i++)
		{
			for(int j = y_begin;j < y_end;j++)
			{
				PrimerPair p(fwd_list[i], rev_list[j], fwd_seq, rev_seq);
				_pairs.push_back(p);
			}
		}
		_OriginalSize = _pairs.size();
	}

	void PrimerPairList::CreateList(DataSequence fwd_seq, DataSequence rev_seq, std::vector<Primer> fwd_list, std::vector<Primer> rev_list)
	{
		for(int i = 0;i < fwd_list.size();i++)
		{
			for(int j = 0;j < rev_list.size();j++)
			{
				PrimerPair p(fwd_list[i], rev_list[j], fwd_seq, rev_seq);
				_pairs.push_back(p);
			}
		}
		_OriginalSize = _pairs.size();
	}

	void PrimerPairList::Erase(int index)
	{
		_pairs.erase(_pairs.begin() + index);
	}

	void PrimerPairList::PushBack(PrimerPair pair)
	{
		_pairs.push_back(pair);
	}

	void PrimerPairList::Sort()
	{
		sort(_pairs.begin(), _pairs.end());
	}

	string PrimerPairList::FilterAmpliconLength()
	{
		string ret;
		int filtercount = 0;
		for(int i = _pairs.size() - 1; i >= 0;i--)
		{
			if(GlobalSettings::GetMinimumAmplicon() > _pairs[i].AmpSize())
			{
				Erase(i);
				filtercount++;
			}
		}
		ret = FilterMessage("FilterAmpLength", filtercount);
		return ret;
	}

	string PrimerPairList::FilterTemperatureDifference()
	{
		string ret;
		int filtercount = 0;
		// We want to filter any primers whose difference in melting temperature
		// is greater than MAX_TEMP_DIFFERENCE.
		for(int i = _pairs.size() - 1;i >= 0;i--)
		{
			if(MAX_TEMP_DIFFERENCE < _pairs[i].TempDiff())
			{
				Erase(i);
				filtercount++;
			}
		}
		ret = FilterMessage("FilterTempDiff", filtercount);
		return ret;
	}

	string PrimerPairList::FilterMessage(std::string func, int filtercount)
	{
		bool start = func == "start";
		bool final = func == "final";
		string ret = "";
		string line = "";
		if(start)
		{
			line = "Filter Name";
			ret += Format(line, 25, Alignment::Center);
			line = "Filtered  ";
			ret += Format(line, 25, Alignment::Right);
			line = "Percent";
			ret += Format(line, STR_FORMAT - 50, Alignment::Left);
		}
		else if(final)
		{
			line = "";
			for(int i = 0;i < STR_FORMAT;i++)
			{
				line += "-";
			}
			ret += line + "\n";
			line = "Remaining";
			ret += Format(line, 25, Alignment::Left);
			line = to_string(size()) + "   ";
			ret += Format(line, 25, Alignment::Right);
			line = "(";
			float percent = 100.0 * ((float)size())/((float)_OriginalSize);
			line += Format(percent, 2) + "%)";
			ret += Format(line, STR_FORMAT - 50, Alignment::Left) + "\n";
		}
		else
		{
			line = "     " + func;
			ret += Format(line, 25, Alignment::Left);
			line = to_string(filtercount) + "   ";
			ret += Format(line, 25, Alignment::Right);
			line = "(";
			float percent = 100.0 * ((float)filtercount)/((float)_OriginalSize);
			line += Format(percent, 2) + "%)";
			ret += Format(line, STR_FORMAT - 50, Alignment::Left);
		}
		return ret;
	}

	int PrimerPairList::FilterAnnealingTemp(DataSequence fwd, DataSequence rev, int ignore)
	{
		// We want to filter any primer pairs if the annealing temperature of
		// the least stable (highest Gibbs @ melting temp) primer is 10 deg C
		// more or less than the NNMeltingTemperature at std conditions of that primer.
		// the return value is the number of primers filtered.
		int filtered = 0;
		for(int i = _pairs.size() - 1;i >= ignore;i--)
		{
			DataSequence primer_fwd = fwd.SubSeq(_pairs[i].GetForward().Index(),_pairs[i].GetForward().Length());
			DataSequence primer_rev = rev.SubSeq(_pairs[i].GetReverse().Index(),_pairs[i].GetReverse().Length());

			float fwd_gibbs = primer_fwd.Gibbs();
			float rev_gibbs = primer_rev.Gibbs();


			// The most stable is the primer with the lowest gibbs
			DataSequence most_stable = (fwd_gibbs < rev_gibbs) ? primer_fwd : primer_rev;
			DataSequence least_stable = (fwd_gibbs < rev_gibbs) ? primer_rev : primer_fwd;
			float temp = (fwd_gibbs < rev_gibbs) ? primer_fwd.NNMeltingTemperature() 
				: primer_rev.NNMeltingTemperature();
			float l_temp = (fwd_gibbs < rev_gibbs) ? primer_rev.NNMeltingTemperature()
				: primer_fwd.NNMeltingTemperature();
			
			DataSequence product = (fwd_gibbs < rev_gibbs) ? fwd.SubSeq(_pairs[i].GetForward().Index(), _pairs[i].AmpSize())
				: rev.SubSeq(_pairs[i].GetReverse().Index(), _pairs[i].AmpSize());
			DataSequence l_product = (fwd_gibbs < rev_gibbs) ? rev.SubSeq(_pairs[i].GetReverse().Index(), _pairs[i].AmpSize())
				: fwd.SubSeq(_pairs[i].GetForward().Index(), _pairs[i].AmpSize());

			float anneal = most_stable.BasicAnneal(product);
			float l_anneal = least_stable.BasicAnneal(l_product);
			float diff = (temp - anneal);
			float l_diff = (l_temp - l_anneal);
			
			if((diff < -10.0 || diff > 10.0) && (l_diff < -10 || l_diff > 10.0))
			{
				Erase(i);
				filtered++;
			}
		}
		return filtered;
	}

	int PrimerPairList::FilterUnique()
	{
		int filtered = 0;
		for(int i = _pairs.size() - 1;i >= 1;i--)
		{
			int j = 0;
			int fwd_index = _pairs[i].GetForward().Index();
			int rev_index = _pairs[i].GetReverse().Index();
			while(i > j)
			{
				int f_index = _pairs[j].GetForward().Index();
				int r_index = _pairs[j].GetReverse().Index();
				if(fwd_index == _pairs[j].GetForward().Index() ||
					rev_index == _pairs[j].GetReverse().Index())
				{
					Erase(i);
					filtered++;
					break;
				}
				else
				{
					j++;
				}
			}
		}
		return filtered;
	}

	bool PrimerPairList::comparator(const PrimerPair& lhs, const PrimerPair& rhs)
	{
		float temp1 = lhs.TempDiff();
		float temp2 = rhs.TempDiff();
		return (temp1 < temp2);
	}

	string PrimerPairList::PrintAll(DataSequence fwd, DataSequence rev)
	{
		string ret = "";
		string line = "";
		for(int i = 0;i < _pairs.size();i++)
		{
			line = "Primer Pair #" + to_string(i + 1);
			ret += Format(line, 25, Alignment::Center);
			ret += _pairs[i].Print(fwd, rev) + "\n";
		}
		return ret;
	}

	string PrimerPairList::PrintAllShort(DataSequence fwd, DataSequence rev)
	{
		string ret = "";
		for(int i =0;i < _pairs.size();i++)
		{
			ret += "Pair #" + to_string(i + 1) + " ";
			ret += _pairs[i].PrintShort(fwd, rev);
		}
		return ret;
	}

	std::vector<PrimerPair> PrimerPairList::GetPairs() const { return _pairs; }

	int PrimerPairList::PartitionCount(int fwd_size, int rev_size) const
	{
		int rect = fwd_size * rev_size;
		const int area = 1600;
		double ratio = (float)rect/(float)area;
		if(ratio <= 1.0)
		{
			return 1;
		}
		else
		{
			int r = (int)ceil(ratio);
			int n = 2;
			while(pow(n,2) < r)
			{
				n++;
			}
			return (pow(n,2));
		}
	}

	int PrimerPairList::size() const
	{
		return (_pairs.size());
	}
} // End of DeGenPrime

