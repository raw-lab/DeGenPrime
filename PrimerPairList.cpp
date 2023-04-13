// PrimerPairList.cpp

#include <algorithm>
#include <iostream>
#include <vector>
#include "PrimerPairList.h"
#include "PrimerPair.h"
#include "Primer.h"
#include "DataSequence.h"
#include "global.h"

using namespace std;

namespace DeGenPrime
{
	PrimerPairList::PrimerPairList() { }
	PrimerPairList::PrimerPairList(DataSequence fwd_seq, DataSequence rev_seq, std::vector<PrimerPair> pair_list)
	{
		_fwd = fwd_seq;
		_rev = rev_seq;
		_pairs = pair_list;
	}

	PrimerPairList PrimerPairList::SubList(int startIndex, int length)
	{
		PrimerPairList sub_list;
		for(int i = startIndex;i < startIndex+length;i++)
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

	void PrimerPairList::CreateList(DataSequence fwd_seq, DataSequence rev_seq, std::vector<Primer> fwd_list, std::vector<Primer> rev_list)
	{
		_fwd = fwd_seq;
		_rev = rev_seq;
		for(int i = 0;i < fwd_list.size();i++)
		{
			for(int j = 0;j < rev_list.size();j++)
			{
				PrimerPair p(fwd_list[i], rev_list[j], fwd_seq, rev_seq);
				_pairs.push_back(p);
			}
		}
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

	void PrimerPairList::FilterAmpliconLength()
	{
		for(int i = _pairs.size() - 1; i >= 0;i--)
		{
			if(MIN_AMPLICON_LENGTH > _pairs[i].AmpSize())
			{
				Erase(i);
			}
		}
	}

	void PrimerPairList::FilterTemperatureDifference()
	{
		// We want to filter any primers whose difference in melting temperature
		// is greater than MAX_TEMP_DIFFERENCE.
		for(int i = _pairs.size() - 1;i >= 0;i--)
		{
			if(MAX_TEMP_DIFFERENCE < _pairs[i].TempDiff())
			{
				Erase(i);
			}
		}
	}

	int PrimerPairList::FilterAnnealingTemp(DataSequence fwd, DataSequence rev)
	{
		// We want to filter any primer pairs if the annealing temperature of
		// the least stable (highest Gibbs @ melting temp) primer is 10 deg C
		// more or less than the NNMeltingTemperature at std conditions of that primer.
		// the return value is the number of primers filtered.
		int filtered = 0;
		for(int i = _pairs.size() - 1;i >= 0;i--)
		{
			// cout << "Inside FilterAnnealingTemp Loop, Iteration[" << i << "]" << endl;
			DataSequence primer_fwd = fwd.SubSeq(_pairs[i].GetForward().Index(),_pairs[i].GetForward().Length());
			DataSequence primer_rev = rev.SubSeq(_pairs[i].GetReverse().Index(),_pairs[i].GetReverse().Length());

			// primer_fwd.Print(i,37.0,50.0,50.0);
			// primer_rev.Print(i,37.0,50.0,50.0);

			float fwd_gibbs = primer_fwd.Gibbs(50.0);
			float rev_gibbs = primer_rev.Gibbs(50.0);

			DataSequence most_stable = (fwd_gibbs > rev_gibbs) ? primer_fwd : primer_rev;
			float temp = (fwd_gibbs > rev_gibbs) ? primer_fwd.NNMeltingTemperature(50.0,50.0) 
				: primer_rev.NNMeltingTemperature(50.0,50.0);

			// cout << "Just after finding the temp, which was: " << temp << endl;
			
			DataSequence product = (fwd_gibbs > rev_gibbs) ? fwd.SubSeq(_pairs[i].GetForward().Index(), _pairs[i].AmpSize())
				: rev.SubSeq(_pairs[i].GetReverse().Index(), _pairs[i].AmpSize());

			// cout << "Just before calling BasicAnneal" << endl;

			float anneal = most_stable.BasicAnneal(product, 50.0,50.0);
			float diff = (temp - anneal);
			if(diff < -10.0 || diff > 10.0)
			{
				Erase(i);
				filtered++;
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
	
	void PrimerPairList::PrintSize()
	{
		cout << "The number of forward-reverse primer pairs in this list is: ";
		cout << size() << endl;
	}

	void PrimerPairList::PrintAll(DataSequence fwd, DataSequence rev)
	{
		for(int i = 0;i < _pairs.size();i++)
		{
			cout << "Primer pair #" << i + 1 << endl;
			_pairs[i].Print(fwd, rev);
		}
	}

	std::vector<PrimerPair> PrimerPairList::GetPairs() const { return _pairs; }

	int PrimerPairList::size() const
	{
		return (_pairs.size());
	}
} // End of DeGenPrime

