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

	void PrimerPairList::FilterCrossDimers()
	{
		for(int i = _pairs.size() - 1;i >= 0;i--)
		{
			bool flag = false;
			int match_count = 0;

			// -	Filter any primers who sequences have more than 3 parallel matches
			DataSequence fwd_primer = _fwd.SubSeq(_pairs[i].GetForward().Index(), _pairs[i].GetForward().Length());
			DataSequence rev_primer = _rev.SubSeq(_pairs[i].GetReverse().Index(), _pairs[i].GetReverse().Length());

			// Invert the reverse primer
			rev_primer = rev_primer.InvSeq();
			
			// consider that primer pairs may contain primers of different sizes.
			int size_difference = fwd_primer.size() - rev_primer.size();
			
			// If size_difference > 0, the forward primer is bigger, need to slide along it
			// If size_difference < 0, the reverse primer is bigger, need to slide along it
			// If size_difference = 0, they are the same and direct comparison is fine.

			// Begin loop at index 3, run through smaller size
			int smallersize = (size_difference >= 0) ? rev_primer.size() : fwd_primer.size();
			for(int j = 3;j < smallersize;j++)
			{
				if(flag)
				{
					break;
				}
				
				// For this algorithm we want to slide the larger primer along the smaller primer.
				// Use the ternary operator on size_difference to find this primer.
				// It does not matter which we use if size_difference = 0.
				
				DataSequence smaller_sub = (size_difference >= 0) ? rev_primer.SubSeq(0,j) : fwd_primer.SubSeq(0,j);
				
				// We need to do some extra slides for size differences when the primers entirely overlap
				if(j == smallersize - 1 && size_difference != 0)
				{
					// Use size_difference as a while loop condition to get all combinations.
					int n = 0;
					int increment = (size_difference >= 0) ? -1 : 1;
					size_difference -= increment; // We need one extra loop
					while(size_difference != 0)
					{
						if(flag)
						{
							break;
						}
						DataSequence larger_sub = (size_difference >= 0) ? fwd_primer.SubSeq(fwd_primer.size() - j-n,j) 
							: rev_primer.SubSeq(rev_primer.size() - j-n,j);
						match_count = larger_sub.CountMatches(smaller_sub);
						flag = (match_count > 3) ? true : false;
						n++;
						size_difference += increment;
					}
				}
				else
				{
					DataSequence larger_sub = (size_difference >= 0) ? fwd_primer.SubSeq(fwd_primer.size() - j,j) 
						: rev_primer.SubSeq(rev_primer.size() - j,j);
					match_count = larger_sub.CountMatches(smaller_sub);
					flag = (match_count > 3) ? true : false;
				}
			}
			
			// This time we run the primers the other way
			size_difference = fwd_primer.size() - rev_primer.size(); // Reset this value
			for(int j = 3;j < smallersize;j++)
			{
				if(flag)
				{
					break;
				}
				
				// For this algorithm we want to slide the larger primer along the smaller primer.
				// Use the ternary operator on size_difference to find this primer.
				// It does not matter which we use if size_difference = 0.
				
				DataSequence smaller_sub = (size_difference >= 0) ? rev_primer.SubSeq(rev_primer.size() - j,j)
					: fwd_primer.SubSeq(fwd_primer.size() - j,j);
				
				// We need to do some extra slides for size differences when the primers entirely overlap
				if(j == smallersize - 1 && size_difference != 0)
				{
					// Use size_difference as a while loop condition to get all combinations.
					int n = 0;
					int increment = (size_difference >= 0) ? -1 : 1;
					size_difference -= increment; // We need one extra loop
					while(size_difference != 0)
					{
						if(flag)
						{
							break;
						}
						DataSequence larger_sub = (size_difference >= 0) ? fwd_primer.SubSeq(n,j) : rev_primer.SubSeq(n,j);
						match_count = larger_sub.CountMatches(smaller_sub);
						flag = (match_count > 3) ? true : false;
						n++;
						size_difference += increment;
					}
				}
				else
				{
					DataSequence larger_sub = (size_difference >= 0) ? fwd_primer.SubSeq(0,j) : rev_primer.SubSeq(0,j);
					match_count = larger_sub.CountMatches(smaller_sub);
					flag = (match_count > 3) ? true : false;
				}
			}
			if(flag)
			{
				Erase(i);
			}
		} // End of Looping through PrimerPairList
	}

	/*
	void PrimerPairList::Sort()
	{
		// Sort using lambda expression
		sort(_pairs.begin(), _pairs.end());
		
	}

	bool PrimerPairList::comparator(const PrimerPair& lhs, const PrimerPair& rhs)
	{
		float temp1 = lhs.TemperatureDifference(_fwd,_rev);
		float temp2 = rhs.TemperatureDifference(_fwd,_rev);
		return (temp1 < temp2);
	}
	*/
	
	void PrimerPairList::PrintSize()
	{
		cout << "The number of forward-reverse primer pairs in this list is: ";
		cout << size() << endl;
	}

	std::vector<PrimerPair> PrimerPairList::GetPairs() const { return _pairs; }

	int PrimerPairList::size() const
	{
		return (_pairs.size());
	}
} // End of DeGenPrime

