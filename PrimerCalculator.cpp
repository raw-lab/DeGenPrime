#include <cctype>
#include <cmath>
#include <iostream>
#include <vector>
#include "Primer.h"
#include "SequenceList.h"
#include "DataSequence.h"
#include "PrimerCalculator.h"
#include "global.h"

using namespace std;

namespace DeGenPrime
{
	PrimerCalculator::PrimerCalculator() { }
	
	void PrimerCalculator::InitializePrimers(DataSequence data, int AmpliconLength)
	{
		int endIndex = data.size() - AmpliconLength;

		// Check AmpliconLength is valid
		if(AmpliconLength < 0)
		{
			endIndex = data.size();
		}
		if(endIndex < 0)
		{
			endIndex = 0;
		}

		for(int i = MIN_PRIMER_LENGTH;i <= MAX_PRIMER_LENGTH;i++) // length
		{
			for(int j = 0;j <= endIndex;j++)
			{
				Primer id(j,i);
				PushBack(id);
			}
		}
		/* Old Version
		// Needs to loop through all primers in the forward direction
		for(int i = MIN_PRIMER_LENGTH;i <= MAX_PRIMER_LENGTH;i++) // length
		{
			for(int j = 0;j < (data.size() - i);j++) // index
			{
				Primer id(j,i);
				PushBack(id);
			}
		}*/
	}

	void PrimerCalculator::FilterDegeneracy(DataSequence data)
	{
		for(int i = _primers.size() - 1; i >= 0;i--)
		{
			// Define the primer
			DataSequence p = data.SubSeq(_primers[i].Index(), _primers[i].Length());

			// Filter sequences with too many degenerate codes
			//	- raise flag if N is found anywhere in the sequence, 
			//	- or any degenerate code within 3 of an end
			//	- or more than 1 B,D,H,V
			//	- or more than 2 others. (count B,D,H,V as 2 any other as 1, ends as 3)
			int degeneracy_count = 0;
			std::vector<char> junk;
			for(int j = 0;j < p.size();j++)
			{
				if(degeneracy_count >= 3)
				{
					break;
				}
				int position_multiplier =  ((j < 3) || (j >= p.size() - 3)) ? 3 : 1;
				char c = toupper(p.GetDataSequence()[j].GetCode());
				switch(c)
				{
					case 'N':
						degeneracy_count += 3;
						junk.push_back(c);
						break;
					case 'B':
					case 'D':
					case 'V':
					case 'H':
						degeneracy_count += (2 * position_multiplier);
						junk.push_back(c);
						break;
					case 'A':
					case 'C':
					case 'G':
					case 'T':
					case '-':
						break;
					default:
						degeneracy_count += position_multiplier;
						junk.push_back(c);
						break;
				}
			}
			if(degeneracy_count >= 3)
			{
				/* Testing
				cout << "These chars increased degeneracy: ";
				for(char degen : junk)
				{
					cout << degen << " ";
				}
				cout << endl;*/
				Erase(i);
			}
		}
	}

	void PrimerCalculator::FilterDeletions(DataSequence data, SequenceList list)
	{
		// data can be reversed in argument, list cannot.
		for(int i = _primers.size() - 1; i >= 0;i--)
		{
			DataSequence p = data.SubSeq(_primers[i].Index(), _primers[i].Length());
			bool flag = false;
			int total_deletions_count = 0;
			// Filter sequences for deletions
			//	- Filter primer with any deletion in any sequence on the first/last 3.
			//	- Filter primer with any three consecutive datasequence chars of '-'
			//	- Filter primer with more than 6 total deletions.
			//	- Filter any primer whose deletions make it smaller than MIN_PRIMER_SIZE
			for(int j = 0;j < p.size();j++)
			{
				if(flag)
				{
					break;
				}
				int consecutive_deletions_count = 0;

				if(j < 3 || j > p.size() - 3)
				{
					for(Sequence seq : list.GetSequenceList())
					{
						if(seq.GetCodes()[_primers[i].Index() + j] == '-')
						{
							flag = true;
							break;
						}
					}
				}
				if(p.GetDataSequence()[j].GetCode() == '-')
				{
					total_deletions_count++;
					consecutive_deletions_count++;	
				}
				if(total_deletions_count > 6 || consecutive_deletions_count > 2)
				{
					flag = true;
				}
				if((p.size() - total_deletions_count) < MIN_PRIMER_LENGTH)
				{
					flag = true;
				}
			}
			if(flag)
			{
				Erase(i);
			}
		}
	}

	void PrimerCalculator::FilterGCContent(DataSequence data)
	{
		for(int i = _primers.size() - 1; i >= 0;i--)
		{

			DataSequence p = data.SubSeq(_primers[i].Index(), _primers[i].Length());
			int gc_end_count = 0;
			int gc_total_count = 0;

			// Filter for GC content.
			//	- filter any primer without 40% minimum or 60% maximum gc content
			//	- filter any primer with more than 3 G or C in last five nucleotides
			//	- count by most common to avoid issues with degenerate codes
			for(int j = 0; j < p.size();j++)
			{
				char c = p.GetDataSequence()[j].GetMostCommon();
				if(c == 'C' || c == 'G')
				{
					gc_total_count++;
					if(j >= p.size() - 5)
					{
						gc_end_count++;
					}
				}
			}
			float end_ratio = (float)(gc_end_count) / 5.0;
			float total_ratio = (float)(gc_total_count) / (float)(p.size());
			if(end_ratio > MAX_GC_EXTREMA_RATIO || total_ratio < MIN_GC_TOTAL_RATIO || total_ratio > MAX_GC_TOTAL_RATIO)
			{
				Erase(i);
			}
		}
	}

	void PrimerCalculator::FilterGibbs(DataSequence data, float temperature, float salt_concentration)
	{
		for(int i = _primers.size() - 1;i >= 0;i--)
		{
			// Define the primer
			DataSequence p = data.SubSeq(_primers[i].Index(), _primers[i].Length());

			// Filter sequences that are likely to form dimers or hairpins by Gibbs energy
			// -	filter any primer with internal gibbs < -6
			// -	filter any primer with a 3' end < -5
			DataSequence InternalSeq = p.SubSeq(0, p.size() - 3);
			DataSequence EndSeq = p.SubSeq(p.size() - 3, 3);

			float InternalGibbs = InternalSeq.Gibbs(temperature, salt_concentration);
			float EndGibbs = EndSeq.Gibbs(temperature, salt_concentration);

			if(InternalGibbs < -6.0 || EndGibbs < -5.0)
			{
				Erase(i);
			}
		}
	}

	void PrimerCalculator::FilterRepeats(DataSequence data)
	{
		for(int i = _primers.size() - 1;i >= 0;i--)
		{	
			// Define the primer
			DataSequence p = data.SubSeq(_primers[i].Index(), _primers[i].Length());
			bool flag = false;

			// Filter sequences with too much repition of nucleotides
			//	- raise flag if four matches are found to any pair, more than one codon, or any higher
			//	- break loop if flag or length * (match count) < remaining checks * match count
			//	- there is no need to check subsequences longer than sqrt the total length.
			for(int length = 2; length < sqrt(p.size());length++)
			{
				if(flag)
				{
					break;
				}

				for(int j = 0;j < p.size() - length;j++)
				{
					if(flag)
					{
						break;
					}
					DataSequence first = p.SubSeq(j,length);
					int match_count = 0;
					for(int k = j+1; k < p.size() - length + 1;k++)
					{
						if(flag)
						{
							break;
						}
						DataSequence second = p.SubSeq(k,length);
						if(first.checkMatch(second))
						{
							match_count++;
							flag = (match_count > TooManyRepeats(length));
							if(flag)
							{
								/* Debug
								cout << "Primer [" << _primers[i].Index() << "] failed FilterRepeats.";
								cout << "[match_count=" << match_count << "][TooManyRepeats(l)=";
								cout << TooManyRepeats(length) << "][length=" << length << "]" << endl;
								*/
							}
						}
					}
				}
			}
			if(flag)
			{
				Erase(i);
			}
		}
	}

	void PrimerCalculator::FilterComplementaryEnds(DataSequence data)
	{
		for(int i = _primers.size() - 1;i >= 0;i--)
		{
			DataSequence p = data.SubSeq(_primers[i].Index(), _primers[i].Length());
			DataSequence firstThree = p.SubSeq(0,3);
			DataSequence lastThree = p.SubSeq(p.size() - 3,3);
			lastThree = lastThree.InvSeq().RevSeq();

			if(firstThree.checkMatch(lastThree))
			{
				Erase(i);
			}
		}
	}

	void PrimerCalculator::FilterTemperature(DataSequence data)
	{
		for(int i = _primers.size() - 1; i >= 0;i--)
		{
			DataSequence p = data.SubSeq(_primers[i].Index(), _primers[i].Length());
			if(p.Temperature() < MIN_PRIMER_TEMP || p.Temperature() > MAX_PRIMER_TEMP)
			{
				Erase(i);
			}
		}
	}
	
	void PrimerCalculator::Erase(int index)
	{
		_primers.erase(_primers.begin() + index);
	} 
	void PrimerCalculator::PushBack(Primer primer)
	{
		_primers.push_back(primer);
	}

	int PrimerCalculator::TooManyRepeats(int size)
	{
		int b = 0;
		switch(size)
		{
			case 2:
				b = 3;
				break;
			case 3:
				b = 2;
				break;
			default:
				b = 1;
				break;
		}
		return b;
	}

	void PrimerCalculator::PrintSize()
	{
		cout << "The number of primers in the list is: " << size() << endl;
	}

	void PrimerCalculator::PrintAll()
	{
		for(Primer p : _primers)
		{
			p.Print();
		}
		cout << endl;
	}

	void PrimerCalculator::SetPrimers(std::vector<Primer> primerList) { _primers = primerList; }

	std::vector<Primer> PrimerCalculator::GetPrimers() const { return _primers; }

	int PrimerCalculator::size() const { return _primers.size(); }

	int PrimerCalculator::AmpliconLength(DataSequence fwd_data, DataSequence rev_data, Primer forward, Primer reverse) const
	{
		DataSequence rev = rev_data.SubSeq(reverse.Index(), reverse.Length());
		rev.RevSeq().InvSeq();
		int EndOfFwdSequence = forward.Index() + forward.Length();
		int BeginningOfReverseSequence = fwd_data.IndexOf(rev);
		return (BeginningOfReverseSequence - EndOfFwdSequence);
	}

	float PrimerCalculator::TemperatureDifference(DataSequence fwd_data, DataSequence rev_data, Primer forward, Primer reverse) const
	{
		float firstTemp  = fwd_data.SubSeq(forward.Index(), forward.Length()).Temperature();
		float secondTemp = rev_data.SubSeq(reverse.Index(), reverse.Length()).Temperature();
		float difference = firstTemp > secondTemp ? (firstTemp - secondTemp) : secondTemp - firstTemp;
		return difference;
	}
} // End of DeGenPrime