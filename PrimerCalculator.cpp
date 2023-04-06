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

	void PrimerCalculator::InitializeBoundedPrimers(DataSequence data, int lowerBound)
	{
		// Check the lowerBound is appropriate
		int loopBound = lowerBound > 0 ? lowerBound : 0;
		if(loopBound + MAX_PRIMER_LENGTH > data.size())
		{
			loopBound = data.size() - MAX_PRIMER_LENGTH;
		}
		
		// Needs to loop through all primers in the forward direction
		for(int i = MIN_PRIMER_LENGTH;i <= MAX_PRIMER_LENGTH;i++)
		{
			for(int j = 0;j < loopBound;j++)
			{
				Primer id(j,i);
				PushBack(id);
			}
		}
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

	void PrimerCalculator::FilterHairpins(DataSequence data)
	{
		for(int i = _primers.size() - 1;i >= 0;i--)
		{
			DataSequence p = data.SubSeq(_primers[i].Index(), _primers[i].Length());
			bool flag = false;
			int match_count = 0;

			// Filter all sequences that would form Hairpins.
			// -	define a fold point within range <3, size - 4> (j)
			// -	split primer into two subsequences
			// -  determine which end of the primer is closest to the fold point
			// -	reverse and invert the end sequence
			// -	count matches between the subsequences
			// -	if matches >= sqrt(sub_length) raise flag, reject primer.
			//		(usually sqrt(sub_length) is 2.xxx, but can be 3 for larger primers.)

			for(int j = 3;j < p.size() - 3;j++)
			{
				if(flag)
				{
					break;
				}
				
				DataSequence begin = p.SubSeq(0, j - 1);
				DataSequence end = p.SubSeq(j + 1, p.size() - j - 1);
				int size_difference = end.size() - begin.size();
				if(size_difference != 0)
				{
					if(size_difference > 0)
					{
						// Get first chars of end sequence
						end = end.SubSeq(0,begin.size());
					}
					else
					{
						// Get last chars of begin sequence
						begin = begin.SubSeq(begin.size() - end.size(), end.size());
					}
				}

				end = end.RevSeq().InvSeq();
				match_count = begin.CountMatches(end);
				if((float)match_count >= sqrt(end.size()))
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

	void PrimerCalculator::FilterDimers(DataSequence data)
	{
		for(int i = _primers.size() - 1;i >= 0;i--)
		{
			DataSequence p = data.SubSeq(_primers[i].Index(), _primers[i].Length());
			bool flag = false;
			int match_count = 0;

			// Filter all sequences that would for self-dimers
			// -	use the 3' end of the primer (first three nucleotides) as a reference check
			// -	self dimers are most likely to form when these nucleotides have matches.
			// -	slide these nucleotides backward along the InvSeq (not RevSeq) of the primer
			// -	start at index p.size() - 4.
			// -	if there are 2 matches, filter the primer

			DataSequence begin = p.SubSeq(0,3);
			DataSequence inv_seq = p.InvSeq();

			for(int j = inv_seq.size() - 4;j >= 0;j--)
			{
				if(flag)
				{
					break;
				}

				DataSequence end = inv_seq.SubSeq(j,3);
				match_count = begin.CountMatches(end);

				flag = (match_count >= 2) ? true : false;
			}
			if(flag)
			{
				Erase(i);
			} 
		}

		/* Old system, was filtering way too much
		for(int i = _primers.size() - 1;i >= 0;i--)
		{
			bool flag = false;
			int match_count = 0;
			// Filter all sequences that would form Self Dimers
			//	- more than 3 matches would likely form self dimer
			//	- Begin with subsequence of length 4 (length 3 covered by FilterComplementaryEnds())
			//	- Slide Window comparing nucleotide matches between Fwd.InvSeq() and Rev.RevSeq()
			//	- Compare forward sliding along rev AND reverse sliding along forward
			DataSequence p = data.SubSeq(_primers[i].Index(), _primers[i].Length());
			DataSequence r = p.InvSeq();
			
			// Slide end of reverse along beginning of forward
			for(int j = 4;j < p.size();j++)
			{
				if(flag)
				{
					break;
				}
				DataSequence p_sub = p.SubSeq(0,j);
				DataSequence r_sub = r.SubSeq(r.size() - j,j);
				match_count = p_sub.CountMatches(r_sub);
				flag = (match_count > 3) ? true : false;
			}
			
			// Slide end of forward along beginning of reverse
			for(int j = 4;j < p.size();j++)
			{
				if(flag)
				{
					break;
				}
				DataSequence p_sub = p.SubSeq(p.size() - j,j);
				DataSequence r_sub = r.SubSeq(0,j);
				match_count = p_sub.CountMatches(r_sub);
				flag = (match_count > 3) ? true : false;
			}
			if(flag)
			{
				Erase(i);
			}
		}
		*/
	}

	void PrimerCalculator::FilterTemperature(DataSequence data, float offset)
	{
		for(int i = _primers.size() - 1; i >= 0;i--)
		{
			DataSequence p = data.SubSeq(_primers[i].Index(), _primers[i].Length());
			if(p.BasicTemperature() < (MIN_PRIMER_TEMP + offset/2.0) || p.BasicTemperature() > (MAX_PRIMER_TEMP - offset/2.0))
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
} // End of DeGenPrime