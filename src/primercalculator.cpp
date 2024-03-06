#include <algorithm>
#include <cctype>
#include <cstring>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "primer.h"
#include "sequencelist.h"
#include "datasequence.h"
#include "primercalculator.h"
#include "globalsettings.h"
#include "format.h"
#include "global.h"

using namespace std;

namespace DeGenPrime
{
	PrimerCalculator::PrimerCalculator() { }

	void PrimerCalculator::InitializeTestPrimer(DataSequence data)
	{
		Primer test(0, data.size());
		test.SetPenalty(data.Penalty());
		PushBack(test);
		_OriginalSize = size();
	}

	void PrimerCalculator::InitializePrimers(DataSequence data)
	{
		int endIndex = data.size() - GlobalSettings::GetMinimumAmplicon();

		// Check AmpliconLength is valid
		if(endIndex < 0)
		{
			endIndex = 0;
			GlobalSettings::SetMinimumAmplicon(data.size());
		}

		for(int i = GlobalSettings::GetMinimumPrimerLength();i <= GlobalSettings::GetMaximumPrimerLength();i++)
		{
			for(int j = 0;j <= endIndex;j++)
			{
				Primer id(j, i);
				DataSequence sub = data.SubSeq(j,i);
				float penalty = sub.Penalty();
				if(penalty > MAX_PENALTY)
				{
					continue;
				}
				else
				{
					id.SetPenalty(penalty);
					PushBack(id);
				}
			}
		}
		_OriginalSize = size();
	}

	void PrimerCalculator::InitializeBoundedPrimers(DataSequence data, int lowerBound)
	{
		int loopBound = lowerBound;

		// Check the loopBound is appropriate
		loopBound = loopBound > 0 ? loopBound : 0;
		int min = GlobalSettings::GetMinimumPrimerLength();
		int max = GlobalSettings::GetMaximumPrimerLength();
		if(loopBound + max > data.size())
		{
			loopBound = data.size() - max;
		}
		
		// Needs to loop through all primers in the forward direction
		for(int i = min;i <= max;i++)
		{
			for(int j = 0;j <= loopBound;j++)
			{
				Primer id(j, i);
				DataSequence sub = data.SubSeq(j,i);
				float penalty = sub.Penalty();
				if(penalty > MAX_PENALTY)
				{
					continue;
				}
				else
				{
					id.SetPenalty(penalty);
					PushBack(id);
				}
			}
		}
		_OriginalSize = size();
	}

	void PrimerCalculator::InitializeFromRegion(std::vector<Primer> region, DataSequence data)
	{
		int min = GlobalSettings::GetMinimumPrimerLength();
		int max = GlobalSettings::GetMaximumPrimerLength();
		const int threshold = max - min;
		for(Primer p : region)
		{
			int region_size = p.Length();
			for(int i = min;i <= max;i++)
			{
				int endIndex = (region_size - i < threshold) ? region_size - i : threshold;
				/*
				bool search = (GlobalSettings::GetSearchFwd() || GlobalSettings::GetSearchRev());
				if(search)endIndex = region_size - i;*/
				for(int j = 0;j <= endIndex;j++)
				{
					Primer pr(p.Index() + j, i);
					DataSequence sub = data.SubSeq(pr.Index(), pr.Length());
					float penalty = sub.Penalty();
					if(penalty > MAX_PENALTY)
					{
						continue;
					}
					else
					{
						pr.SetPenalty(penalty);
						PushBack(pr);
					}
				}
			}
		}
		_OriginalSize = size();
	}

	void PrimerCalculator::Sort()
	{
		sort(_primers.begin(), _primers.end());
	}

	string PrimerCalculator::FilterAll(DataSequence data)
	{
		string ret = "";
		ret += FilterMessage("before", 0);
		ret += FilterMessage("start", 0);
		ret += FilterDegeneracy(data);
		ret += FilterDeletions(data);
		ret += FilterGCContent(data);
		ret += FilterRepeats(data);
		ret += FilterComplementaryEnds(data);
		ret += FilterTemperature(data, 0.0);
		ret += FilterMessage("final", _OriginalSize - size());
		return ret;
	}

	string PrimerCalculator::FilterDegeneracy(DataSequence data)
	{
		string ret = "";
		int filtercount = 0;
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
						break;
					case 'B':
					case 'D':
					case 'V':
					case 'H':
						degeneracy_count += (2 * position_multiplier);
						break;
					case 'A':
					case 'C':
					case 'G':
					case 'T':
					case '-':
						break;
					default:
						degeneracy_count += position_multiplier;
						break;
				}
			}
			if(degeneracy_count >= 3)
			{
				Erase(i);
				filtercount++;
			}
		}
		ret = FilterMessage("FilterDegeneracy", filtercount);
		return ret;
	}

	string PrimerCalculator::FilterDeletions(DataSequence data)
	{
		string ret = "";
		int filtercount = 0;
		for(int i = _primers.size() - 1; i >= 0;i--)
		{
			DataSequence p = data.SubSeq(_primers[i].Index(), _primers[i].Length());
			bool flag = false;
			int total_deletions_count = 0;
			// Filter sequences for deletions
			for(int j = 0;j < p.size();j++)
			{
				if(flag)
				{
					break;
				}
				if(p.GetDataSequence()[j].GetCode() == '-')
				{
					total_deletions_count++;
					if(j < 3 || j > (p.size() - 4))
					{
						flag = true;
					}
				}
				if(total_deletions_count > 3)
				{
					flag = true;
				}
				if((p.size() - total_deletions_count) < GlobalSettings::GetMinimumPrimerLength())
				{
					flag = true;
				}
			}
			if(flag)
			{
				Erase(i);
				filtercount++;
			}
		}
		ret += FilterMessage("FilterDeletions", filtercount);
		return ret;
	}

	string PrimerCalculator::FilterGCContent(DataSequence data)
	{
		string ret = "";
		int filtercount = 0;
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
					if(j > (p.size() - 6))
					{
						gc_end_count++;
					}
				}
			}
			float end_ratio = ((float)(gc_end_count) / 5.0);
			float total_ratio = (float)(gc_total_count) / (float)(p.size());
			bool ending = end_ratio > MAX_GC_EXTREMA_RATIO;
			bool total = (total_ratio < MIN_GC_TOTAL_RATIO || total_ratio > MAX_GC_TOTAL_RATIO);
			if(ending || total)
			{
				Erase(i);
				filtercount++;
			}
		}
		ret = FilterMessage("FilterGCContent", filtercount);
		return ret;
	}

	string PrimerCalculator::FilterRepeats(DataSequence data)
	{
		string ret;
		int filtercount = 0;
		for(int i = _primers.size() - 1; i >= 0;i--)
		{
			// Define the primer
			DataSequence p = data.SubSeq(_primers[i].Index(), _primers[i].Length());
			bool flag = false;
			// Filter sequences with too much repition of nucleotides
			//	- raise flag if four matches are found to any pair, more than one codon, or any higher
			//	- break loop if flag or length * (match count) < remaining checks * match count
			//	- there is no need to check subsequences longer than sqrt the total length.
			for(int j = 0;j < p.size() - 6;j++)
			{
				if(flag)
				{
					break;
				}
				// Define pattern matching subsequences.
				DataSequence two = p.SubSeq(j,2);
				DataSequence three = p.SubSeq(j,3);
				DataSequence four = p.SubSeq(j,4);

				// Define match counters
				int double_match_count = 0;
				int triple_match_count = 0;
				int quadra_match_count = 0;

				// Choose k starting at index j + 1 and run through size - 4
				int k = j + 1;
				while(k < p.size() - 5)
				{
					if(flag)
					{
						break;
					}
					if(two.GetDataSequence()[j].GetMostCommon() == 
						p.GetDataSequence()[k].GetMostCommon())
					{
						double_match_count += two.checkMatch(p.SubSeq(k,2)) ? 1 : 0;
						triple_match_count += three.checkMatch(p.SubSeq(k,3)) ? 1 : 0;
						quadra_match_count += four.checkMatch(p.SubSeq(k,4)) ? 1 : 0;
						flag = double_match_count > TooManyRepeats(2) ||
							triple_match_count > TooManyRepeats(3) ||
							quadra_match_count > TooManyRepeats(4);
					}
					k++;
				}
			}
			if(flag)
			{
				Erase(i);
				filtercount++;
			}
		}
		ret = FilterMessage("FilterRepeats", filtercount);
		return ret;
	}

	string PrimerCalculator::FilterComplementaryEnds(DataSequence data)
	{
		string ret;
		int filtercount = 0;
		for(int i = _primers.size() - 1;i >= 0;i--)
		{
			DataSequence p = data.SubSeq(_primers[i].Index(), _primers[i].Length());
			DataSequence firstThree = p.SubSeq(0,3);
			DataSequence lastThree = p.SubSeq(p.size() - 3,3);
			lastThree = lastThree.InvSeq().RevSeq();

			if(firstThree.checkMatch(lastThree))
			{
				Erase(i);
				filtercount++;
			}
		}
		ret = FilterMessage("FilterCompEnds", filtercount);
		return ret;
	}

	string PrimerCalculator::FilterHairpins(DataSequence data)
	{
		string ret;
		int filtercount = 0;
		for(int i = _primers.size() - 1;i >= 0;i--)
		{
			DataSequence p = data.SubSeq(_primers[i].Index(), _primers[i].Length());
			bool flag = false;
			
			// Filter all sequences that would form hairpins loops
			// -	for triloops,
			//	compare the code five nucleotides away for complement
			//	if match, check 2nd nucleotide for 'G'
			//	check 2nd to last nucleotide for 'A'
			//	filter if those conditions are met
			// -	for tetraloops, compare the nucleotide six nucleotides away for complement
			//	if match filter

			int k_mer = 4;
			while(k_mer < ((p.size() - 3.0)/2.0))
			{
				if(flag)
				{
					break;
				}
				for(int j = 0;j + 3 < p.size() - k_mer;j++)
				{
					if(flag)
					{
						break;
					}
					
					DataSequence first = p.SubSeq(j,k_mer);
					if(k_mer == 4 && first.GCRatio() < 1.0)
					{
						continue;
					}
					for(int k = j + k_mer + 3;k < p.size() - k_mer;k++)
					{
						if(flag)
						{
							break;
						}
						DataSequence second = p.SubSeq(k,k_mer);
						second = second.InvSeq().RevSeq();
						int matches = first.CountMatches(second);
						if((k_mer == 4 || k_mer == 5) && k_mer == matches)
						{
							flag = true;
						}
						else if((k_mer == 6 || k_mer == 7) && matches >= 5)
						{
							flag = true;
						}
						else if((k_mer == 8 || k_mer == 9) && matches >= 6)
						{
							flag = true;
						}
						else if(matches >= 7)
						{
							flag = true;
						}
					}
				}
				k_mer++;
			}
			if(flag)
			{
				Erase(i);
				filtercount++;
			}
		}
		ret = FilterMessage("FilterHairpins", filtercount);
		return ret;
	}

	string PrimerCalculator::FilterDimers(DataSequence data)
	{
		string ret;
		int filtercount = 0;
		float value = 0.0;
		for(int i = _primers.size() - 1;i >= 0;i--)
		{
			DataSequence p = data.SubSeq(_primers[i].Index(), _primers[i].Length());

			// Filter all sequences that would for self-dimers
			// -	use the 3' end of the primer (last 5 nucleotides) as a reference check
			// -	self dimers are most likely to form when these nucleotides have too low of enthalpy
			// -	filter primers with 3' end < -3.0

			// Pull Deletions off the sequence
			for(int j = 0;j < p.size();j++)
			{
				if(p.GetDataSequence()[j].GetMostCommon() == '-')
				{
					p.Erase(j);
				}
			}
			DataSequence ending = (p.size() > 5) ? p.SubSeq(p.size() - 6, 5) : p;
			value = ending.isEmpty() ? 0 : ending.Gibbs();
			if(value < GlobalSettings::GetDeltaG())
			{
				Erase(i);
				filtercount++;
			} 
		}
		ret = FilterMessage("FilterDimers", filtercount);
		return ret;
	}

	string PrimerCalculator::FilterTemperature(DataSequence data, float offset)
	{
		string ret;
		int filtercount = 0;
		for(int i = _primers.size() - 1; i >= 0;i--)
		{
			DataSequence p = data.SubSeq(_primers[i].Index(), _primers[i].Length());
			if(p.NNMeltingTemperature() < (GlobalSettings::GetMinimumTemperature() + offset/2.0) || 
				p.NNMeltingTemperature() > (GlobalSettings::GetMaximumTemperature() - offset/2.0))
			{
				Erase(i);
				filtercount++;
			}
		}
		ret = FilterMessage("FilterTemp", filtercount);
		return ret;
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
				b = 4;
				break;
			case 3:
				b = 3;
				break;
			default:
				b = 2;
				break;
		}
		return b;
	}

	void PrimerCalculator::PrintSize()
	{
		cout << "The number of primers in the list is: " << size() << endl;
	}

	std::string PrimerCalculator::PrintAll()
	{
		string ret = "";
		for(Primer p : _primers)
		{
			ret += p.Print();
		}
		return ret;
	}

	void PrimerCalculator::SetPrimers(std::vector<Primer> primerList) { _primers = primerList; }

	std::vector<Primer> PrimerCalculator::GetPrimers() const { return _primers; }
	std::string PrimerCalculator::FilterMessage(std::string func, int filtercount)
	{
		bool start = func == "start";
		bool before = func == "before";
		bool final = func == "final";
		string ret = "";
		string line = "";
		if(start)
		{
			line = "Filter Name";
			ret += Format(line, 25, Alignment::Center);
			line = "Filtered  ";
			ret += Format(line, 25, Alignment::Right);
			line = "   Percent";
			ret += Format(line, STR_FORMAT - 50, Alignment::Left) + "\n";
			return ret;
		}
		if(before)
		{
			ret = "Before filters [" + to_string(size());
			ret += "] total primers.\n\n";
			return ret;
		}
		line += final ? "After all filters" : ("     " + func);
		ret += Format(line, 25, Alignment::Left);
		line = to_string(filtercount) + "   ";
		ret += Format(line, 25, Alignment::Right);
		line = "(";
		float percent = 100.0 * ((float)filtercount)/((float)_OriginalSize);
		line += Format(percent, 2);
		line += "% of total.)";
		ret += Format(line, STR_FORMAT - 50, Alignment::Left) + "\n";
		if(final)
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
			percent = 100.0 * ((float)size())/((float)_OriginalSize);
			line += Format(percent,2) + "% of total.)";
			ret += Format(line, STR_FORMAT - 50, Alignment::Left) + "\n";
		}
		return ret;
	}

	int PrimerCalculator::IndexOf(DataSequence data, string str) const
	{
		int index = -1;
		bool flag = false;

		// Build a DataSequence of the chars in str
		DataSequence s;
		for(char c : str)
		{
			DataNode node(c,c,1.0);
			s.PushBack(node);
		}

		for(int i = 0;i < _primers.size();i++)
		{
			if(flag)
			{
				break;
			}
			DataSequence p = data.SubSeq(_primers[i].Index(), _primers[i].Length());
			if(p.checkMatch(s))
			{
				index = i;
				flag = true;
			}
		}
		return index;
	}

	int PrimerCalculator::size() const { return _primers.size(); }
} // End of DeGenPrime