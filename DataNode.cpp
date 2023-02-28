// DataNode.cpp
#include <iostream>
#include <vector>
#include "DataNode.h"

using namespace std;

namespace DeGenPrime
{
	DataNode::DataNode(std::vector<char> char_list)
	{
		int Count[5] = {0,0,0,0,0}; // {A,C,G,T,N}
		for(char c : char_list)
		{
			switch(c)
			{
				case 'a':
				case 'A':
					Count[0]++;
					break;
				case 'c':
				case 'C':
					Count[1]++;
					break;
				case 'g':
				case 'G':
					Count[2]++;
					break;
				case 't':
				case 'T':
					Count[3]++;
					break;
				default:
					Count[4]++;
					break;
			}
		}
		ChooseCode(Count, char_list.size());
		int Index = MostCommonIndex(Count);
		switch(Index)
		{
			case 0:
				_most_common = 'A';
				break;
			case 1:
				_most_common = 'C';
				break;
			case 2:
				_most_common = 'G';
				break;
			case 3:
				_most_common = 'T';
				break;
			default:
				_most_common = 'N';
				break;
		}
		_ratio = ((float)Count[Index]) / (float)char_list.size();
	}

	// Test Function
	void DataNode::Print()
	{
		cout << "Code: [" << _code;
		cout << "] MC: [" << _most_common;
		cout << "] Ratio: [" << _ratio;
		cout << "%] Weighted: [" << WeightedRatio();
		cout << "]";
	}
	
	void DataNode::ChooseCode(int Count[5], int Size)
	{
		// This first check is going to automatically assign 'N' as
		// the code if there are more than the argument's ratio
		// of 'N' codes.  This could be improved with a global
		// variable setting.
		
		// Get the Value of 'N'
		int N = Count[4];
		float ratio = (float)N / (float)Size;
		if(ratio >= 0.1)
		{
			_code = 'N';
			return;
		}
		
		if(Count[0] != 0) // There is at least one 'A'
		{
			if(Count[1] != 0) // There is at least one 'C'
			{
				if(Count[2] != 0) // There is at least one 'G'
				{
					if(Count[3] != 0) // There is at least one 'T'
					{
						_code = 'N'; // A,C,G,T
					}
					else // No 'T'
					{
						_code = 'V'; // A,C,G
					}
				}
				else // No 'G'
				{
					if(Count[3] != 0) // There is at least one 'T'
					{
						_code = 'H'; // A,C,T
					}
					else // No 'T'
					{
						_code = 'M';
					}
				}
			}
			else // No 'C'
			{
				if(Count[2] != 0) // There is at least one 'G'
				{
					if(Count[3] != 0) // There is at least one 'T'
					{
						_code = 'D'; // A,G,T
					}
					else // No 'T'
					{
						_code = 'R'; // A,G
					}
				}
				else // No 'G'
				{
					if(Count[3] != 0) // There is at least one 'T'
					{
						_code = 'W'; // A,T
					}
					else // No 'T'
					{
						_code = 'A';
					}
				}
			}
		}
		else // No 'A'
		{
			if(Count[1] != 0) // There is at least one 'C'
			{
				if(Count[2] != 0) // There is at least one 'G'
				{
					if(Count[3] != 0) // There is at least one 'T'
					{
						_code = 'B'; // C,G,T
					}
					else // No 'T'
					{
						_code = 'S'; // C,G
					}
				}
				else // No 'G'
				{
					if(Count[3] != 0) // There is at least one 'T'
					{
						_code = 'Y'; // C,T
					}
					else // No 'T'
					{
						_code = 'C';
					}
				}
			}
			else // No 'C'
			{
				if(Count[2] != 0) // There is at least one 'G'
				{
					if(Count[3] != 0) // There is at least one 'T'
					{
						_code = 'K'; // G,T
					}
					else // No 'T'
					{
						_code = 'G';
					}
				}
				else // No 'G'
				{
					if(Count[3] != 0) // There is at least one 'T'
					{
						_code = 'T';
					}
					else // No 'T'
					{
						_code = 'N'; // Default Case
					}
				}
			}
		}
		return;
	}

	int DataNode::MostCommonIndex(int Count[5])
	{
		int index = 0;
		int value = Count[0];
		for(int i = 1;i<5;i++)
		{
			if(value < Count[i])
			{
				value = Count[i];
				index = i;
			}
		}
		return index;
	}

	char DataNode::GetCode() const { return _code; }
	char DataNode::GetMostCommon() const { return _most_common; }
	float DataNode::Ratio() const { return _ratio; }

	float DataNode::WeightedRatio() const
	{
		float weighted = _ratio;
		switch(_code)
		{
			case 'N':
				weighted *= weighted;
			case 'D':
			case 'V':
			case 'B':
			case 'H':
				weighted *= weighted;
			case 'R':
			case 'Y':
			case 'M':
			case 'K':
			case 'S':
			case 'W':
				weighted *= weighted;
				break;
			default:
				break;
		}
		return weighted;
	}

} // End of DeGenPrime