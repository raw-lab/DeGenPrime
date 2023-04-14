// DataNode.cpp
#include <iostream>
#include <vector>
#include "DataNode.h"
#include "global.h"

using namespace std;

namespace DeGenPrime
{
	DataNode::DataNode(std::vector<char> char_list)
	{
		int Count[6] = {0,0,0,0,0,0}; // {A,C,G,T,-,N}
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
				case '-':
					Count[4]++;
					break;
				default:
					Count[5]++;
					break;
			}
		}
		ChooseCode(Count, char_list.size());
		int Index = MostCommonIndex(Count);
		switch(Index) // Set Most Common
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
			case 4:
				_most_common = '-';
				break;
			default:
				_most_common = 'N';
				break;
		}
		_ratio = ((float)Count[Index]) / (float)char_list.size(); // Set ratio
		EvaluateCode(); // Double check code
	}

	DataNode::DataNode(char code, char mc, float ratio) // build node directly
	{
		_code = code;
		_most_common = mc;
		_ratio = ratio;
	}

	DataNode DataNode::InvNode()
	{
		char inv_code;
		char inv_most_common;
		switch(_code)
		{
			case 'A':
				inv_code = 'T';
				break;
			case 'C':
				inv_code = 'G';
				break;
			case 'G':
				inv_code = 'C';
				break;
			case 'T':
				inv_code = 'A';
				break;
			case 'M':
				inv_code = 'K';
				break;
			case 'K':
				inv_code = 'M';
				break;
			case 'R':
				inv_code = 'Y';
				break;
			case 'Y':
				inv_code = 'R';
				break;
			case 'B':
				inv_code = 'V';
				break;
			case 'V':
				inv_code = 'B';
				break;
			case 'D':
				inv_code = 'H';
				break;
			case 'H':
				inv_code = 'D';
				break;
			default:
				inv_code = _code;
				break;
		}
		switch(_most_common)
		{
			case 'A':
				inv_most_common = 'T';
				break;
			case 'C':
				inv_most_common = 'G';
				break;
			case 'G':
				inv_most_common = 'C';
				break;
			case 'T':
				inv_most_common = 'A';
				break;
			default:
				inv_most_common = _most_common;
				break;
		}
		return DataNode(inv_code, inv_most_common, _ratio);
	}

	void DataNode::Print()
	{
		cout << "Code: [" << _code;
		cout << "] MC: [" << _most_common;
		cout << "] Ratio: [" << _ratio;
		cout << "%]";
	}
	
	void DataNode::ChooseCode(int Count[6], int Size)
	{
		// This first check is going to automatically assign 'N' as
		// the code if there are more than the argument's ratio
		// of 'N' codes.  This could be improved with a global
		// variable setting.
		
		// {A,C,G,T,-,N}

		// Get the Value of 'N'
		int N = Count[5];
		float ratio = static_cast<float>(N) / static_cast<float>(Size);
		if(ratio > MAX_N_TOLERANCE)
		{
			_code = 'N';
			return;
		}

		// Get the Value of '-'
		N = Count[4];
		ratio = static_cast<float>(N) / static_cast<float>(Size);
		if(ratio > MAX_DASH_VERTICAL_TOLERANCE)
		{
			_code = '-';
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

	void DataNode::EvaluateCode()
	{
		switch(_code)
		{
			case 'A':
			case 'C':
			case 'G':
			case 'T':
			case '-':
				return;
				break;
			default:
				if(1.0 - _ratio >= MIN_DEGENERATE_THRESHOLD)
				{
					return; // In this case, keep the degenerate code
				}
				else
				{
					_code = _most_common; // Not enough degeneracy
				}
				break;
		}
		return;
	}

	int DataNode::MostCommonIndex(int Count[6]) // {A,C,G,T,-,N}
	{
		int index = 0;
		int value = Count[0];
		for(int i = 1;i<6;i++)
		{
			if(value < Count[i])
			{
				value = Count[i];
				index = i;
			}
		}
		return index;
	}

	float DataNode::Enthalpy(DataNode node) const
	{
		// Return zero enthalpy for a '-' most common code.
		// Ethalpy is symmetric (except GC and CG)
		//	- typecast char to int and multiply ascii codes
		//	- run through switch statement to set enthalpy
		float enthalpy = 0.0;
		int node_mc = (node.GetMostCommon() == '-') ? 0 : (int)node.GetMostCommon();
		int this_mc = (GetMostCommon() == '-') ? 0 : (int)GetMostCommon();
		switch(node_mc * this_mc)
		{
			// 'A' = 65, 'C' = 67, 'G' = 71, 'T' = 84
			// Values from: https://www.pnas.org/doi/10.1073/pnas.95.4.1460
			// Values given in units of kcal/mol
			case 4225:	// 'AA'
			case 7056:	// 'TT'
				enthalpy =  -7.9;
				break;
			case 4355:	// 'AC', 'CA'
				enthalpy = -8.5;
				break;
			case 4615:	// 'AG', 'GA'
				enthalpy = -8.2;
				break;
			case 5460:	// 'AT', 'TA'
				enthalpy = -7.2;
				break;
			case 4489:	// 'CC'
			case 5041:	// 'GG'
				enthalpy = -8.0;
				break;
			case 4757:	// 'CG', 'GC'
				enthalpy = (this_mc == 67) ? -10.6 : -9.8;
				break;
			case 5628:	// 'CT', 'TC'
				enthalpy = -7.8;
				break;
			case 5964:	// 'GT', 'TG'
				enthalpy = -8.4;
				break;
			default:	// either contains '-'
				break;
		}
		return enthalpy;
	}

	float DataNode::Entropy(DataNode node) const
	{
		// Return zero entropy for a '-' char
		//	- typecast most common char to int
		//	- multiply this value and run through switch
		float entropy = 0.0;
		int node_mc = (node.GetMostCommon() == '-') ? 0 : (int)node.GetMostCommon();
		int this_mc = (GetMostCommon() == '-') ? 0 : (int)GetMostCommon();
		switch(node_mc * this_mc)
		{
			// 'A' = 65, 'C' = 67, 'G' = 71, 'T' = 84
			// Values from: https://www.pnas.org/doi/10.1073/pnas.95.4.1460
			// Values given in units of cal/(k * mol)
			case 4225:	// 'AA'
			case 7056:	// 'TT'
				entropy =  -22.2;
				break;
			case 4355:	// 'AC', 'CA'
				entropy = -22.7;
				break;
			case 4615:	// 'AG', 'GA'
				entropy = -22.2;
				break;
			case 5460:	// 'AT', 'TA'
				entropy = (this_mc == 65) ? -20.4 : -21.3;
				break;
			case 4489:	// 'CC'
			case 5041:	// 'GG'
				entropy = -19.9;
				break;
			case 4757:	// 'CG', 'GC'
				entropy = (this_mc == 67) ? -27.2 : -24.4;
				break;
			case 5628:	// 'CT', 'TC'
				entropy = -21.0;
				break;
			case 5964:	// 'GT', 'TG'
				entropy = -22.4;
				break;
			default:	// either contains '-'
				break;
		}
		return entropy;
	}

	float DataNode::Gibbs(DataNode node, float temperature) const
	{
		// Delta G = Delta H - Temperature * Delta S
		// Temperature for the Equation is in Celsius, but will be given in celsius
		// Delta S is in cal/K * mol, Delta H is kcal/mol
		float gibbs = Enthalpy(node);
		gibbs -= (temperature + 273.15) * ( Entropy(node) / 1000.0 );

		char node_char = node.GetMostCommon();
		char this_char = GetMostCommon();

		// cout << "Calculating Gibbs of [" << this_char << "-" << node_char << "]: " << gibbs << endl;
		return gibbs;
		
		/* Old Method
		float gibbs = 0.0;

		
		
		int node_mc = (node_char == '-') ? 0 : (int)node_char;
		int this_mc = (this_char == '-') ? 0 : (int)this_char;
		switch(node_mc * this_mc)
		{
			// 'A' = 65, 'C' = 67, 'G' = 71, 'T' = 84
			// Values from: https://www.pnas.org/doi/10.1073/pnas.95.4.1460
			// Values given in units of kcal/(k * mol)
			// Values were taken at 37 degrees celsius.
			case 4225:	// 'AA'
			case 7056:	// 'TT'
				gibbs =  -1.0;
				break;
			case 4355:	// 'AC', 'CA'
				gibbs = -1.45;
				break;
			case 4615:	// 'AG', 'GA'
				gibbs = -1.3;
				break;
			case 5460:	// 'AT', 'TA'
				gibbs = (this_mc == 65) ? -0.88 : -0.58;
				break;
			case 4489:	// 'CC'
			case 5041:	// 'GG'
				gibbs = -1.84;
				break;
			case 4757:	// 'CG', 'GC'
				gibbs = (this_mc == 67) ? -2.17 : -2.24;
				break;
			case 5628:	// 'CT', 'TC'
				gibbs = -1.28;
				break;
			case 5964:	// 'GT', 'TG'
				entropy = -1.44;
				break;
			default:	// either contains '-'
				break;
		}

		cout << gibbs << endl;
		
		return gibbs;
		*/
	}

	char DataNode::GetCode() const { return _code; }
	char DataNode::GetMostCommon() const { return _most_common; }
	float DataNode::Ratio() const { return _ratio; }
} // End of DeGenPrime