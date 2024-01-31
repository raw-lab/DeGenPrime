// **************************** datanode.h ****************************
// Purpose: Class DataNode is a class containing data from
//		a SequenceList at a particular intersection
//		across all Sequences in the list.
// Constructors:	This object can be constructed from the
//			vector of chars across the intersection.  It can
//			also be constructed from representative char code and mc
//			with float ratio.
// Mutators: InvNode(): switches the code and mc of the node to
//				its inverse.  A -- T and C -- G are inverses.
// Acesssors:	char GetCode() - returns _code
//			char GetMostCommon() - returns _most_common.
//			float Ratio() - returns _percentage.
// Enthalpy(DataNode): Returns the nearest neighbor enthalpy between
//		caller and the argument according to values from:
//		https://www.pnas.org/doi/10.1073/pnas.95.4.1460
// Entropy(DataNode): Returns the nearest neighbor entropy between
//		caller and the argument according to values from:
//		https://www.pnas.org/doi/10.1073/pnas.95.4.1460
// Gibbs(DataNode): Returns the nearest neighbor gibbs between
//		caller and the argument according to values from:
//		https://www.pnas.org/doi/10.1073/pnas.95.4.1460
// Private Methods: ChooseCode(int Count[5], int Size)
//			Called in constuctor to choose a Code
//			  MostCommonIndex(int Count[5]).
//			Called in constructor. Returns the index of
//			the greatest value in the ine argument.
//			EvaluateCode() runs a few additional parameters
//				to make sure the chosen code is correct.
// Private Members: char _code: The chosen code of the DataNode.
//			  char _most_common: The most common code from
//				the data used to construct the object.
//			  float _ratio: This is a ratio of the
//				occurence of the most common char in 
//				the list to the size of the list.
//				Range: 0.0 - 1.0
// ***********************************************************************
#ifndef DATA_NODE
#define DATA_NODE

#include <string>
#include <vector>

namespace DeGenPrime
{	
	class DataNode
	{
	public:
		DataNode(std::vector<char> char_list);
		DataNode(char code, char mc, float ratio);

		DataNode InvNode();
		void Print(); // Test Function
		std::string NodeInfo();

		float Enthalpy(DataNode node) const;
		float Entropy(DataNode node) const;
		float Gibbs(DataNode node) const;

		char GetCode() const;
		char GetMostCommon() const;
		float Ratio() const;
	private:
		void ChooseCode(int Count[6], int Size);
		void EvaluateCode();
		int MostCommonIndex(int Count[5]);

		char _code;
		char _most_common;
		float _ratio;
	};
} // End of DeGenPrime
#endif // DATA_NODE