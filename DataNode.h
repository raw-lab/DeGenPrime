// ************************** class DataNode **************************
// Purpose: Class DataNode is a class containing data from
//		a SequenceList at a particular intersection
//		across all Sequences in the list.
// Constructor:	This object can be constructed from the
//			vector of chars across the intersection.
// Acesssors:	char GetCode() - returns _code
//			char GetMostCommon() - returns _most_common.
//			float Ratio() - returns _percentage.
//			float WeightedRatio()
//				returns a weighted ratio based on
//				the code and the ratio.  It will
//				always be <= the _ratio value.
//				For Codes(A,T,C,G):
//				  weighted = Ratio^1 (same value)
//				For Codes(R,Y,M,K,S,W):
//				  weighted = Ratio^2
//				For Codes(H,B,V,D):
//				  weighted = Ratio^3
//				For Code (N): weighted = Ratio^4
//				This will produce a lower weighted score
//				for codes that have lower ratio.
//				Example: DataNode[0] is Code S with
//					   9 'C's and 1 'G'. (ratio = .9)
//					   Weighted = .9^2 = .81
//					   DataNode[1] is Code S with
//					   5 'C's and 5 'G's. (ratio = .5)
//					   Weighted = .5^2 = .25
//				In this example, DataNode 0 and 1 both have the 
//				same code, but the more consistent DataNode has
//				a much higher weighted score.
// Private Methods: ChooseCode(int Count[5], int Size)
//			Called in constuctor to choose a Code
//			  MostCommonIndex(int Count[5])
//			  In the case of two codes being equally common
//			  (unlikely in large data lists) it will return the
//			  lowest index.
//			Called in constructor. Returns the index of
//			the greatest value in the ine argument.
// Private Members: char _code: The chosen code of the DataNode.
//			  char _most_common: The most common code from
//				the data used to construct the object.
//			  float _ratio: This is a ratio of the
//				occurence of the most common char in 
//				the list to the size of the list.
//				Range: 0.0 - 1.0
// Potential
// Improvements:	-	Adding a Global Threshold variable setting for the
//				ratio of 'N' across a node to make it count as N.
//			-	Implement an enum for SequenceType as a second
//				argument to the constructor would allow this class
//				to process protein codes as well.
// ***********************************************************************
#ifndef DATA_NODE
#define DATA_NODE

#include <vector>

namespace DeGenPrime
{	
	class DataNode
	{
	public:
		DataNode(std::vector<char> char_list);

		void Print(); // Test Function

		char GetCode() const;
		char GetMostCommon() const;
		float Ratio() const;
		float WeightedRatio() const;
	private:
		void ChooseCode(int Count[5], int Size);
		int MostCommonIndex(int Count[5]);
		char _code;
		char _most_common;
		float _ratio;
	};
} // End of DeGenPrime
#endif // DATA_NODE