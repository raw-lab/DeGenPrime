// *********************** class DataSequence ***********************
// Purpose:	class DataSequence is a list class containing only
//		a list of DataNodes and designed to perform specific
//		opeations on the list.
// Constructor:	- Default constructor to initialize a blank object
// Mutators: void SetList(vector<DataNode> catalog) standard mutator.
// Functions: PushBack(DataNode node): Appends node to list.
//		  		    PopBack(): Removes last node of list.
// SubSeq(int startIndex, int length): Returns a length list of DataNodes
//						   beginning at the start index.
//						   Return is designed to be used
//						   with SetList to make a new
//						   DataSequence.
// IndecesOfNonNSubsequences(int len): Returns a list of all indeces
//						   where a subsequence of length len
//						   would contain no DataNodes with
//						   code 'N'.
// Accessors: GetDataSequence(): Returns _list.
//				 size(): Returns the size of the list.
//		     AverageRatio(): Returns the average of all ratios in the list.
// Private members: vector<DataNode> list: the list of data

#ifndef DATA_SEQUENCE
#define DATA_SEQUENCE

#include <vector>
#include "DataNode.h"
#include "global.h"

namespace DeGenPrime
{
	class DataSequence
	{
	public:
		DataSequence();
		
		void SetList(std::vector<DataNode> catalog);
		void PushBack(DataNode node);
		void PopBack();

		void Print(int id); // Test function

		DataSequence SubSeq(int startIndex,int length);
		DataSequence InvSeq();
		DataSequence RevSeq();

		float Enthalpy() const;
		float Entropy(float salt_concentration) const;
		float Gibbs(float temperature, float salt_concentration) const;
		
		std::vector<DataNode> GetDataSequence() const;
		int size() const;
		float AverageRatio() const;

		int IndexOf(DataSequence data) const;
		float Temperature() const;
		bool checkMatch(DataSequence data) const;
	private:
		std::vector<DataNode> _list;
	};
} // End of DeGenPrime
#endif // DATA_SEQUENCE