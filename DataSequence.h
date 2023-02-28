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
//		  AverageWeighted(): Returns the average weighted ratio of the list.
// Private members: vector<DataNode> list: the list of data
// Potential
// Improvements:	- SubSeq might be able to return another DataSequence object.
//			- This class needs to implement a way to check for repetitive
//			  codes within the subsequence.

#ifndef DATA_SEQUENCE
#define DATA_SEQUENCE

#include <vector>
#include "DataNode.h"

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
		void PrintNonNSequences(); // Test Function

		std::vector<DataNode> SubSeq(int startIndex,int length);
		std::vector<int> IndecesOfNonNSubsequences(int length);
		
		std::vector<DataNode> GetDataSequence() const;
		int size() const;
		float AverageRatio() const;
		float AverageWeighted() const;
	private:
		std::vector<DataNode> _list;
	};
} // End of DeGenPrime
#endif // DATA_SEQUENCE