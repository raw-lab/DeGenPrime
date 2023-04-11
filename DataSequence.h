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

		void Print(int id, float temperature, float salt_conc, float mg_conc, float primer_conc); // Test function

		DataSequence SubSeq(int startIndex,int length);
		DataSequence InvSeq();
		DataSequence RevSeq();

		float AverageRatio() const;
		float Enthalpy() const;
		float Entropy() const;
		float Gibbs(float temperature) const;
		float BasicTemperature() const;
		float RTlnK(float temperature, float primer_conc) const;
		float RlnK(float primer_conc) const;
		float AdvancedTemperature(float salt_conc, float mg_conc, float primer_conc) const;
		float BasicAnneal(DataSequence source, int startIndex) const;
		
		std::vector<DataNode> GetDataSequence() const;

		int CountMatches(DataSequence data) const;
		int RevIndex(int index) const;
		int IndexOf(DataSequence data) const;
		int size() const;

		bool checkMatch(DataSequence data) const;
	private:
		std::vector<DataNode> _list;
	};
} // End of DeGenPrime
#endif // DATA_SEQUENCE