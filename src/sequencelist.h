// *************************** class SequenceList *****************************	//
// Purpose: class SequenceList is a list class containing only a list of		//
//		sequence objects and designed to perform specific operations			//
//		on the list.															//
// Constructor:	- Default constructor to initialize a blank object.				//
// Mutators: SetList(vector<Sequence> catalog): standard mutator.				//
// Functions: PushBack(Sequence seq): This function checks the name of			//
//			the sequence argument.  If there is already a sequence				//
//			in the list with that name.  It appends the chars in				//
//			the argument's codes to the sequence already in the					//
//			list.  If there is not a sequence with that name in					//
//			the list.  It appends the Sequence to the list.						//
//		  Erase(int index): Erases the Sequence in the list at index.			//
//		  PopBack(): This removes the last sequence in the list.				//
// 		  ProcessList():	Condense the longitudinal data across all 			//
//					the stringsin the list into a data sequence.				//
//					It Uses an iterator to check move through 					//
//					every element in a column of a list and 					//
//					stores it's char value into a char list.  					//
//					The char list is used to generate a DataNode.  				//
//					Pushes the data note onto a data sequence 					//
//					and return the value after cycling through 					//
//					the entire sequence list.									//
//	   CharsAt(int index):	Returns the list of chars at the given index		//
//					of all sequences in the list.								//
// Accessors:	GetSequenceList(): Returns the list								//
//				     size(): Returns the size of the list.						//
// Private members: vector<Sequence> _list: the list of Sequences.				//
// ****************************************************************************	//
#ifndef SEQUENCE_LIST
#define SEQUENCE_LIST

#include <iostream>
#include <vector>
#include "sequence.h"
#include "datasequence.h"

namespace DeGenPrime
{
	class SequenceList
	{
	public:
		SequenceList();

		SequenceList InvRevList();
		DataSequence ProcessList();
		
		void SetList(std::vector<Sequence> catalog);
		void Clear();
		void Erase(int index);
		void PushBack(Sequence seq); // Do not use reg push_back
		void PopBack();	// Do not use reg pop_back
		void Sort();
		SequenceList FilterDashes();
		void RemoveDashes();

		std::string PrintSequenceNames() const;
		std::string DecodeProteins() const;
		std::string CreateFasta() const;
		std::string Section(int index, int length) const;

		bool TestAlignment() const;

		std::vector<char> CharsAt(int index) const;
		std::vector<Sequence> GetSequenceList() const;
		int size() const;
		int IndexOf(std::string name);
	private:
		std::vector<Sequence> _list;
		std::string Codon(char c) const;
	};

} // End of DeGenPrime
#endif // SEQUENCE_LIST