// ********************** class SequenceReader **********************
// Purpose:	class SequenceReader is a class designed to read a file
//		and build data classes based on what it reads.
// Constructor:	- Default constructor to initialize a blank object.
// Methods:	CreateList(ifstream ifs): The main function of the class.
//			Finds the ends of the header of the file, then
//			begins to read the file, line by line, creating
//			Sequence objects of each line that passes a parse
//			check.  The Sequences are added to a SequenceList
//			object which is returned at the end of the function
//		CreateSequence(string): The argument contains both the
//			Sequence name and its data.  This function splits
//			them and builds a Sequence object.  The name of the
//			sequence currently always ends with a ' ' char.
//		ParseLine(string): Returns true if the argument is
//			able to be converted into a Sequence object.
//		FindEndOfHeader(ifstream ifs): Finds the end of the header.
//		ReturnToBeginning(ifstream ifs): Sets the read position
//			back to the end of the header location.
// Private Members: _start_position: The end of the header location.
// Possible
// Improvements:	- Currently this reader is only designed to read the
//			  .clust file format.  Other file formats could be
//			  implemented.
//			- The typecast in Parseline may be misreading some
//			  .clust data newlines.
//			- The ' ' is always saved at the end of Sequence names
//			- A global variable that flags if the first Sequence
//			  is a reference sequence could be established.
// ***********************************************************************
#ifndef SEQUENCE_READER
#define SEQUENCE_READER

#include <fstream>
#include <string>
#include "Sequence.h"
#include "SequenceList.h"

namespace DeGenPrime
{
		
	class SequenceReader
	{
	public:
		SequenceReader();
		SequenceList CreateList(std::ifstream &ifs);
		
	private:
		int _start_position;

		Sequence CreateSequence(std::string line);
		bool ParseLine(std::string line);

		void FindEndOfHeader(std::ifstream &ifs);
		void ReturnToBeginning(std::ifstream &ifs);
	};
} // DeGenPrime end
#endif // SEQUENCE_READER