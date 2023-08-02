// *********************** class SequenceReader ***********************	//
// Purpose:	class SequenceReader is a class designed to read a file		//
//		and build data classes based on what it reads.					//
// Constructor:	- Default constructor to initialize a blank object.		//
// Methods:	CreateList(ifstream ifs): The main function of the class.	//
//			Finds the ends of the header of the file, then				//
//			begins to read the file, line by line, creating				//
//			Sequence objects of each line that passes a parse			//
//			check.  The Sequences are added to a SequenceList			//
//			object which is returned at the end of the function			//
//		CreateSequence(string): The argument contains both the			//
//			Sequence name and its data.  This function splits			//
//			them and builds a Sequence object.  The name of the			//
//			sequence currently always ends with a ' ' char.				//
//		ParseLine(string): Returns true if the argument is				//
//			able to be converted into a Sequence object.				//
//		FindEndOfHeader(ifstream ifs): Finds the end of the header.		//
//		ReturnToBeginning(ifstream ifs): Sets the read position			//
//			back to the end of the header location.						//
// Private Members: _start_position: The end of the header location.	//
// ******************************************************************** //
#ifndef SEQUENCE_READER
#define SEQUENCE_READER

#include <fstream>
#include <string>
#include "sequence.h"
#include "sequencelist.h"

namespace DeGenPrime
{
	enum FileType
	{
		clust,
		fasta
	};
		
	class SequenceReader
	{
	public:
		SequenceReader();
		
		FileType IdentifyFileType(std::ifstream &ifs);
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