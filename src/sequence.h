// *********************** class Sequence *************************	//
// Purpose: A class to directly store data read from a				//
//		file into meaningful data.									//
// Constructors:	- Default and String argument to				//
//			  initialize the name of the sequence.					//
// Mutators: SetName(string name): set name of sequence.			//
//		 SetList(vector<char> list): set the list of chars.			//
// Functions: PushBack(char c): Append char to the list.			//
//		  PopBack(): Remove the last char of the list.				//
// Accessors: GetName(): Returns the name							//
//		  GetCodes(): Returns the list of chars.					//
//		  size(): Returns the size of the list of chars.			//
// Private members: string _name: name of the list.					//
//			  vector<char> _codes: the list of chars.				//
// ****************************************************************	//
#ifndef SEQUENCE
#define SEQUENCE

#include "datasequence.h"
#include <string>
#include <vector>

namespace DeGenPrime
{
	class Sequence
	{
	public:
		Sequence();
		Sequence(std::string name);

		void SetName(std::string name);
		void SetList(std::vector<char> list);
		void CalculateScore(DataSequence data);
		void RemoveDashes();

		void Erase(int index);
		void PushBack(char c);
		void PushBack(std::string str);
		void PopBack();

		void Invert();
		void Reverse();
		
		std::string GetName() const;
		std::string Fasta() const;
		std::vector<char> GetCodes() const;

		int Score() const;
		int size() const;

		bool operator <(const Sequence& rhs) const;
	private:
		std::vector<char> _codes;
		std::string _name;
		int _score;
	};

} // End of DeGenPrime
#endif // SEQUENCE