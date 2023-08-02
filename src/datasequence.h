// ****************************** class DataSequence ******************************
// Purpose:	class DataSequence is a list class containing only a list of 
//		DataNodes and designed to perform specific opeations on the list.
// Constructor:	- Default constructor to initialize a blank object.  Another
//		constructor to generate DataSequence directly from a string.
// Mutators: void SetList(vector<DataNode> catalog) standard mutator.
// Functions: PushBack(DataNode node): Appends node to list.
//		PopBack(): Removes last node of list.
//		Erase(int i): Removes the node at the argument index.
//		string Print(): Returns a string of the primer information.
//	    string Codes(): Returns the list of the Codes in the sequence to cout.
//		string MC(): Reutrns the list of the Most common nucleotide chars.
// 		SubSeq(start,length): Returns a length list of DataNodes beginning
//				 of the start index.  Return is designed to be used
//				 with SetList to make a newDataSequence.
// 		InvSeq(): Returns a new DataSequence where all of the codes
//				have been changed to the inverse of their code.
//		RevSeq(): Returns a new DataSequence in the reverse order
//				of the calling DataSequence.
//		Enthalpy(): Returns the Enthalpy of the DataSequence (float).
//			Values from https://www.pnas.org/doi/10.1073/pnas.95.4.1460
//		Entropy(): Returns the Entropy of the DataSequence (float).
//			Values from https://www.pnas.org/doi/10.1073/pnas.95.4.1460
//		GCRatio(): Returns the float ratio of G and C nucleotides in the sequence.
//	    Gibbs(): Returns the Gibbs of the DataSequence at temp (float).
//			Values from https://www.pnas.org/doi/10.1073/pnas.95.4.1460
//  BasicTemperature(): Returns the basic melting temperature of the DataSequence
//				(float) using the formula:
//					64.9 + (41.0 * (gc - 16.4) / (float)(gc + at) )
//				where gc is the number of G or C nucleotides
//				and at is the number of A or T nucleotides
//				http://biotools.nubic.northwestern.edu/OligoCalc.html
//   RlnK(primer_conc): Takes the natural log of the argument, multiplies it by
//				the ideal gas constant in cal/(mol * K) (1.9872) and
//				returns this value (float).
// MonoIonMod(saltcon): Returns the log of the argument times 16.6 (float).  This
//				value modifyings entropy for monovalent cations.
// NNMeltingTemperature(salt_conc, primer_conc): This function uses a nearest
//				neighbor formula to compute a more accurate melting temp
//				of the primer.  The formula can be found at:
// 	https://academic.oup.com/bioinformatics/article/21/6/711/199347?login=true
//				there is one adjustment where salt concentration modifies
//				the entropy in the denominator of this equation using
//				0.368 * N * ln(salt_conc*10^-3) where N is the nucleotide
//				length.  https://www.pnas.org/doi/10.1073/pnas.95.4.1460
//	ProductMelt(): Returns the melting point of the product using Santalucia formula.
//				https://www.pnas.org/doi/10.1073/pnas.95.4.1460. 
//	BasicAnneal(product, salt_conc, primer_conc): finds the Annealing temp
//				of the primer.  It uses the formula:
//				Ta = 0.3 * Tm(primer) + 0.7 * Tm(product) - 14.9.
//			http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html
//	Penalty(): Returns a floating point number showing the penalty of the primer.
//		Primers get more penalty by being outside temperature ranges, having
//		degenerate nucleotides, having too much repetition, and compelementary ends.
// Accessors: GetDataSequence(): Returns _list.
//			size(): Returns the size of the list.
//			ActualSize(): Returns the size of the list without deletions.
//   	CountMatches(DataSequence): Returns the number of nucleotides that match the
//			arguments's most common at each index.
//		RevIndex(index): Convert a forward index to a reverse index.
//	    IndexOf(DataSequence): Searches caller for a subsequence that matches
//					   the argument and returns a value of that index.
//					   Returns -1 if no match was found.
//	 	checkMatch(DataSequence): Returns true if argument has the same most common
//					   nucleotides and size as the calling Datasequence.
//		isEmpty(): Returns true if the primer is entirely deletions or size zero.
// Private members: vector<DataNode> list: the list of data
// ********************************************************************************

#ifndef DATA_SEQUENCE
#define DATA_SEQUENCE

#include <fstream>
#include <vector>
#include "datanode.h"
#include "global.h"

namespace DeGenPrime
{
	class DataSequence
	{
	public:
		DataSequence();
		DataSequence(std::string str);
		
		void SetList(std::vector<DataNode> catalog);
		void PushBack(DataNode node);
		void PopBack();
		void Erase(int index);

		std::string Print();
		std::string Codes() const;
		std::string Consensus(std::vector<int> fwd_ind, std::vector<int> fwd_len, bool cons);
		std::string MC();

		DataSequence SubSeq(int startIndex,int length);
		DataSequence InvSeq();
		DataSequence RevSeq();

		float Enthalpy() const;
		float Entropy() const;
		float Gibbs() const;
		float GCRatio() const;
		float BasicTemperature() const;
		float RlnK() const;
		float MonoIonMod() const;
		float ProductMelt() const;
		float NNMeltingTemperature() const;
		float BasicAnneal(DataSequence product);
		float Penalty() const;
		
		std::vector<DataNode> GetDataSequence() const;

		int CountMatches(DataSequence data) const;
		int RevIndex(int index) const;
		int IndexOf(DataSequence data) const;
		int size() const;
		int ActualSize() const;

		bool checkMatch(DataSequence data) const;
		bool isEmpty() const;
	private:
		std::vector<DataNode> _list;
	};
} // End of DeGenPrime
#endif // DATA_SEQUENCE