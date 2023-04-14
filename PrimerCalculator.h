// ************************* PrimerCalculator *************************	//
// Purpose: Create an object to build lists of candidate primers and	//
//		evaluate those primers based on a list of critera.		//
// Constructors: Default: Create Empty PrimerCalculator.			//
// Mutators: InitializePrimers(data, ampliconlength): Used to build	//
//			list of primers based on a datasequence and a given	//
//			amplicon length.  						//
//		 InitializeBoundedPrimers(data, lowerbound): Used to build	//
//			list of primers based on a datasequence and a max	//
//			bound for the end of the primer.				//
//		 SetPrimers(): Standard Mutuator.					//
// Methods: Erase(): Erase a primer from the list.				//
//		PushBack(): Add a primer to the end of the list.		//
//		PrintSize(): Print the number of primers in the list.		//
//		PrintAll(): Print all primers in the list.			//
//		TooManyRepeats(): Private function used to check for too	//
//			much repetition within the primer.				//
// Filters:	FilterDegeneracy(DataSequence data): Filter primers with	//
//			too much degeneracy.						//
//		FilterDeletions(DataSequence data, SequenceList list):	//
//			Filter primers with too many deletions.			//
//		FilterGCContent(DataSequence data): Filter primers with	//
//			too much or not enough GC Content.				//
//		FilterRepeats(DataSequence data): Filter primers with too	//
//			much repetition within the sequence.			//
//		FilterComplementaryEnds(DataSequence data): Filter primers	//
//			that have complementary ends.					//
//		FilterHairpins(DataSequence data): Filter primers likely	//
//			to form hairpins.							//
//		FilterDimers(DataSequence data): Filter primers likely to	//
//			form self or cross dimers.					//
//		FilterTemperature(DataSequence data, float offset): Filter	//
//			primers whose melting temperature is outside range.	//
// Accessors:	GetPrimers(): Returns the primer list.			//
//			size(): Returns the size of the primer list.		//
// ******************************************************************** //

#ifndef PRIMER_CALCULATOR
#define PRIMER_CALCULATOR

#include <vector>
#include "Primer.h"
#include "SequenceList.h"
#include "DataSequence.h"

namespace DeGenPrime
{
	class PrimerCalculator
	{
	public:
		PrimerCalculator();
		
		void InitializePrimers(DataSequence data, int AmpliconLength);
		void InitializeBoundedPrimers(DataSequence data, int lowerBound);

		void FilterDegeneracy(DataSequence data);
		void FilterDeletions(DataSequence data, SequenceList list);
		void FilterGCContent(DataSequence data);
		void FilterRepeats(DataSequence data);
		void FilterComplementaryEnds(DataSequence data);
		void FilterHairpins(DataSequence data);
		void FilterDimers(DataSequence data);
		void FilterTemperature(DataSequence data, float offset);

		void Erase(int index);
		void PushBack(Primer primer);

		void PrintSize();
		void PrintAll(); // Test Function

		void SetPrimers(std::vector<Primer> primerList);
		
		std::vector<Primer> GetPrimers() const;
		int size() const;
	private:
		int TooManyRepeats(int size);
		std::vector<Primer> _primers;
	};

} // End of DeGenPrime
#endif // PRIMER_CALCULATOR