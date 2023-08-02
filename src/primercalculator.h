// ************************* PrimerCalculator *************************	//
// Purpose: Create an object to build lists of candidate primers and	//
//		evaluate those primers based on a list of critera.				//
// Constructors: Default: Create Empty PrimerCalculator.				//
// Mutators: InitializePrimers(data): Used to build list of primers 	//
//			based on a datasequence.									//
//		 InitializeBoundedPrimers(data): Used to build list of 			//
//			primers based on a datasequence.							//
//		 SetPrimers(): Standard Mutuator.								//
// Methods: Erase(): Erase a primer from the list.						//
//		PushBack(): Add a primer to the end of the list.				//
//		PrintSize(): Print the number of primers in the list.			//
//		PrintAll(): Print all primers in the list.						//
//		TooManyRepeats(): Private function used to check for too		//
//			much repetition within the primer.							//
// Filters:	FilterAll(Datasequence data, SequenceList list): Runs all	//
//			filters on the dataseqence/list and returns a string		//
//			with the filter information.								//
//		FilterDegeneracy(DataSequence data): Filter primers with		//
//			too much degeneracy.										//
//		FilterDeletions(DataSequence data, SequenceList list):			//
//			Filter primers with too many deletions.						//
//		FilterGCContent(DataSequence data): Filter primers with			//
//			too much or not enough GC Content.							//
//		FilterRepeats(DataSequence data): Filter primers with too		//
//			much repetition within the sequence.						//
//		FilterComplementaryEnds(DataSequence data): Filter primers		//
//			that have complementary ends.								//
//		FilterHairpins(DataSequence data): Filter primers likely		//
//			to form hairpins.											//
//		FilterDimers(DataSequence data): Filter primers likely to		//
//			form self or cross dimers.									//
//		FilterTemperature(DataSequence data, float offset): Filter		//
//			primers whose melting temperature is outside range.			//
// Accessors:	GetPrimers(): Returns the primer list.					//
//			size(): Returns the size of the primer list.				//
// ******************************************************************** //

#ifndef PRIMER_CALCULATOR
#define PRIMER_CALCULATOR

#include <vector>
#include "primer.h"
#include "sequencelist.h"
#include "datasequence.h"

namespace DeGenPrime
{
	class PrimerCalculator
	{
	public:
		PrimerCalculator();
		
		void InitializeTestPrimer(DataSequence data);
		void InitializePrimers(DataSequence data);
		void InitializeBoundedPrimers(DataSequence data, int lowerBound);
		void InitializeFromRegion(std::vector<Primer> region, DataSequence data);

		std::string FilterAll(DataSequence data);
		std::string FilterDegeneracy(DataSequence data);
		std::string FilterDeletions(DataSequence data);
		std::string FilterGCContent(DataSequence data);
		std::string FilterRepeats(DataSequence data);
		std::string FilterComplementaryEnds(DataSequence data);
		std::string FilterHairpins(DataSequence data);
		std::string FilterDimers(DataSequence data);
		std::string FilterTemperature(DataSequence data, float offset);

		void Erase(int index);
		void PushBack(Primer primer);

		void PrintSize();
		std::string PrintAll();
		void Sort();

		void SetPrimers(std::vector<Primer> primerList);
		
		std::vector<Primer> GetPrimers() const;
		int IndexOf(DataSequence data, std::string str) const;
		int size() const;
	private:
		int _OriginalSize;
		int TooManyRepeats(int size);
		std::vector<Primer> _primers;
		std::string FilterMessage(std::string func, int filtercount);
	};

} // End of DeGenPrime
#endif // PRIMER_CALCULATOR