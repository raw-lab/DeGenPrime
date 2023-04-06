// ************************* PrimerCalculator *************************	//
// Purpose: Create an object to build lists of candidate primers and	//
//		evaluate those primers based on a list of critera.		//
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