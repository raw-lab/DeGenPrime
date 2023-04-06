// PrimerPairList.h

#ifndef PRIMER_PAIR_LIST
#define PRIMER_PAIR_LIST

#include <vector>
#include "PrimerPair.h"
#include "Primer.h"
#include "DataSequence.h"

namespace DeGenPrime
{
	class PrimerPairList
	{
	public:
		PrimerPairList();
		PrimerPairList(DataSequence fwd_seq, DataSequence rev_seq, std::vector<PrimerPair> pair_list);
		
		void CreateList(DataSequence fwd_seq, DataSequence rev_seq, std::vector<Primer> fwd_list, std::vector<Primer> rev_list);
		void Erase(int index);

		void FilterAmpliconLength();
		void FilterTemperatureDifference();
		void FilterCrossDimers();

		void Sort();
		void PrintSize();

		std::vector<PrimerPair> GetPairs() const;

		int size() const;
	private:
		bool comparator(const PrimerPair& lhs, const PrimerPair& rhs);

		std::vector<PrimerPair> _pairs;
		DataSequence _fwd;
		DataSequence _rev;
	};
} // End of DeGenPrime
#endif // PRIMER_PAIR_LIST