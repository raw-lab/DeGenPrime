// ******************************* PrimerPairList.h ******************************* //
// Purpose: To have a list of all pairs of forward and reverse primers to run		//
//		sorting and filtering operations.											//
// Constructors: Default: Create an empty PrimerPairList.							//
//			PrimerPairList(fwd_seq, rev_seq, pair_list): Build list from			//
//				a given list of primer pairs, forward and reverse data.				//
// Mutators: CreateList(fwd_seq, rev_seq, fwd_list, rev_list) creates list of		//
//				primer pairs from forward and reverse primer lists.					//
// Methods: Append(list): Adds the argument to the end of the list of primerpairs.	//
//		Erase():	Erases a primerpair from the list.								//
//		SubList(start, length): Returns a sublist starting at start and going		//
//			length distance.														//
//		PushBack(): Adds a primerpair to the end of the list.						//
//		Sort():	Sort the list from lowest to highest temperature difference			//
//			between the forward and reverse primers.								//
// Filters:	FilterAmpliconLength(): Filters all primers whose amplicon length is	//
//			less than the specified minimum.										//
//		FilterTemperatureDifference(): Filters all primers whose temperature		//
//			difference is higher than the specified minimum.						//
//		FilterAnnealingTemperature(): Filters all primers whose annealing			//
//			temperature is 10 degrees celsius more or less than the melting			//
//			temperature of the most stable primer.  WARNING: This function is		//
//			slow and should only be called on a small list of primerpairs.			//
// Functions: PrintAll(): Prints the list of primer pairs.							//
//		  PrintSize(): Prints the size of the list of primer pairs.					//
// Accessors: GetPairs(): Returns the list of primer pairs.							//
//		  size(): Returns the size of the list of primer pairs.						//
// ******************************************************************************** //

#ifndef PRIMER_PAIR_LIST
#define PRIMER_PAIR_LIST

#include <string>
#include <vector>
#include "primerpair.h"
#include "primer.h"
#include "datasequence.h"

namespace DeGenPrime
{
	class PrimerPairList
	{
	public:
		PrimerPairList();
		PrimerPairList(DataSequence fwd_seq, DataSequence rev_seq, std::vector<PrimerPair> pair_list);
		PrimerPairList SubList(int startIndex, int length);

		void Append(PrimerPairList list);
		
		void CreateFromRange(DataSequence fwd_seq, DataSequence rev_seq,
			std::vector<Primer> fwd_list, std::vector<Primer> rev_list, 
			int fwd_begin, int fwd_end, int rev_begin, int rev_end);
		std::string CreateList(DataSequence fwd_seq, DataSequence rev_seq, std::vector<Primer> fwd_list, std::vector<Primer> rev_list);
		void Erase(int index);
		void PushBack(PrimerPair pair);

		std::string FilterMessage(std::string func, int filtercount);
		std::string FilterAmpliconLength();
		std::string FilterTemperatureDifference();
		int FilterAnnealingTemp(DataSequence fwd, DataSequence rev, int ignore);
		int FilterUnique();
		int PartitionCount() const;
		int PartitionCount(int fwd_size, int rev_size) const;

		void Sort();

		std::string PrintAll(DataSequence fwd, DataSequence rev);
		std::string PrintAllShort(DataSequence fwd, DataSequence rev);
		std::string CreateCSV(DataSequence fwd, DataSequence rev);

		std::vector<PrimerPair> GetPairs() const;

		int size() const;
	private:
		bool comparator(const PrimerPair& lhs, const PrimerPair& rhs);
		int _OriginalSize;
		std::vector<PrimerPair> _pairs;
	};
} // End of DeGenPrime
#endif // PRIMER_PAIR_LIST