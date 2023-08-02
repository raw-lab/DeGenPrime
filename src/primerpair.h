// *********************************** PrimerPair.h ***********************************	//
// Purpose:	Group forward and reverse primers together so operations involving filters	//
//		and output can easily know which primers to relate.								//
// Constructors: Default: empty constructor												//
//	PrimerPair(forward, reverse, fwd_data,rev_data): Stores forward and reverse as		//
//			the primers and uses fwd_data and rev_data to get some of the				//
//			other private member values saved by the primerpair.						//
// Accessors:	GetForward():	Returns forward primer									//
//			GetReverse():	Returns reverse primer										//
//			AmpSize():		Return the size of the amplicon between primers.			//
//			TempDiff():		Returns the temperature difference between primers.			//
// Functions:	Print(fwd_data, rev_data): Prints Datasequence data from primers		//
//				using fwd_data and rev_data as the template for the data.				//
// Operator:	< :	Ranks primerpairs by least temperature difference.					//
// ************************************************************************************	//
#ifndef PRIMER_PAIR
#define PRIMER_PAIR

#include "datasequence.h"
#include "primer.h"

namespace DeGenPrime
{
	class PrimerPair
	{
	public:
		PrimerPair();
		PrimerPair(Primer forward, Primer reverse, DataSequence fwd_data, DataSequence rev_data);

		Primer GetForward() const;
		Primer GetReverse() const;

		std::string Print(DataSequence fwd_data, DataSequence rev_data);
		std::string PrintShort(DataSequence fwd_data, DataSequence rev_data);

		int AmpSize() const;
		float TempDiff() const;

		bool operator <(const PrimerPair& rhs) const;
	private:
		Primer _fwd;
		Primer _rev;

		int _length;
		float _temp;
		float _anneal;

		int AmpliconLength(DataSequence data);
		float TemperatureDifference(DataSequence fwd_data, DataSequence rev_data);
	};
} // End of DeGenPrime
#endif // PRIMER_PAIR