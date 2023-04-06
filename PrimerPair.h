// PrimerPair.h
#ifndef PRIMER_PAIR
#define PRIMER_PAIR

#include "DataSequence.h"
#include "Primer.h"

namespace DeGenPrime
{
	class PrimerPair
	{
	public:
		PrimerPair();
		PrimerPair(Primer forward, Primer reverse, DataSequence fwd_data, DataSequence rev_data);

		Primer GetForward() const;
		Primer GetReverse() const;

		int AmpSize() const;
		float TempDiff() const;

		bool operator <(const PrimerPair& rhs) const;
	private:
		Primer _fwd;
		Primer _rev;

		int _length;
		float _temp;

		int AmpliconLength(DataSequence data);
		float TemperatureDifference(DataSequence fwd_data, DataSequence rev_data);
	};
} // End of DeGenPrime
#endif // PRIMER_PAIR