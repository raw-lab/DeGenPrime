// ***************************** PrimerID ***************************** //
// Purpose:	Provide a data structure for information about a primer.	//
// Members:	Length:	How many nucleotides are in the primer.				//
//		Index:	Where the primer begins.								//
//		Penalty: The penalty of the primer.								//
// Functions:	Print():	Print data about the primer.				//
// Operator: <: Used to sort primers by lowest penalty first.			//
// ********************************************************************	//

#ifndef PRIMER_CLASS
#define PRIMER_CLASS

#include <string>
#include "datasequence.h"

namespace DeGenPrime
{
	class Primer
	{
	public:
		Primer();
		Primer(int index, int length);

		void SetPenalty(float penalty);

		int Length() const;
		int Index() const;
		float Penalty() const;

		std::string Print();

		bool operator <(const Primer& rhs) const;

	private:
		int _Length;
		int _Index;
		float _Penalty;
	};
} // End of DeGenPrime
#endif // PRIMER_CLASS