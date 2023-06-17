// ***************************** PrimerID ***************************** //
// Purpose:	Provide a data structure for information about a primer.	//
// Members:	Length:	How many nucleotides are in the primer.		//
//		Index:	Where the primer begins.				//
// Functions:	Print():	Print data about the primer.			//
// ********************************************************************	//

#ifndef PRIMER_CLASS
#define PRIMER_CLASS

#include "DataSequence.h"

namespace DeGenPrime
{
	class Primer
	{
	public:
		Primer();
		Primer(int index, int length);
		//Primer(int index, int length, DataSequence src_data);

		void Print();
		void SetQuality(float quality);

		int Length() const;
		int Index() const;
		float Quality() const;

		bool operator <(const Primer& rhs) const;

	private:
		int _Length;
		int _Index;
		float _Quality;
	};
} // End of DeGenPrime
#endif // PRIMER_CLASS