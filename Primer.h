// ***************************** PrimerID ***************************** //
// Purpose:	Provide a data structure for information about a primer.	//
// Members:	Length:	How many nucleotides are in the primer.		//
//		Index:	Where the primer begins.				//
//		isForward:	Is this a forward or reverse primer. (bool)	//
// ********************************************************************	//

#ifndef PRIMER_CLASS
#define PRIMER_CLASS

namespace DeGenPrime
{
	class Primer
	{
	public:
		Primer(int index, int length);

		void Print();

		int Length() const;
		int Index() const;

	private:
		int _Length;
		int _Index;
	};
} // End of DeGenPrime
#endif // PRIMER_CLASS