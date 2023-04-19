// PrimerPair.cpp
#include <iostream>
#include <vector>
#include "PrimerPair.h"
#include "DataSequence.h"

using namespace std;

namespace DeGenPrime
{
	PrimerPair::PrimerPair() { }
	PrimerPair::PrimerPair(Primer forward, Primer reverse, DataSequence fwd_data, DataSequence rev_data)
	{
		_fwd = forward;
		_rev = reverse;
		_length = AmpliconLength(fwd_data);
		_temp = TemperatureDifference(fwd_data, rev_data);
	}

	Primer PrimerPair::GetForward() const { return _fwd; }
	Primer PrimerPair::GetReverse() const { return _rev; }

	int PrimerPair::AmpSize() const { return _length; }
	float PrimerPair::TempDiff() const { return _temp; }

	void PrimerPair::Print(DataSequence fwd_data, DataSequence rev_data)
	{
		cout << "Forward Primer: Index[" << _fwd.Index();
		cout << "] Length: [" << _fwd.Length() << "] ";
		DataSequence fwdSub = fwd_data.SubSeq(_fwd.Index(), _fwd.Length());
		cout << "Codes: [";
		for(DataNode node : fwdSub.GetDataSequence())
		{
			cout << node.GetCode();
		}
		cout << "] Tm: [" << fwdSub.NNMeltingTemperature(50.0, 50.0) << "]" << endl;
		cout << "Reverse Primer: Index[" << _rev.Index();
		cout << "] Length: [" << _rev.Length() << "] ";
		DataSequence revSub = rev_data.SubSeq(_rev.Index(), _rev.Length());
		cout << "Codes: [";
		for(DataNode node : revSub.GetDataSequence())
		{
			cout << node.GetCode();
		}
		cout << "] Tm: [" << revSub.NNMeltingTemperature(50.0, 50.0) << "]" << endl;
		cout << "Temperature Difference: [" << TempDiff();
		cout << "] Amplicon Length: [" << AmpSize() << "] ";
		DataSequence mostStable = (fwdSub.Gibbs(50.0) > revSub.Gibbs(50.0)) ? fwdSub : revSub;
		DataSequence product = fwdSub.Gibbs(50.0) > revSub.Gibbs(50.0) ? fwd_data.SubSeq(_fwd.Index(), AmpSize())
			: rev_data.SubSeq(_rev.Index(), AmpSize());
		cout << "Annealing Temp: [" << mostStable.BasicAnneal(product, 50.0, 50.0) << "]\n" << endl;
	}

	int PrimerPair::AmpliconLength(DataSequence data)
	{
		// The amplicon runs from the beginning of the forward index to the beginning of the reverse index.
		// Translate the reverse index to a forward index and subtract.
		// Either the forward or reverse DataSequence is sufficient for this argument.
		int begin = _fwd.Index();
		int end = data.RevIndex(_rev.Index());
		return (end - begin);
	}

	float PrimerPair::TemperatureDifference(DataSequence fwd_data, DataSequence rev_data)
	{
		float firstTemp  = fwd_data.SubSeq(_fwd.Index(), _fwd.Length()).NNMeltingTemperature(50.0, 50.0);
		float secondTemp = rev_data.SubSeq(_rev.Index(), _rev.Length()).NNMeltingTemperature(50.0, 50.0);
		float difference = firstTemp > secondTemp ? (firstTemp - secondTemp) : secondTemp - firstTemp;
		return difference;
	}

	
	bool PrimerPair::operator <(const PrimerPair& rhs) const
	{
		return (_temp < rhs.TempDiff());
	}
	
} // End of DeGenPrime