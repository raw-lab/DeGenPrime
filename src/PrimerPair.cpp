// PrimerPair.cpp
#include <iostream>
#include <string>
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

	std::string PrimerPair::Print(DataSequence fwd_data, DataSequence rev_data)
	{
		string ret = "";
		ret += "Forward Primer: Index[" + to_string(_fwd.Index());
		ret += "] Reverse Index[" + to_string(fwd_data.RevIndex(_fwd.Index()));
		ret += "] Length: [" + to_string(_fwd.Length()) + "] ";
		DataSequence fwdSub = fwd_data.SubSeq(_fwd.Index(), _fwd.Length());
		ret += "Codes: [" + fwdSub.Codes();
		ret += "]\nTm(NN): [" + to_string(fwdSub.NNMeltingTemperature()) + "] ";
		ret += "Tm(Basic): [" + to_string(fwdSub.BasicTemperature()) + "]\n";
		ret += "Reverse Primer: Index[" + to_string(_rev.Index());
		ret += "] Reverse Index[" + to_string(rev_data.RevIndex(_rev.Index()));
		ret += "] Length: [" + to_string(_rev.Length()) + "] ";
		DataSequence revSub = rev_data.SubSeq(_rev.Index(), _rev.Length());
		ret += "Codes: [" + revSub.Codes();
		ret += "]\nTm(NN): [" + to_string(revSub.NNMeltingTemperature()) + "] ";
		ret += "Tm(Basic): [" + to_string(revSub.BasicTemperature()) + "]\n";
		ret += "Temperature Difference: [" + to_string(TempDiff());
		ret += "] Amplicon Length: [" + to_string(AmpSize()) + "] ";
		DataSequence mostStable = (fwdSub.Gibbs() > revSub.Gibbs()) ? fwdSub : revSub;
		DataSequence product = fwdSub.Gibbs() > revSub.Gibbs() ? fwd_data.SubSeq(_fwd.Index(), AmpSize())
			: rev_data.SubSeq(_rev.Index(), AmpSize());
		ret += "Annealing Temp: [" + to_string(mostStable.BasicAnneal(product)) + "]\n\n";
		return ret;
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
		float firstTemp  = fwd_data.SubSeq(_fwd.Index(), _fwd.Length()).NNMeltingTemperature();
		float secondTemp = rev_data.SubSeq(_rev.Index(), _rev.Length()).NNMeltingTemperature();
		float difference = firstTemp > secondTemp ? (firstTemp - secondTemp) : secondTemp - firstTemp;
		return difference;
	}

	
	bool PrimerPair::operator <(const PrimerPair& rhs) const
	{
		return (_temp < rhs.TempDiff());
	}
	
} // End of DeGenPrime