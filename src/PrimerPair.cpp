// PrimerPair.cpp
#include <iostream>
#include <string>
#include <vector>
#include "PrimerPair.h"
#include "DataSequence.h"
#include "GlobalSettings.h"

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
		DataSequence fwdSub = fwd_data.SubSeq(_fwd.Index(), _fwd.Length());
		DataSequence revSub = rev_data.SubSeq(_rev.Index(), _rev.Length());
		int fwdSDiff = fwdSub.size() - fwdSub.ActualSize();
		int revSDiff = revSub.size() - revSub.ActualSize();
		ret += "Forward Primer:  Fwd Index[" + to_string(_fwd.Index());
		ret += "] Length: [" + to_string(fwdSub.ActualSize()) + "]";
		ret += (fwdSDiff > 0) ? (" (" + to_string(fwdSDiff) + " deletions)\n") :
			("\n");
		ret += "Codes: [" + fwdSub.Codes();
		ret += "]\nTm(NN): [" + to_string(fwdSub.NNMeltingTemperature()) + "] ";
		ret += "Tm(Basic): [" + to_string(fwdSub.BasicTemperature()) + "]\n";
		ret += "Penalty: [" + to_string(fwdSub.Penalty()) + "]\n";
		ret += "Reverse Primer: Rev Index[" + to_string(rev_data.RevIndex(_rev.Index()));
		ret += "] Length: [" + to_string(revSub.ActualSize()) + "]";
		ret += (revSDiff > 0) ? (" (" + to_string(revSDiff) + " deletions)\n") :
			("\n");
		ret += "Codes: [" + revSub.Codes();
		ret += "]\nTm(NN): [" + to_string(revSub.NNMeltingTemperature()) + "] ";
		ret += "Tm(Basic): [" + to_string(revSub.BasicTemperature()) + "]\n";
		ret += "Penalty: [" + to_string(revSub.Penalty()) + "]\n";
		ret += "Temperature Difference: [" + to_string(TempDiff());
		ret += "] Amplicon Length: [" + to_string(AmpSize()) + "]\n\n";
		return ret;
	}

	std::string PrimerPair::PrintShort(DataSequence fwd_data, DataSequence rev_data)
	{
		string ret = "";
		DataSequence fwdSub = fwd_data.SubSeq(_fwd.Index(), _fwd.Length());
		DataSequence revSub = rev_data.SubSeq(_rev.Index(), _rev.Length());
		ret += "fwd:\"" + fwdSub.Codes();
		ret += "\" rev:\"" + revSub.Codes();
		ret += "\"\n";
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
		/*if(GlobalSettings::GetSortByTemp())
		{
			return (_temp < rhs.TempDiff());
		}
		else
		{*/
			float q_l1 = _fwd.Penalty();
			float q_l2 = _rev.Penalty();
			float q_r1 = rhs.GetForward().Penalty();
			float q_r2 = rhs.GetReverse().Penalty();
			return((q_l1 * q_l2) < (q_r1 * q_r2));
		//}
	}
	
} // End of DeGenPrime