// PrimerPair.cpp
#include <iostream>
#include <string>
#include <vector>
#include "primerpair.h"
#include "datasequence.h"
#include "globalsettings.h"
#include "format.h"

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
		string line = "";
		DataSequence fwdSub = fwd_data.SubSeq(_fwd.Index(), _fwd.Length());
		DataSequence revSub = rev_data.SubSeq(_rev.Index(), _rev.Length());
		
		line = "Amplicon Size: [" + to_string(AmpSize()) + "]";
		ret += Format(line, 25, Alignment::Center);
		line = "Temp Diff: [" + Format(TempDiff(), 2) + "]";
		ret += Format(line, STR_FORMAT - 50, Alignment::Center) + "\n";
		line = "Forward";
		ret += Format(line, STR_FORMAT / 2, Alignment::Center);
		line = "Reverse";
		ret += Format(line, STR_FORMAT / 2, Alignment::Center) + "\n";
		line = "[";
		line += fwdSub.Codes() + "]";
		ret += Format(line, STR_FORMAT / 2, Alignment::Center);
		line = "[";
		line += revSub.Codes() + "]";
		ret += Format(line, STR_FORMAT / 2, Alignment::Center) + "\n";
		line = "Basic Information";
		ret += Format(line, STR_FORMAT / 2, Alignment::Center);
		ret += Format(line, STR_FORMAT / 2, Alignment::Center) + "\n";
		line = "Range: ";
		ret += Format(line, STR_FORMAT / 4, Alignment::Right);
		line = "[";
		int digs1 = digits(_fwd.Index());
		int digs2 = digits(_fwd.Index() + _fwd.Length() - 1);
		int max_dig = (digs2 > digs1) ? digs2 : digs1;
		line += Format(_fwd.Index(), max_dig) + "-";
		line += Format(_fwd.Index() + _fwd.Length() - 1, max_dig) + "]";
		ret += Format(line, STR_FORMAT / 4, Alignment::Left);
		line = "Range: ";
		ret += Format(line, STR_FORMAT / 4, Alignment::Right);
		line = "[";
		digs2 = digits(rev_data.RevIndex(_rev.Index()));
		digs1 = digits(rev_data.RevIndex(_rev.Index() + _rev.Length()));
		max_dig = (digs2 > digs1) ? digs2 : digs1;
		line += Format(rev_data.RevIndex(_rev.Index()), max_dig) + "-";
		line += Format(rev_data.RevIndex(_rev.Index() + _rev.Length()), max_dig) + "]";
		ret += Format(line, STR_FORMAT / 4, Alignment::Left) + "\n";
		line = "Length: ";
		ret += Format(line, STR_FORMAT / 4, Alignment::Right);
		line = "[" + to_string(fwdSub.ActualSize()) + "]";
		ret += Format(line, STR_FORMAT / 4, Alignment::Left);
		line = "Length: ";
		ret += Format(line, STR_FORMAT / 4, Alignment::Right);
		line = "[" + to_string(revSub.ActualSize()) + "]";
		ret += Format(line, STR_FORMAT / 4, Alignment::Left) + "\n";
		line = "Penalty: ";
		ret += Format(line, STR_FORMAT / 4, Alignment::Right);
		line = "[" + Format(fwdSub.Penalty(), 2) + "]";
		ret += Format(line, STR_FORMAT / 4, Alignment::Left);
		line = "Penalty: "; 
		ret += Format(line, STR_FORMAT / 4, Alignment::Right);
		line = "[" + Format(revSub.Penalty(), 2) + "]";
		ret += Format(line, STR_FORMAT / 4, Alignment::Left) + "\n";
		line = "Melting Temperature";
		ret += Format(line, STR_FORMAT / 2, Alignment::Center);
		ret += Format(line, STR_FORMAT / 2, Alignment::Center) + "\n";
		line = "Basic: ";
		ret += Format(line, STR_FORMAT / 4, Alignment::Right);
		line = "[" + Format(fwdSub.BasicTemperature(), 2) + "]";
		ret += Format(line, STR_FORMAT / 4, Alignment::Left);
		line = "Basic: ";
		ret += Format(line, STR_FORMAT / 4, Alignment::Right);
		line = "[" + Format(revSub.BasicTemperature(), 2) + "]";
		ret += Format(line, STR_FORMAT / 4, Alignment::Left) + "\n";
		line = "NN: ";
		ret += Format(line, STR_FORMAT / 4, Alignment::Right);
		line = "[" + Format(fwdSub.NNMeltingTemperature(), 2) + "]";
		ret += Format(line, STR_FORMAT / 4, Alignment::Left);
		line = "NN: ";
		ret += Format(line, STR_FORMAT / 4, Alignment::Right);
		line = "[" + Format(revSub.NNMeltingTemperature(), 2) + "]";
		ret += Format(line, STR_FORMAT / 4, Alignment::Left) + "\n";
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
		float q_l1 = _fwd.Penalty();
		float q_l2 = _rev.Penalty();
		float q_r1 = rhs.GetForward().Penalty();
		float q_r2 = rhs.GetReverse().Penalty();
		return((q_l1 * q_l2) < (q_r1 * q_r2));
	}
	
} // End of DeGenPrime