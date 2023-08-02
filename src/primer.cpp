// Primer.cpp
#include <iostream>
#include <string>
#include "datasequence.h"
#include "primer.h"

using namespace std;

namespace DeGenPrime
{
	Primer::Primer() { }
	Primer::Primer(int index, int length)
	{
		_Length = length;
		_Index = index;
	}

	std::string Primer::Print()
	{
		string ret = "[index=" + to_string(_Index);
		ret += "][length=" + to_string(_Length);
		ret += "]\n";
		return ret;
	}

	void Primer::SetPenalty(float penalty) 
	{
		_Penalty = penalty;
	}

	bool Primer::operator <(const Primer& rhs) const
	{
		return (_Penalty < rhs.Penalty());
	}

	int Primer::Length() const { return _Length; }
	int Primer::Index() const { return _Index; }
	float Primer::Penalty() const { return _Penalty; }
} // End of DeGenPrime