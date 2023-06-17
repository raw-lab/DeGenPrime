// Primer.cpp
#include <iostream>
#include <string>
#include "DataSequence.h"
#include "Primer.h"

using namespace std;

namespace DeGenPrime
{
	Primer::Primer() { }
	Primer::Primer(int index, int length)
	{
		_Length = length;
		_Index = index;
	}
	/*
	Primer::Primer(int index, int length, DataSequence src_data)
	{
		_Length = length;
		_Index = index;
		DataSequence sub = src_data.SubSeq(index, length);
		_Quality = sub.Quality();
	}*/

	void Primer::Print()
	{
		cout << "[index=" << _Index;
		cout << "][length=" << _Length;
		cout << "]\t";
	}

	void Primer::SetQuality(float quality) {_Quality = quality;}

	bool Primer::operator <(const Primer& rhs) const
	{
		return (_Quality < rhs.Quality());
	}

	int Primer::Length() const { return _Length; }
	int Primer::Index() const { return _Index; }
	float Primer::Quality() const { return _Quality; }
} // End of DeGenPrime