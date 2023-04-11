// Primer.cpp
#include <iostream>
#include <string>
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

	void Primer::Print()
	{
		cout << "[index=" << _Index;
		cout << "][length=" << _Length;
		cout << "]\t";
	}

	int Primer::Length() const { return _Length; }
	int Primer::Index() const { return _Index; }
} // End of DeGenPrime