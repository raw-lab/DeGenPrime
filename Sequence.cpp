#include <iostream>
#include <string>
#include <vector>
#include "Sequence.h"

using namespace std;

namespace DeGenPrime
{
	Sequence::Sequence() { }
	Sequence::Sequence(string name) { _name = name; }

	void Sequence::SetName(string name) { _name = name; }
	void Sequence::SetList(std::vector<char> list) {_codes = list; }

	void Sequence::PushBack(char c) { _codes.push_back(c); }
	void Sequence::PopBack() { _codes.pop_back(); }

	string Sequence::GetName() const { return _name; }
	std::vector<char> Sequence::GetCodes() const { return _codes; }

	int Sequence::size() const { return _codes.size(); }
} // End of DeGenPrime