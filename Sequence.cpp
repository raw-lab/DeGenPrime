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
	void Sequence::PushBack(string str)
	{
		for(char c : str)
		{
			if(c != '\0')
			{
				_codes.push_back(c);
			}
		}
	}
	
	void Sequence::PopBack() { _codes.pop_back(); }

	void Sequence::Invert()
	{
		for(int i = 0; i < _codes.size();i++)
		{
			switch(_codes[i])
			{
				case 'A':
				case 'a':
					_codes[i] = 'T';
					break;
				case 'C':
				case 'c':
					_codes[i] = 'G';
					break;
				case 'G':
				case 'g':
					_codes[i] = 'C';
					break;
				case 'T':
				case 't':
					_codes[i] = 'A';
					break;
				default:
					break;
			}
		}
	}

	void Sequence::Reverse()
	{
		std::vector<char> rev_seq;
		for(int i = _codes.size() - 1;i >= 0;i--)
		{
			rev_seq.push_back(_codes[i]);
		}
		_codes = rev_seq;
	}

	string Sequence::GetName() const { return _name; }
	std::vector<char> Sequence::GetCodes() const { return _codes; }

	int Sequence::size() const { return _codes.size(); }
} // End of DeGenPrime