#include <iostream>
#include <string>
#include <vector>
#include "datasequence.h"
#include "sequence.h"

using namespace std;

namespace DeGenPrime
{
	Sequence::Sequence() { }
	Sequence::Sequence(string name) { _name = name; }

	void Sequence::SetName(string name) { _name = name; }
	void Sequence::SetList(std::vector<char> list) 
	{	
		_codes.clear();
		for(char c : list)
		{
			_codes.push_back(c);
		}
	}
	void Sequence::CalculateScore(DataSequence data)
	{
		int count = 0; // number of mismatches
		string MC = data.MC();
		for(int i = 0;i < data.size();i++)
		{
			if(MC[i] != _codes[i])count++;
		}
		_score = count;
	}

	void Sequence::Erase(int index)
	{
		_codes.erase(_codes.begin() + index);
	}
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

	void Sequence::RemoveDashes()
	{
		std::vector<char> chars;
		for(int i = 0;i < _codes.size();i++)
		{
			char c = _codes[i];
			if(c != '-')
			{
				chars.push_back(c);
			}
		}
		_codes.clear();
		for(char c : chars)
		{
			_codes.push_back(c);
		}
	}

	std::string Sequence::GetName() const { return _name; }
	std::string Sequence::Fasta() const
	{
		string fasta = ">" + _name + "\n";
		string codes(_codes.begin(), _codes.end());
		fasta += codes + "\n";
		return fasta;
	}
	std::vector<char> Sequence::GetCodes() const { return _codes; }

	int Sequence::Score() const { return _score; }
	int Sequence::size() const { return _codes.size(); }

	bool Sequence::operator <(const Sequence& rhs) const
	{
		return (Score() < rhs.Score());
	}
} // End of DeGenPrime