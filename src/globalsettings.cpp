// GlobalSettings.cpp
#include <string>
#include "globalsettings.h"
#include "global.h"

using namespace std;

namespace DeGenPrime
{
	GlobalSettings::GlobalSettings() { }

	void GlobalSettings::SetMinimumAmplicon(int amplicon) { _ampLength = (0 < amplicon) ? amplicon : 0; }
	void GlobalSettings::SetBeginningNucleotide(int begin){ _beginningNucleotide = (0 < begin) ? begin : 0; }
	void GlobalSettings::SetDeltaG(float g) {_deltag = g;}
	void GlobalSettings::SetEndingNucleotide(int ending) {_endingNucleotide = ending; }
	void GlobalSettings::SetNonDegenerate(bool non) {_nonDegenerate = non; }
	void GlobalSettings::SetMeasureByAmpliconSize(bool size) {_measureByAmpliconSize = size; }
	void GlobalSettings::SetProteinSequence(bool seq) {_proteinSequence = seq; }
	void GlobalSettings::SetMinimumTemperature(float temp) {_minTemp = (MIN_PRIMER_TEMP < temp) ? temp : MIN_PRIMER_TEMP;}
	void GlobalSettings::SetMaximumTemperature(float temp) {_maxTemp = (MAX_PRIMER_TEMP > temp) ? temp : MAX_PRIMER_TEMP;}
	void GlobalSettings::SetMaximumPrimerLength(int len) {_maxLen = (MAX_PRIMER_LENGTH > len) ? len : MAX_PRIMER_LENGTH;}
	void GlobalSettings::SetMinimumPrimerLength(int len) {_minLen = (MIN_PRIMER_LENGTH < len) ? len : MIN_PRIMER_LENGTH;}
	void GlobalSettings::SetPrimerConcentration(float primer_conc) {_primerConcentration = (MIN_PRIMER_CONC < primer_conc) ? primer_conc : MIN_PRIMER_CONC;}
	void GlobalSettings::SetMonoIonConcentration(float salt_conc) {_monovalentIonConcentration = (MIN_SALT_CONC < salt_conc) ? salt_conc : MIN_SALT_CONC; }
	void GlobalSettings::SetMaximumReturnPrimers(int max) {_maxPrimers = (MAX_PRIMER_RETURNS > max) ? max : MAX_PRIMER_RETURNS;}
	void GlobalSettings::SetThermodynamicTemperature(float temp) {_thermodynamicTemperature = temp;}
	void GlobalSettings::SetBeginFlag(bool begin) {_beginflag = begin;}
	void GlobalSettings::SetEndFlag(bool ending) {_endflag = ending;}
	void GlobalSettings::SetRunTest(bool test) {_testRun = test;}
	void GlobalSettings::SetRunInvRev(bool test) {_invRevRun = test;}
	void GlobalSettings::SetSearchFwd(bool search) {_SearchFwd = search;}
	void GlobalSettings::SetSearchRev(bool search) {_SearchRev = search;}
	void GlobalSettings::SetTestValue(std::string str) {_testStr = str;}
	void GlobalSettings::SetInvRevValue(std::string str) {_invRevValue = str;}
	void GlobalSettings::SetSearchFwdArg(std::string str) {_searchFwdArg = str;}
	void GlobalSettings::SetSearchRevArg(std::string str) {_searchRevArg = str;}
	void GlobalSettings::SetSortByTemp(bool temp) {_sortbytemp = temp;}
	void GlobalSettings::SetUserTemp(bool temp) {_userTemp = temp;}

	int GlobalSettings::GetMinimumAmplicon() { return _ampLength; }
	int GlobalSettings::GetBeginningNucleotide() { return _beginningNucleotide; }
	float GlobalSettings::GetDeltaG() { return _deltag; }
	int GlobalSettings::GetEndingNucleotide() { return _endingNucleotide; }
	bool GlobalSettings::GetNonDegenerate() { return _nonDegenerate; }
	bool GlobalSettings::GetMeasureByAmpliconSize() { return _measureByAmpliconSize; }
	bool GlobalSettings::GetProteinSequence() { return _proteinSequence; }
	bool GlobalSettings::GetBeginFlag() { return _beginflag; }
	bool GlobalSettings::GetEndFlag() { return _endflag; }
	float GlobalSettings::GetMinimumTemperature() { return _minTemp; }
	float GlobalSettings::GetMaximumTemperature() { return _maxTemp; }
	int GlobalSettings::GetMaximumPrimerLength() { return _maxLen; }
	int GlobalSettings::GetMinimumPrimerLength() { return _minLen; }
	float GlobalSettings::GetPrimerConcentration() { return _primerConcentration; }
	float GlobalSettings::GetMonoIonConcentration() { return _monovalentIonConcentration; }
	int GlobalSettings::GetMaximumReturnPrimers() { return _maxPrimers; }
	float GlobalSettings::GetThermodynamicTemperature() { return _thermodynamicTemperature; }
	bool GlobalSettings::GetRunTest() { return _testRun; }
	bool GlobalSettings::GetRunInvRev() { return _invRevRun; }
	bool GlobalSettings::GetSearchFwd() { return _SearchFwd; }
	bool GlobalSettings::GetSearchRev() { return _SearchRev; }
	bool GlobalSettings::GetSortByTemp() { return _sortbytemp; }
	bool GlobalSettings::GetUserTemp() { return _userTemp; }
	string GlobalSettings::GetTestValue() { return _testStr; }
	string GlobalSettings::GetInvRevValue() { return _invRevValue; }
	string GlobalSettings::GetSearchFwdArg() { return _searchFwdArg; }
	string GlobalSettings::GetSearchRevArg() { return _searchRevArg; }
} // End of DeGenPrime