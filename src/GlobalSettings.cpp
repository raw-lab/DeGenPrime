// GlobalSettings.cpp
#include "GlobalSettings.h"
#include "global.h"

using namespace std;

namespace DeGenPrime
{
	// GlobalSettings::GlobalSettings() { }

	void GlobalSettings::SetMinimumAmplicon(int amplicon) { _ampLength = (0 < amplicon) ? amplicon : 0; }
	void GlobalSettings::SetBeginningNucleotide(int begin){ _beginningNucleotide = (0 < begin) ? begin : 0; }
	void GlobalSettings::SetEndingNucleotide(int ending) {_endingNucleotide = ending; }
	void GlobalSettings::SetMeasureByAmpliconSize(bool size) {_measureByAmpliconSize = size; }
	void GlobalSettings::SetMinimumTemperature(float temp) {_minTemp = (MIN_PRIMER_TEMP < temp) ? temp : MIN_PRIMER_TEMP;}
	void GlobalSettings::SetMaximumTemperature(float temp) {_maxTemp = (MAX_PRIMER_TEMP > temp) ? temp : MAX_PRIMER_TEMP;}
	void GlobalSettings::SetPrimerConcentration(float primer_conc) {_primerConcentration = primer_conc;}
	void GlobalSettings::SetMonoIonConcentration(float salt_conc) {_monovalentIonConcentration = salt_conc; }
	void GlobalSettings::SetMaximumReturnPrimers(int max) {_maxPrimers = (MAX_PRIMER_RETURNS > max) ? max : MAX_PRIMER_RETURNS;}

	int GlobalSettings::GetMinimumAmplicon() { return _ampLength; }
	int GlobalSettings::GetBeginningNucleotide() { return _beginningNucleotide; }
	int GlobalSettings::GetEndingNucleotide() { return _endingNucleotide; }
	bool GlobalSettings::GetMeasureByAmpliconSize() { return _measureByAmpliconSize; }
	float GlobalSettings::GetMinimumTemperature() { return _minTemp; }
	float GlobalSettings::GetMaximumTemperature() { return _maxTemp; }
	float GlobalSettings::GetPrimerConcentration() { return _primerConcentration; }
	float GlobalSettings::GetMonoIonConcentration() { return _monovalentIonConcentration; }
	int GlobalSettings::GetMaximumReturnPrimers() { return _maxPrimers; }
} // End of DeGenPrime