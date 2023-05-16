// ******************************* GlobalSettings.h *******************************	//
// Purpose: Create a static class to evaluate and store the user settings.		//
// Mutators: Set<setting>: 	Static members to initialize a particular value of 	//
//					the class.  Compares that value to global	settings 	//
//					file to ensure a reasonable value is saved.		//
// Accessors: Get<setting>:	Returns the specified value.					//
// ********************************************************************************	//
#ifndef GLOBAL_SETTINGS
#define GLOBAL_SETTINGS

namespace DeGenPrime
{
	class GlobalSettings
	{
	public:
		// GlobalSettings();

		static void SetMinimumAmplicon(int amplicon);
		static void SetBeginningNucleotide(int begin);
		static void SetEndingNucleotide(int ending);
		static void SetMeasureByAmpliconSize(bool size);
		static void SetMinimumTemperature(float temp);
		static void SetMaximumTemperature(float temp);
		static void SetPrimerConcentration(float primer_conc);
		static void SetMonoIonConcentration(float salt_conc);
		static void SetMaximumReturnPrimers(int max);
	
		static int GetMinimumAmplicon();
		static int GetBeginningNucleotide();
		static int GetEndingNucleotide();
		static bool GetMeasureByAmpliconSize();
		static float GetMinimumTemperature();
		static float GetMaximumTemperature();
		static float GetPrimerConcentration();
		static float GetMonoIonConcentration();
		static int GetMaximumReturnPrimers();
	private:
		static int _ampLength;
		static int _beginningNucleotide;
		static int _endingNucleotide;
		static bool _measureByAmpliconSize;
		static float _minTemp;
		static float _maxTemp;
		static float _primerConcentration;
		static float _monovalentIonConcentration;
		static int _maxPrimers;
	};
} // end of DeGenPrime
#endif // GLOBAL_SETTINGS
