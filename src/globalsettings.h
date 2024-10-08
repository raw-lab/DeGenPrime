// ******************************* globalsettings.h *******************************	//
// Purpose: Create a static class to evaluate and store the user settings.			//
// Mutators: Set<setting>: 	Static members to initialize a particular value of 		//
//					the class.  Compares that value to global	settings 			//
//					file to ensure a reasonable value is saved.						//
// Accessors: Get<setting>:	Returns the specified value.							//
// ********************************************************************************	//
#ifndef GLOBAL_SETTINGS
#define GLOBAL_SETTINGS

#include <fstream>

namespace DeGenPrime
{
	class GlobalSettings
	{
	public:
		GlobalSettings();

		static void SetMinimumAmplicon(int amplicon);
		static void SetBeginningNucleotide(int begin);
		static void SetDeltaG(float g);
		static void SetEndingNucleotide(int ending);
		static void SetNonDegenerate(bool non);
		static void SetMeasureByAmpliconSize(bool size);
		static void SetProteinSequence(bool seq);
		static void SetMinimumTemperature(float temp);
		static void SetMaximumTemperature(float temp);
		static void SetMaximumPrimerLength(int len);
		static void SetMinimumPrimerLength(int len);
		static void SetPrimerConcentration(float primer_conc);
		static void SetMonoIonConcentration(float salt_conc);
		static void SetMaximumReturnPrimers(int max);
		static void SetThermodynamicTemperature(float temp);
		static void SetSeqFilterLengthMin(int min);
		static void SetBeginFlag(bool begin);
		static void SetEndFlag(bool ending);
		static void SetRunTest(bool test);
		static void SetRunInvRev(bool test);
		static void SetSearchFwd(bool search);
		static void SetSearchRev(bool search);
		static void SetSortByTemp(bool temp);
		static void SetUserTemp(bool temp);
		static void SetFilterSeqs(bool fil);
		static void SetDetail(bool det);
		static void SetTestValue(std::string str);
		static void SetInvRevValue(std::string str);
		static void SetSearchFwdArg(std::string str);
		static void SetSearchRevArg(std::string str);
		static void SetInputFile(std::string str);
		static void SetOutputFile(std::string str);
		static void SetDetailFile(std::string str);
	
		static int GetMinimumAmplicon();
		static int GetBeginningNucleotide();
		static float GetDeltaG();
		static int GetEndingNucleotide();
		static bool GetNonDegenerate();
		static bool GetMeasureByAmpliconSize();
		static bool GetProteinSequence();
		static bool GetBeginFlag();
		static bool GetEndFlag();
		static float GetMinimumTemperature();
		static float GetMaximumTemperature();
		static int GetMaximumPrimerLength();
		static int GetMinimumPrimerLength();
		static float GetPrimerConcentration();
		static float GetMonoIonConcentration();
		static int GetMaximumReturnPrimers();
		static float GetThermodynamicTemperature();
		static int GetSeqFilterLengthMin();
		static bool GetRunTest();
		static bool GetRunInvRev();
		static bool GetSearchFwd();
		static bool GetSearchRev();
		static bool GetSortByTemp();
		static bool GetUserTemp();
		static bool GetFilterSeqs();
		static bool Detail();
		static std::string GetTestValue();
		static std::string GetInvRevValue();
		static std::string GetSearchFwdArg();
		static std::string GetSearchRevArg();
		static std::string GetInputFile();
		static std::string GetOutputFile();
		static std::string GetDetailFile();

	private:
		static int _ampLength;
		static int _beginningNucleotide;
		static float _deltag;
		static int _endingNucleotide;
		static bool _measureByAmpliconSize;
		static bool _proteinSequence;
		static bool _beginflag;
		static bool _endflag;
		static bool _nonDegenerate;
		static float _minTemp;
		static float _maxTemp;
		static int _maxLen;
		static int _minLen;
		static float _primerConcentration;
		static float _monovalentIonConcentration;
		static int _maxPrimers;
		static float _thermodynamicTemperature;
		static int _seqfiltermin;
		static bool _testRun;
		static bool _invRevRun;
		static bool _SearchFwd;
		static bool _SearchRev;
		static bool _sortbytemp;
		static bool _userTemp;
		static bool _filterseqs;
		static bool _detailOut;
		static std::string _testStr;
		static std::string _invRevValue;
		static std::string _searchFwdArg;
		static std::string _searchRevArg;
		static std::string _inputfile;
		static std::string _outputfile;
		static std::string _detailfile;
	};
} // end of DeGenPrime
#endif // GLOBAL_SETTINGS