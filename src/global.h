// ********************** GLOBAL VARIABLES ********************** //
// MAX_PRIMER_LENGTH (Default 24)	The largest number of 		//
//	nucleotides the program allows a primer to have.		//
// MIN_PRIMER_LENGTH (Default 18)	The smallest number of 		//
//	nucleotides the program allows a primer to have.		//
// MAX_N_TOLERANCE (Default 0.0) The maximum allowed ratio of N	//
//	degenerates allowed in a data node.  Any more than this	//
//	and the node is defaulted to 'N'.					//
// MAX_DASH_HORIZONTAL_TOLERANCE (Default 0.5) The maximum 		//
//	ratio of '-'s in a sequence.  More '-'s and the sequence	//
//	is rejected from the list.						//
// MAX_DASH_VERTICAL_TOLERANCE (Default 0.3) The maximum allowed	//
//	ratio os '-'s in a datanode.  Any more '-'s and the node	//
//	is defaulted to being a '-' node.					//
// MIN_AMPLICON_LENGTH (Default 0) The minimum amplicon size.	//
// MAX_GC_EXTREMA_RATIO (Default 0.6) The maximum ratio of G or C	//
//	allowed in the last five nucleotides.				//
// MIN_GC_TOTAL_RATIO (Default 0.4) The minimum ratio of G or C	//
//	allowed in any primer.							//
// MAX_GC_TOTAL_RATIO (Default 0.6) The maximum ratio of G or C	//
//	allowed in any primer.							//
// MIN_DEGENERATE_THRESHOLD (Default 0.3) The minimum amount of	//
//	degeneracy within a node to determine degeneracy.		//
// MIN_PRIMER_TEMP (Default 50.0) The minimum allowed temperature	//
//	of any primer in degrees celsius.					//
// MAX_PRIMER_TEMP (Default 60.0) The maximum allowed temperature	//
//	of any primer in degrees celsius.					//
// MAX_TEMP_DIFFERENCE (Default 1.0) Maximum allowed temperature	//
//	difference between forward and reverse primer.			//
// MAX_PRIMER_RETURNS (Default 5) Maximum number of primers to	//
//	return in the program output.						//
//				EXIT CODES						//
// PROGRAM_SUCCESS 0								//
// FILE_MISALIGNED 1								//
// SETTINGS_FILE_NOT_FOUND 2							//
// NO_PRIMERS_FOUND 3								//
// IMPROPER_FILE_FORMAT 4							//
// BAD_INPUT_FILE									//
// ************************************************************** //

#ifndef GLOBAL_DEGENPRIME
#define GLOBAL_DEGENPRIME

#define MIN_AMPLICON_LENGTH 0
#define MAX_PRIMER_LENGTH 22
#define MIN_PRIMER_LENGTH 18
#define MAX_GC_EXTREMA_RATIO 0.6
#define MIN_GC_TOTAL_RATIO 0.4
#define MAX_GC_TOTAL_RATIO 0.6
#define MIN_DEGENERATE_THRESHOLD 0.3
#define MAX_N_TOLERANCE 0.0
#define MAX_DASH_HORIZONTAL_TOLERANCE 0.5
#define MAX_DASH_VERTICAL_TOLERANCE 0.3
#define MIN_PRIMER_TEMP 50.0
#define MAX_PRIMER_TEMP 60.0
#define MAX_TEMP_DIFFERENCE 1.0
#define MAX_PRIMER_RETURNS 5

#define PROGRAM_SUCCESS 0
#define FILE_MISALIGNED 1
#define SETTINGS_FILE_NOT_FOUND 2
#define NO_PRIMERS_FOUND 3
#define IMPROPER_FILE_FORMAT 4
#define BAD_INPUT_FILE 5

#endif // GLOBAL_DEGENPRIME