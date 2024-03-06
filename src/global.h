// ********************** GLOBAL VARIABLES ************************ //
// MAX_PRIMER_LENGTH (Default 25)	The largest number of 		    //
//	nucleotides the program allows a primer to have.	        	//
// MIN_PRIMER_LENGTH (Default 18)	The smallest number of 	    	//
//	nucleotides the program allows a primer to have.	        	//
// MAX_N_TOLERANCE (Default 0.0) The maximum allowed ratio of N	    //
//	degenerates allowed in a data node.  Any more than this	        //
//	and the node is defaulted to 'N'.			            		//
// MAX_DASH_HORIZONTAL_TOLERANCE (Default 0.5) The maximum 		    //
//	ratio of '-'s in a sequence.  More '-'s and the sequence        //
//	is rejected from the list.					                	//
// MAX_DASH_VERTICAL_TOLERANCE (Default 0.3) The maximum allowed	//
//	ratio os '-'s in a datanode.  Any more '-'s and the node	    //
//	is defaulted to being a '-' node.				            	//
// MAX_PENALTY (Default 100.0) The maximum allowed penalty score    //
//  for primers.                                                    //
// MAX_GC_EXTREMA_RATIO (Default 0.6) The maximum ratio of G or C	//
//	allowed in the last five nucleotides.			               	//
// MIN_GC_TOTAL_RATIO (Default 0.4) The minimum ratio of G or C	    //
//	allowed in any primer.						                	//
// MAX_GC_TOTAL_RATIO (Default 0.6) The maximum ratio of G or C	    //
//	allowed in any primer.						                	//
// MIN_DEGENERATE_THRESHOLD (Default 60.0) The minimum amount of 	//
//	degeneracy within a node to determine degeneracy.		        //
// MIN_PRIMER_TEMP (Default 50.0) The minimum allowed temperature	//
//	of any primer in degrees celsius.			            		//
// MAX_PRIMER_TEMP (Default 60.0) The maximum allowed temperature	//
//	of any primer in degrees celsius.				            	//
// MAX_TEMP_DIFFERENCE (Default 1.0) Maximum allowed temperature	//
//	difference between forward and reverse primer.			        //
// MAX_PRIMER_RETURNS (Default 10) Maximum number of primers to	    //
//	return in the program output.						            //
// MIN_SALT_CONC (Default 50.0) Minimum salt concentration in mM	//
//	allowed by the program						                	//
// MIN_PRIMER_CONC (Default 50.0) Minimum primer concentration in	//
//	nM allowed by the program.				                		//
//									                        		//
//				EXIT CODES					                    	//
// PROGRAM_SUCCESS 0					                			//
// FILE_MISALIGNED 1					                			//
// SETTINGS_FILE_NOT_FOUND 2			               				//
// NO_PRIMERS_FOUND 3						                		//
// IMPROPER_FILE_FORMAT 4						                	//
// BAD_INPUT_FILE	5							                	//
//				DEFAULT SETTINGS					                //
// DEFAULT_AMPLICON_LENGTH (Default 0) Default minimum amplicon	    //
// DEFAULT_BEGIN_NUCLEOTIDE (Default 0) Default minimum		        //
//	nucleotide for the forward primer. Note that range	        	//
//	measuring of primers is false by default.			        	//
// DEFAULT_END_NUCLEOTIDE (Default 0) Default maximum nucleotide	//
//	for the end of the forward primer. Note that range		        //
//	measuring of primers is false by default.				        //
// DEFAULT_MEASURE_BY_AMPLICON (Default true) Method by which	    //
//	primers are chosen.  Can either by by minimum amplicon	        //
//	length (true) or range (false).					                //
// DEFAULT_PROTEIN_SEQUENCE (Default false) Method by which		    //
//	program decides if an input file is a protein sequence.     	//
// DEFAULT_MIN_PRIMER_LENGTH (Default 20) Default minimum size of   //
//  primer.  User can lower this to at least 18.                    //
// DEFAULT_MAX_PRIMER_LENGTH (Default 22) Default maximum size of   //
//  primer.  User can raise this to at max 25.                      //
// DEFAULT_MIN_TEMP (Default 50.0) Default minimum melting		    //
//	temperature to find primers.						            //
// DEFAULT_MAX_TEMP (Default 60.0) Default maximum melting	    	//
//	temperature to find primers.					            	//
// DEFAULT_PRIMER_CONC (Default 50.0) Default concentration of  	//
//	primers for the PCR reaction.  Measured in nanoMolar		    //
// DEFAULT_SALT_CONC (Default 50.0) Default concentration of	    //
//	monovalent ions in the PCR reaction. Measured in milliMolar	    //
// DEFAULT_MAX_PRIMERS (Default 5) Maximum number of primers to	    //
//	return in the program output.					            	//
// DEFAULT_THERMODYNAMIC_TEMPERATURE (Default 37.0) The default	    //
//	temperature in celsius which thermodynamic calculations     	//
//	are made.									                    //
// **************************************************************** //

#ifndef GLOBAL_DEGENPRIME
#define GLOBAL_DEGENPRIME

#define MAX_PRIMER_LENGTH 25
#define MIN_PRIMER_LENGTH 18
#define MAX_GC_EXTREMA_RATIO 0.7
#define MIN_GC_TOTAL_RATIO 0.4
#define MAX_GC_TOTAL_RATIO 0.6001
#define MIN_DEGENERATE_THRESHOLD 60
#define MAX_N_TOLERANCE 0.0
#define MAX_DASH_HORIZONTAL_TOLERANCE 0.5
#define MAX_DASH_VERTICAL_TOLERANCE 0.3
#define MAX_PENALTY 100.0
#define MIN_PRIMER_TEMP 50.0
#define MAX_PRIMER_TEMP 65.0
#define MAX_TEMP_DIFFERENCE 1.0
#define MAX_PRIMER_RETURNS 5
#define MIN_SALT_CONC 10.0
#define MIN_PRIMER_CONC 25.0
#define STR_FORMAT 76

#define PROGRAM_SUCCESS 0
#define FILE_MISALIGNED 1
#define SETTINGS_FILE_NOT_FOUND 2
#define NO_PRIMERS_FOUND 3
#define IMPROPER_FILE_FORMAT 4
#define BAD_INPUT_FILE 5
#define TEST_MODE 6

#define DEFAULT_AMPLICON_LENGTH 100
#define DEFAULT_BEGIN_NUCLEOTIDE 0
#define DEFAULT_DELTA_G -4.0
#define DEFAULT_END_NUCLEOTIDE 2147483647
#define DEFAULT_MEASURE_BY_AMPLICON true
#define DEFAULT_PROTEIN_SEQUENCE false
#define DEFAULT_MAX_PRIMER_LENGTH 22
#define DEFAULT_MIN_PRIMER_LENGTH 20
#define DEFAULT_MIN_TEMP 50.0
#define DEFAULT_MAX_TEMP 60.0
#define DEFAULT_PRIMER_CONC 50.0
#define DEFAULT_SALT_CONC 50.0
#define DEFAULT_MAX_PRIMERS 5
#define DEFAULT_THERMODYNAMIC_TEMPERATURE 37.0
#define DEFAULT_BEGIN_FLAG false
#define DEFAULT_END_FLAG false
#define DEFAULT_RUN_TEST false
#define DEFAULT_TEST_VAL false

#endif // GLOBAL_DEGENPRIME