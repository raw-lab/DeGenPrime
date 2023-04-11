// ********************** GLOBAL VARIABLES ********************** //
// MAX_PRIMER_LENGTH (Default 24)	The largest number of 		//
//	nucleotides the program allows a primer to have.		//
// MIN_PRIMER_LENGTH (Default 18)	The smallest number of 		//
//	nucleotides the program allows a primer to have.		//
// MAX_N_TOLERANCE (Default 0.0) The maximum allowed ratio of N	//
//	degenerates allowed in a data node.  Any more than this	//
//	and the node is defaulted to 'N'.					//
// MAX_-_HORIZONTAL_TOLERANCE (Default 0.5) The maximum allowed	//
//	ratio of '-'s in a sequence.  More '-'s and the sequence	//
//	is rejected from the list.						//
// MAX_-_VERTICAL_TOLERANCE (Default 0.3) The maximum allowed	//
//	ratio os '-'s in a datanode.  Any more '-'s and the node	//
//	is defaulted to being a '-' node.					//
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

#endif // GLOBAL_DEGENPRIME