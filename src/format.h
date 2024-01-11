// ******************************* format.h *******************************	//
// Purpose: Streamline output of strings and formatting	                    //
// Enum: Alignment: How to align a string; Left, Right, or Center        	//
// Functions: string Format(string str, int total_size, Alignment): returns //
//          a string with size equal to total_size,  aligns the text in the //
//          first argument to the alignment in the third argument.  Works   //
//          by filling in spaces.                                           //
//      Format(int, int): Returns a with size of the second argument.  The  //
//          integer in the first argument has leading zeros added until     //
//          it reaches that size.                                           //
//      Format(float, int): Returns a string with the digits before and     //
//          after the decimal point in the first argument are equal to the  //
//          second argument.  Adds leading zeros ahead of integer and       //
//          rounds the number after the decimal to a whole number.          //
//      digits(int): Returns the number of digits in the argument.          //
// ************************************************************************ //

#ifndef FORMAT
#define FORMAT
#include <string>

namespace DeGenPrime
{

    enum Alignment
    {
        Left,
        Right,
        Center
    };

    std::string Format(std::string str, int size, Alignment align);
    std::string Format(int iStr, int size);
    std::string Format(float f, int precision);

    int digits(int num);

} // End of DeGenPrime
#endif // FORMAT