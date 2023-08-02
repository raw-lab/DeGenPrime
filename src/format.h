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
    std::string PrintRuler(int size);

    int digits(int num);

} // End of DeGenPrime
#endif // FORMAT