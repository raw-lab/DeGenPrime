#include <cmath>
#include <string>
#include "format.h"
#include "global.h"

using namespace std;

namespace DeGenPrime
{
    string Format(std::string str, int size, Alignment align)
    {
        string ret = str;
        int len = str.length();
        if(size > STR_FORMAT)
        {
            size = STR_FORMAT;
        }
        if(len >= size)
        {
            return str.substr(0,size);
        }

        int diff = size - len;
        switch(align)
        {
            case Alignment::Left:
            case Alignment::Right:
            {
                string fill = "";
                for(int i = 0;i < diff;i++)
                {
                    fill += " ";
                }
                if(align == Alignment::Left)
                {
                    ret.append(fill);
                }
                else
                {
                    fill.append(ret);
                    ret = fill;
                }
                break;
            }
            case Alignment::Center:
            {
                string fwd_fill = "";
                string rev_fill = "";
                int half = diff / 2;
                bool odd = (diff % 2 == 1) ? true : false;
                for(int i = 0;i < half;i++)
                {
                    fwd_fill += " ";
                    rev_fill += " ";
                }
                if(odd)fwd_fill += " ";
                ret.append(rev_fill);
                fwd_fill.append(ret);
                ret = fwd_fill;
                break;
            }
            default:
                return str;
        }
        return ret;
    }

    string Format(int iStr, int size)
    {
        // Add zeros to the front of the digit until
        // The size of the string equals what is desired.
        string ret = "";
        int zeros = size - digits(iStr);
        if(zeros > 0)
        {
            for(int i = 0;i < zeros;i++)
            {
                ret += "0";
            }
        }
        ret += to_string(iStr);
        return ret;
    }

    string Format(float f, int precision)
    {
        // Floating numbers format
        // if less than zero, include leading zero
        // precision is number after decimal
        //  '0.blah'
        bool neg = (f < 0);
        string ret = neg ? "-" : "";
        int precision2 = neg ? precision - 1 : precision;
        float f2 = abs(f);
        int mult = pow(10, precision);
        int temp = round(f2 * mult);
        int before_dec = temp / mult;
        int after_dec = temp % mult;
        string before = Format(before_dec, precision2);
        string after = Format(after_dec, precision);
        ret += before + "." + after;
        return ret;
    }

    int digits(int num)
    {
        int ret = 0;
        int temp = num;
        if(temp == 0)
        {
            ret = 1;
        }
        else
        {
            while(temp > 0)
            {
                temp /= 10;
                ret++;
            }
        }
        return ret;
    }

} // End of DeGenPrime