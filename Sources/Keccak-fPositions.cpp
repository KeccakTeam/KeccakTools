/*
KeccakTools

The Keccak sponge function, designed by Guido Bertoni, Joan Daemen,
Michaël Peeters and Gilles Van Assche. For more information, feedback or
questions, please refer to our website: http://keccak.noekeon.org/

Implementation by the designers,
hereby denoted as "the implementer".

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/
*/

#include <iomanip>
#include "Keccak-fPositions.h"

BitPosition::BitPosition(unsigned int ax, unsigned int ay, unsigned int az)
: x(ax), y(ay), z(az)
{
}

BitPosition::BitPosition(const BitPosition& aPoint) 
: x(aPoint.x) , y(aPoint.y) , z(aPoint.z)
{
}

ostream& operator<<(ostream& fout, const BitPosition& point)
{
    fout << "(" << dec << point.x << "," << point.y << ","<< std::setw(2) << point.z << ") ";
    return fout;
}

bool BitPosition::operator==(const BitPosition& otherPoint) const
{
    return (x==otherPoint.x && y==otherPoint.y && z==otherPoint.z);
}

bool BitPosition::operator<(const BitPosition& otherPoint) const
{
    if (x <  otherPoint.x) return true;
    else if (x == otherPoint.x) {
        if (y <  otherPoint.y) return true;
        else if (y == otherPoint.y) {
            if (z < otherPoint.z) return true;
        }
    }
    return false;
}

bool BitPosition::nextXY()
{
    if(y<4) { 
        y++;
        return true;
    }
    if (x<4){
        y=0;
        x++;
        return true;
    }
    return false;
}

void BitPosition::set(unsigned int ax, unsigned int ay, unsigned int az)
{
    x = ax;
    y = ay;
    z = az;
}

void BitPosition::xTranslate(int dx)
{
    x = ((((int)x+dx)%5)+5)%5;
}

void BitPosition::yTranslate(int dy)
{
    y = ((((int)y+dy)%5)+5)%5;
}


void BitPosition::zTranslate(int dz, unsigned int aLaneSize)
{
    z = ((((int)z+dz)%aLaneSize)+aLaneSize)%aLaneSize;
}

bool operator<(const ColumnPosition& aCP, const ColumnPosition& bCP)
{
    return (aCP.getXplus5Z() < bCP.getXplus5Z());
}

ostream& operator<<(ostream& fout, const ColumnPosition& aCP)
{
    fout << "(" << aCP.x << ",-," << aCP.z << ") ";
    return fout;
}

bool operator<(const RowPosition& aRP, const RowPosition& bRP)
{
    return (aRP.getYplus5Z() < bRP.getYplus5Z());
}

ostream& operator<<(ostream& fout, const RowPosition& aRP)
{
    fout << "(-," << aRP.y << "," << aRP.z << ") ";
    return fout;
}
