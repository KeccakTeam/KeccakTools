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

#include <sstream>
#include "Keccak-f.h"

using namespace std;

KeccakException::KeccakException(const string& aReason)
    : reason(aReason)
{
}


KeccakF::KeccakF(unsigned int aWidth, unsigned int aNrRounds)
{
    width = aWidth;
    laneSize = width/25;
    switch(width) {
    case 25:
        nominalNrRounds = 12;
        break;
    case 50:
        nominalNrRounds = 14;
        break;
    case 100:
        nominalNrRounds = 16;
        break;
    case 200:
        nominalNrRounds = 18;
        break;
    case 400:
        nominalNrRounds = 20;
        break;
    case 800:
        nominalNrRounds = 22;
        break;
    case 1600:
        nominalNrRounds = 24;
        break;
    default:
        throw(KeccakException("The width must be 25 times a power of two between 1 and 64."));
    }
    nrRounds = (aNrRounds == 0) ? nominalNrRounds : aNrRounds;
    mask = (LaneValue(~0)) >> (64-laneSize);
    initializeRhoOffsets();
    initializeRoundConstants();
}

unsigned int KeccakF::getWidth() const
{
    return width;
}

unsigned int KeccakF::getLaneSize() const
{
    return laneSize;
}

unsigned int KeccakF::index(int x, int y)
{
    x %= 5;
    if (x<0) x += 5;
    y %= 5;
    if (y<0) y += 5;

    return(x + (5*y));
}

unsigned int KeccakF::index(int x)
{
    x %= 5;
    if (x<0) x += 5;
    return x;
}

void KeccakF::pi(unsigned int x, unsigned int y, unsigned int& X, unsigned int& Y)
{
    X = (0*x + 1*y)%5;
    Y = (2*x + 3*y)%5;
}

void KeccakF::inversePi(unsigned int X, unsigned int Y, unsigned int& x, unsigned int& y)
{
    x = (1*X + 3*Y)%5;
    y = (1*X + 0*Y)%5;
}

void KeccakF::ROL(LaneValue& L, int offset) const
{
    offset %= (int)laneSize;
    if (offset < 0) offset += laneSize;

    if (offset != 0) {
        L &= mask;
        L = (((LaneValue)L) << offset) ^ (((LaneValue)L) >> (laneSize-offset));
    }
    L &= mask;
}

void KeccakF::fromBytesToLanes(const UINT8 *in, vector<LaneValue>& out) const
{
    out.resize(25);
    if ((laneSize == 1) || (laneSize == 2) || (laneSize == 4) || (laneSize == 8)) {
        for(unsigned int i=0; i<25; i++)
            out[i] = (in[i*laneSize/8] >> ((i*laneSize) % 8)) & mask;
    }
    else if ((laneSize == 16) || (laneSize == 32) || (laneSize == 64)) {
        for(unsigned int i=0; i<25; i++) {
            out[i] = 0;
            for(unsigned int j=0; j<(laneSize/8); j++)
                out[i] |= LaneValue((UINT8)(in[i*laneSize/8+j])) << (8*j);
        }
    }
}

void KeccakF::fromLanesToBytes(const vector<LaneValue>& in, UINT8 *out) const
{
    if ((laneSize == 1) || (laneSize == 2) || (laneSize == 4) || (laneSize == 8)) {
        for(unsigned int i=0; i<(25*laneSize+7)/8; i++)
            out[i] = 0;
        for(unsigned int i=0; i<25; i++)
            out[i*laneSize/8] |= in[i] << ((i*laneSize) % 8);
    }
    else if ((laneSize == 16) || (laneSize == 32) || (laneSize == 64)) {
        for(unsigned int i=0; i<25; i++)
            for(unsigned int j=0; j<(laneSize/8); j++)
                out[i*(laneSize/8)+j] = (UINT8)((in[i] >> (8*j)) & 0xFF);
    }
}

void KeccakF::operator()(UINT8 * state) const
{
    vector<LaneValue> A(25);
    fromBytesToLanes(state, A);
    forward(A);
    fromLanesToBytes(A, state);
}

void KeccakF::inverse(UINT8 * state) const
{
    vector<LaneValue> A(25);
    fromBytesToLanes(state, A);
    inverse(A);
    fromLanesToBytes(A, state);
}

string KeccakF::getDescription() const
{
    stringstream a;
    a << "Keccak-f[" << dec << width;
    if (nrRounds != nominalNrRounds)
        a << ", nr=" << dec << nrRounds;
    a << "]";
    return a.str();
}

string KeccakF::getName() const
{
    stringstream a;
    a << "KeccakF-" << dec << width << "-" << nrRounds;
    return a.str();
}

string KeccakF::buildFileName(const string& prefix, const string& suffix) const
{
    stringstream fileName;
    fileName << prefix << getName() << suffix;
    return fileName.str();
}

bool LFSR86540(UINT8& state)
{
    bool result = (state & 0x01) != 0;
    if ((state & 0x80) != 0)
        state = (state << 1) ^ 0x71;
    else
        state <<= 1;
    return result;
}

void KeccakF::initializeRoundConstants()
{
    roundConstants.clear();
    UINT8 LFSRstate = 0x01;
    for(unsigned int i=0; i<nrRounds; i++) {
        LaneValue c = 0;
        for(unsigned int j=0; j<7; j++) {
            unsigned int bitPosition = (1<<j)-1; //2^j-1
            if (LFSR86540(LFSRstate))
                c ^= (LaneValue)1<<bitPosition;
        }
        roundConstants.push_back(c & mask);
    }
}

void KeccakF::initializeRhoOffsets()
{
    rhoOffsets.resize(25);
    rhoOffsets[index(0, 0)] = 0;
    unsigned int x = 1;
    unsigned int y = 0;
    for(unsigned int t=0; t<24; t++) {
        rhoOffsets[index(x, y)] = ((t+1)*(t+2)/2) % laneSize;
        unsigned int newX = (0*x + 1*y)%5;
        unsigned int newY = (2*x + 3*y)%5;
        x = newX;
        y = newY;
    }
}

string KeccakF::buildBitName(const string& prefixSymbol, unsigned int laneSize, unsigned int z)
{
    stringstream s;
    s << prefixSymbol;
    if (laneSize > 1) {
        s.fill('0');
        if (laneSize <= 10)
            s.width(1);
        else if (laneSize <= 100)
            s.width(2);
        s << z;
    }
    return s.str();
}

string KeccakF::bitName(const string& prefix, unsigned int x, unsigned int y, unsigned int z) const
{
    return buildBitName(laneName(prefix, x, y), laneSize, z);
}

string KeccakF::laneName(const string& prefix, unsigned int x, unsigned int y)
{
    return prefix + "bgkms"[y] + "aeiou"[x];
}

string KeccakF::sheetName(const string& prefix, unsigned int x)
{
    return prefix + "aeiou"[x];
}
