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

#include <fstream>
#include <math.h>
#include "Keccak-f25LUT.h"

using namespace std;

KeccakF25LUT::KeccakF25LUT(unsigned int aNrRounds)
    : KeccakFfirstRounds(25, aNrRounds)
{
    generateLUT();
}

KeccakF25LUT::KeccakF25LUT()
    : KeccakFfirstRounds(25)
{
    generateLUT();
}

void KeccakF25LUT::generateLUT()
{
   if (!retrieveLUT()) {
        cout << "Generating the look-up table..." << flush;
        vector<LaneValue> lanes(25, 0);
        LUT.resize(1<<25);
        for (SliceValue sliceIn=0; sliceIn< (1<<25) ; sliceIn++) {
            if ((sliceIn & 0xfffff) == 0) cout << " " << floor(sliceIn*100.0/(1<<25)) << "%" << flush;
            setSlice(lanes, sliceIn);
            forward(lanes);
            SliceValue sliceOut = getSlice(lanes);
            LUT[sliceIn] = sliceOut;
        }
        cout << " done, now saving to disk..." << flush;
        saveLUT();
        cout << "and saved." << endl;
    }
}

void KeccakF25LUT::saveLUT() const
{
    string fileName = buildFileName("", ".LUT");
    ofstream fout(fileName.c_str(), ios::binary);
    for (SliceValue i=0 ; i<(1<<25); i++) {
        static unsigned char tmp[4];
        tmp[0] =  LUT[i]&0xFF;
        tmp[1] = (LUT[i]>>8)&0xFF;
        tmp[2] = (LUT[i]>>16)&0xFF;
        tmp[3] = (LUT[i]>>24)&0xFF;
        fout.write((char *)tmp, 4);
    }
}

bool KeccakF25LUT::retrieveLUT()
{
    string fileName = buildFileName("", ".LUT");
    ifstream fin(fileName.c_str(), ios::binary);
    if (!fin) return false;
    LUT.resize(1<<25);
    for (SliceValue i=0 ; i<(1<<25); i++) {
        static unsigned char tmp[4];
        fin.read((char *)tmp, 4);
        LUT[i]  = tmp[3];  LUT[i] <<= 8;
        LUT[i] ^= tmp[2];  LUT[i] <<= 8;
        LUT[i] ^= tmp[1];  LUT[i] <<= 8;
        LUT[i] ^= tmp[0];
    }
    return true;
}
