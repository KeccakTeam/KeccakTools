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

#include "Keccak-fParity.h"

RowValue getParity(SliceValue slice)
{
    RowValue parity = 0;
    for(unsigned int y=0; y<nrRowsAndColumns; y++)
        parity ^= getRowFromSlice(slice, y);
    return parity;
}

PackedParity getParity(const vector<SliceValue>& state)
{
    PackedParity parities = 0;
    for(unsigned int z=0; z<state.size(); z++)
        parities ^= getPackedParityFromParity(getParity(state[z]), z);
    return parities;
}

PackedParity packParity(const vector<RowValue>& parity)
{
    PackedParity result = 0;
    for(unsigned int z=0; z<parity.size(); z++)
        result ^= getPackedParityFromParity(parity[z], z);
    return result;
}

void unpackParity(PackedParity packedParity, vector<RowValue>& parity, unsigned int laneSize)
{
    if (parity.size() != laneSize)
        parity.resize(laneSize);
    for(unsigned int z=0; z<laneSize; z++)
        parity[z] = getParityFromPackedParity(packedParity, z);
}

void getParity(const vector<SliceValue>& state, vector<RowValue>& parity)
{
    parity.resize(state.size());
    for(unsigned int z=0; z<state.size(); z++)
        parity[z] = getParity(state[z]);
}

void getParity(const vector<LaneValue>& state, vector<LaneValue>& parity)
{
    if (parity.size() != nrRowsAndColumns)
        parity.resize(nrRowsAndColumns);
    for(unsigned int x=0; x<nrRowsAndColumns; x++) {
        parity[x] = state[KeccakF::index(x,0)]; 
        for(unsigned int y=1; y<5; y++) 
            parity[x] ^= state[KeccakF::index(x,y)];
    }
}

void fromSlicesToSheetsParity(const vector<RowValue>& paritySlices, vector<LaneValue>& paritySheet)
{
    paritySheet.assign(nrRowsAndColumns, 0);
    for(unsigned int z=0; z<paritySlices.size(); z++)
        for(unsigned int x=0; x<nrRowsAndColumns; x++) {
            if ((paritySlices[z] & ((RowValue)1 << x)) != 0)
                paritySheet[x] ^= (LaneValue)1 << z;
        }
}

void fromSheetsToSlicesParity(const vector<LaneValue>& paritySheet, vector<RowValue>& paritySlices)
{
    for(unsigned int z=0; z<paritySlices.size(); z++) {
        paritySlices[z] = 0;
        for(unsigned int x=0; x<nrRowsAndColumns; x++) {
            if ((paritySheet[x] & ((LaneValue)1 << z)) != 0)
                paritySlices[z] ^= (RowValue)1 << x;
        }
    }
}

string getDisplayOfParity(RowValue C, RowValue D)
{
    string result;
    const int offset = 2;
    for(unsigned int sx=0; sx<5; sx++) {
        unsigned int x = KeccakF::index(sx-offset);
        bool affected = (D & (1 << x)) != 0;
        bool odd = (C & (1 << x)) != 0;
        if (affected)
            if (odd)
                result += "!";
            else
                result += "|";
        else
            if (odd)
                result += ".";
            else
                result += "-";
    }
    return result;
}

void displayParity(ostream& fout, RowValue C, RowValue D)
{
    fout << getDisplayOfParity(C, D) << endl;
}

void displayParity(ostream& fout, const vector<RowValue>& C, const vector<RowValue>& D)
{
    unsigned int z = 0;
    unsigned int laneSize = C.size();
    while(z < laneSize) {
        unsigned int zeroes = 0;
        while((z < laneSize) && (C[z] == 0) && (D[z] == 0)) {
            z++;
            zeroes++;
        }
        if (zeroes >= 2)
            fout << "  z^" << dec << zeroes << endl;
        else {
            for(unsigned int iz=z-zeroes; iz<z; iz++)
                displayParity(fout, C[iz], D[iz]);
        }
        if (z < laneSize) {
            displayParity(fout, C[z], D[z]);
            z++;
        }
    }
}

void writeParity(ostream& out, const vector<RowValue>& C)
{
    out << hex << C.size() << " ";
    for(unsigned int z=0; z<C.size(); z++)
        out << hex << (int)C[z] << " ";
    out << endl;
}

void readParity(istream& in, vector<RowValue>& C)
{
    unsigned int laneSize;
    in >> hex >> laneSize;
    C.resize(laneSize);
    for(unsigned int z=0; z<C.size(); z++) {
        int data;
        in >> hex >> data;
        C[z] = data;
    }
}
