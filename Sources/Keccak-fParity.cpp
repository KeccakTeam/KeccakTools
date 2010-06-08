/*
Tools for the Keccak sponge function family.
Authors: Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche

This code is hereby put in the public domain. It is given as is, 
without any guarantee.

For more information, feedback or questions, please refer to our website:
http://keccak.noekeon.org/
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
