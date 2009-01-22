/*
Tools for the Keccak sponge function family.
Authors: Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche

This code is hereby put in the public domain. It is given as is, 
without any guarantee.

For more information, feedback or questions, please refer to our website:
http://keccak.noekeon.org/
*/

#include <iostream>
#include "Keccak-fParts.h"

using namespace std;

SliceValue getSliceValue(RowValue row0, RowValue row1, RowValue row2, RowValue row3, RowValue row4)
{
    return
        getSliceFromRow(row0, 0)
      ^ getSliceFromRow(row1, 1)
      ^ getSliceFromRow(row2, 2)
      ^ getSliceFromRow(row3, 3)
      ^ getSliceFromRow(row4, 4);
}

unsigned int getNrActiveRows(const SliceValue& slice)
{
    return
        (getRowFromSlice(slice, 0) != 0 ? 1 : 0)
      + (getRowFromSlice(slice, 1) != 0 ? 1 : 0)
      + (getRowFromSlice(slice, 2) != 0 ? 1 : 0)
      + (getRowFromSlice(slice, 3) != 0 ? 1 : 0)
      + (getRowFromSlice(slice, 4) != 0 ? 1 : 0);
}

#if 0
PackedActiveRows getActiveRowPattern(const vector<SliceValue>& state)
{
    PackedActiveRows result = 0;
    for(unsigned int z=0; z<state.size(); z++)
        for(unsigned int y=0; y<5; y++)
            if (getRowFromSlice(state[z], y) != 0)
                result |= (PackedActiveRows)1 << (y+5*z);
    return result;
}

PackedActiveRows translateActiveRowPatternAlongZ(PackedActiveRows a, unsigned int dz, unsigned int numberOfSlices)
{
    if (dz != 0) {
        PackedActiveRows mask = (UINT32(~0)) >> (32-5*numberOfSlices);
        a &= mask;
        a = (((UINT32)a) << (5*dz)) ^ (((UINT32)a) >> (5*(numberOfSlices-dz)));
        return a & mask;
    }
    else
        return a;
}
#endif

RowValue translateRowSafely(RowValue row, int dx)
{
    dx %= 5;
    if (dx < 0) dx += 5;
    return translateRow(row, dx);
}

static const unsigned char HammingWeightTable[256] = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
};

unsigned int getHammingWeightRow(RowValue row)
{
    return HammingWeightTable[row];
}

SliceValue translateSlice(SliceValue slice, unsigned int dx, unsigned int dy)
{
    return getSliceValue(
        translateRow(getRowFromSlice(slice, (5-dy)%5), dx),
        translateRow(getRowFromSlice(slice, (6-dy)%5), dx),
        translateRow(getRowFromSlice(slice, (7-dy)%5), dx),
        translateRow(getRowFromSlice(slice, (8-dy)%5), dx),
        translateRow(getRowFromSlice(slice, (9-dy)%5), dx));
}

SliceValue translateSliceSafely(SliceValue slice, int dx, int dy)
{
    dx %= 5;
    if (dx < 0) dx += 5;
    dy %= 5;
    if (dy < 0) dy += 5;
    return translateSlice(slice, dx, dy);
}

unsigned int getHammingWeightSlice(SliceValue slice)
{
    return HammingWeightTable[slice & 0xFF]
        + HammingWeightTable[(slice >>  8) & 0xFF]
        + HammingWeightTable[(slice >> 16) & 0xFF]
        + HammingWeightTable[(slice >> 24) & 0xFF];
}

unsigned int getHammingWeight(const vector<SliceValue>& state)
{
    unsigned int result = 0;
    for(unsigned int z=0; z<state.size(); z++)
        result += getHammingWeightSlice(state[z]);
    return result;
}

void displayPlane(ostream& fout, const vector<SliceValue>& state, int offset, unsigned int y)
{
    for(unsigned int z=0; z<state.size(); z++) {
        RowValue row = getRowFromSlice(state[z], y);
        for(unsigned int sx=0; sx<5; sx++) {
            unsigned int x = KeccakF::index(sx-offset);
            if ((row & (1<<x)) != 0)
                fout << "X";
            else
                if ((x == 0) && (y == 0) && (z == 0))
                    fout << "+";
                else
                    fout << ".";
        }
        if (z < (state.size()-1))
            fout << "   ";
    }
}

void displayNothing(ostream& fout, const vector<SliceValue>& state)
{
    for(unsigned int z=0; z<state.size(); z++) {
        fout << "     "; // 5 spaces
        if (z < (state.size()-1))
            fout << "   ";
    }
}

void displayParity(ostream& fout, const vector<SliceValue>& state, const int& offset)
{
    for(unsigned int z=0; z<state.size(); z++) {
        RowValue parity = 0;
        for(unsigned int y=0; y<5; y++)
            parity ^= getRowFromSlice(state[z], y);
        for(unsigned int sx=0; sx<5; sx++) {
            unsigned int x = KeccakF::index(sx-offset);
            if ((parity & (1<<x)) != 0)
                fout << "O";
            else
                fout << "-";
        }
        if (z < (state.size()-1))
            fout << "   ";
    }
}

void displayState(ostream& fout, const vector<SliceValue>& state, bool showParity)
{
    const int offset = 2;
    for(unsigned int sy=0; sy<5; sy++) {
        unsigned int y = KeccakF::index(-1-sy-offset);
        displayPlane(fout, state, offset, y);
        fout << endl;
    }
    if (showParity) {
        displayParity(fout, state, offset);
        fout << endl;
    }
}

void displaySlice(ostream& fout, SliceValue slice)
{
    vector<SliceValue> slices;
    slices.push_back(slice);
    displayState(fout, slices);
}

void displayStates(ostream& fout,
                   const vector<SliceValue>& state1, bool showParity1,
                   const vector<SliceValue>& state2, bool showParity2)
{
    const int offset = 2;
    for(unsigned int sy=0; sy<5; sy++) {
        unsigned int y = KeccakF::index(-1-sy-offset);
        displayPlane(fout, state1, offset, y);
        fout << "  |  "; 
        displayPlane(fout, state2, offset, y);
        fout << endl;
    }
    if (showParity1 || showParity2) {
        if (showParity1)
            displayParity(fout, state1, offset);
        else
            displayNothing(fout, state1);
        if (showParity2) {
            fout << "     ";
            displayParity(fout, state2, offset);
        }
        fout << endl;
    }
}

void displayStates(ostream& fout,
                   const vector<SliceValue>& state1, bool showParity1,
                   const vector<SliceValue>& state2, bool showParity2,
                   const vector<SliceValue>& state3, bool showParity3)
{
    const int offset = 2;
    for(unsigned int sy=0; sy<5; sy++) {
        unsigned int y = KeccakF::index(-1-sy-offset);
        displayPlane(fout, state1, offset, y);
        fout << "  |  "; 
        displayPlane(fout, state2, offset, y);
        fout << "  |  ";
         displayPlane(fout, state3, offset, y);
        fout << endl;
    }
    if (showParity1 || showParity2 || showParity3) {
        if (showParity1)
            displayParity(fout, state1, offset);
        else
            displayNothing(fout, state1);
        fout << "     ";
        if (showParity2)
            displayParity(fout, state2, offset);
        else
            displayNothing(fout, state2);
        if (showParity3) {
            fout << "     ";
            displayParity(fout, state3, offset);
        }
        fout << endl;
    }
}

void translateStateAlongZ(vector<SliceValue>& state, unsigned int dz)
{
    if (dz != 0) {
        vector<SliceValue> stateBis(state);
        for(unsigned int z=0; z<stateBis.size(); z++) {
            unsigned int newZ = (z+dz)%stateBis.size();
            state[newZ] = stateBis[z];
        }
    }
}

RowValue getRow(const vector<UINT64>& lanes, unsigned int y, unsigned int z)
{
    RowValue result = 0;
    for(unsigned int x=0; x<5; x++) {
        if ((lanes[KeccakF::index(x, y)] & ((UINT64)1 << z)) != 0)
            result ^= (1 << x);
    }
    return result;
}

void setRow(vector<UINT64>& lanes, RowValue row, unsigned int y, unsigned int z)
{
    for(unsigned int x=0; x<5; x++)
        if ((row & (1 << x)) != 0)
            lanes[KeccakF::index(x, y)] |= (UINT64)1 << z;
        else
            lanes[KeccakF::index(x, y)] &= ~((UINT64)1 << z);
}

SliceValue getSlice(const vector<UINT64>& lanes, unsigned int z)
{
    SliceValue result = 0;
    for(unsigned int y=0; y<5; y++)
        result ^= getSliceFromRow(getRow(lanes, y, z), y);
    return result;
}

void setSlice(vector<UINT64>& lanes, SliceValue slice, unsigned int z)
{
    for(unsigned int y=0; y<5; y++)
        setRow(lanes, getRowFromSlice(slice, y), y, z);
}

void fromLanesToSlices(vector<SliceValue>& slices, const vector<UINT64>& lanes, unsigned int laneSize)
{
    slices.resize(laneSize);
    for(unsigned int z=0; z<laneSize; z++)
        slices[z] = getSlice(lanes, z);
}

void fromSlicesToLanes(vector<UINT64>& lanes, const vector<SliceValue>& slices)
{
    lanes.assign(25, 0);
    for(unsigned int z=0; z<slices.size(); z++)
        setSlice(lanes, slices[z], z);
}
