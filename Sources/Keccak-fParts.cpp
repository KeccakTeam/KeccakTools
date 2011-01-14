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

#include <iostream>
#include <sstream>
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

unsigned int getNrActiveRows(const vector<SliceValue>& slices)
{
    unsigned int a = 0;
    for(unsigned int i=0; i<slices.size(); i++)
        a += getNrActiveRows(slices[i]);
    return a;
}

unsigned int getNrActiveRows(const vector<LaneValue>& lanes)
{
    unsigned int result = 0;
    for(unsigned int y=0; y<5; y++) {
        LaneValue m = lanes[KeccakF::index(0, y)] 
            | lanes[KeccakF::index(1, y)] 
            | lanes[KeccakF::index(2, y)] 
            | lanes[KeccakF::index(3, y)] 
            | lanes[KeccakF::index(4, y)];
        result += getHammingWeightLane(m);
    }
    return result;
}


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

unsigned int getHammingWeightLane(LaneValue lane)
{
    return HammingWeightTable[lane & 0xFF]
        + HammingWeightTable[(lane >>  8) & 0xFF]
        + HammingWeightTable[(lane >> 16) & 0xFF]
        + HammingWeightTable[(lane >> 24) & 0xFF]
        + HammingWeightTable[(lane >> 32) & 0xFF]
        + HammingWeightTable[(lane >> 40) & 0xFF]
        + HammingWeightTable[(lane >> 48) & 0xFF]
        + HammingWeightTable[(lane >> 56) & 0xFF];
}

unsigned int getHammingWeight(const vector<LaneValue>& state)
{
    unsigned int result = 0;
    for(unsigned int z=0; z<state.size(); z++)
        result += getHammingWeightLane(state[z]);
    return result;
}

void getDisplayMap(const vector<SliceValue>& state, vector<unsigned int>& displayMap)
{
    bool optimizeEmptySlices = (state.size() >= 8);
    if (optimizeEmptySlices) {
        unsigned int z=0;
        while(z < state.size()) {
            while((z < state.size()) && (state[z] != 0))
                z++;
            displayMap.push_back(z);
            while((z < state.size()) && (state[z] == 0))
                z++;
            displayMap.push_back(z);
        }
    }
    else {
        displayMap.push_back(state.size());
        displayMap.push_back(state.size());
    }
}

inline unsigned int getDisplayNumberOfSpaces(const unsigned int& delta)
{
    return (delta <= 3) ? delta : ((delta < 10) ? 3 : 4);
}

void displayPlane(ostream& fout, const vector<SliceValue>& state, int offset, unsigned int y, const vector<unsigned int>& displayMap)
{
    unsigned int z=0;
    for(unsigned int i=0; i<displayMap.size(); i+=2) {
        for( ; z<displayMap[i]; z++) {
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
            if (z < (displayMap[i]-1))
                fout << "   ";
            else if (z < (state.size()-1))
                fout << " ";
        }
        if (displayMap[i] < displayMap[i+1]) {
            int delta = displayMap[i+1] - displayMap[i];
            if (y == 0) {
                ostringstream sout;
                if (delta == 1)
                    sout << "z";
                else if (delta == 2)
                    sout << "zz";
                else
                    sout << "z^" << dec << delta;
                fout << sout.str();
            }
            else {
                unsigned int numberOfSpaces = getDisplayNumberOfSpaces(displayMap[i+1] - displayMap[i]);
                for(unsigned int j=0; j<numberOfSpaces; j++)
                    fout << " ";
            }
            z = displayMap[i+1];
            if (z < state.size())
                fout << " ";
        }
    }
}

void displayNothing(ostream& fout, const vector<SliceValue>& state, const vector<unsigned int>& displayMap)
{
    unsigned int z=0;
    for(unsigned int i=0; i<displayMap.size(); i+=2) {
        for( ; z<displayMap[i]; z++) {
            fout << "     "; // 5 spaces
            if (z < (displayMap[i]-1))
                fout << "   ";
            else if (z < (state.size()-1))
                fout << " ";
        }
        if (displayMap[i] < displayMap[i+1]) {
            unsigned int numberOfSpaces = getDisplayNumberOfSpaces(displayMap[i+1] - displayMap[i]);
            for(unsigned int j=0; j<numberOfSpaces; j++)
                fout << " ";
            z = displayMap[i+1];
            if (z < state.size())
                fout << " ";
        }
    }
}

void displayParity(ostream& fout, const vector<SliceValue>& state, const int& offset, const vector<unsigned int>& displayMap)
{
    unsigned int z=0;
    for(unsigned int i=0; i<displayMap.size(); i+=2) {
        for( ; z<displayMap[i]; z++) {
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
            if (z < (displayMap[i]-1))
                fout << "   ";
            else if (z < (state.size()-1))
                fout << " ";
        }
        if (displayMap[i] < displayMap[i+1]) {
            unsigned int numberOfSpaces = getDisplayNumberOfSpaces(displayMap[i+1] - displayMap[i]);
            for(unsigned int j=0; j<numberOfSpaces; j++)
                fout << " ";
            z = displayMap[i+1];
            if (z < state.size())
                fout << " ";
        }
    }
}

void displayState(ostream& fout, const vector<SliceValue>& state, bool showParity)
{
    const int offset = 2;
    vector<unsigned int> displayMap;
    getDisplayMap(state, displayMap);
    for(unsigned int sy=0; sy<5; sy++) {
        unsigned int y = KeccakF::index(-1-sy-offset);
        displayPlane(fout, state, offset, y, displayMap);
        fout << endl;
    }
    if (showParity) {
        displayParity(fout, state, offset, displayMap);
        fout << endl;
    }
}

void displayStates(ostream& fout,
                   const vector<SliceValue>& state1, bool showParity1,
                   const vector<SliceValue>& state2, bool showParity2)
{
    const int offset = 2;
    vector<unsigned int> displayMap1, displayMap2;
    getDisplayMap(state1, displayMap1);
    getDisplayMap(state2, displayMap2);
    for(unsigned int sy=0; sy<5; sy++) {
        unsigned int y = KeccakF::index(-1-sy-offset);
        displayPlane(fout, state1, offset, y, displayMap1);
        fout << "  |  "; 
        displayPlane(fout, state2, offset, y, displayMap2);
        fout << endl;
    }
    if (showParity1 || showParity2) {
        if (showParity1)
            displayParity(fout, state1, offset, displayMap1);
        else
            displayNothing(fout, state1, displayMap1);
        if (showParity2) {
            fout << "     ";
            displayParity(fout, state2, offset, displayMap2);
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
    vector<unsigned int> displayMap1, displayMap2, displayMap3;
    getDisplayMap(state1, displayMap1);
    getDisplayMap(state2, displayMap2);
    getDisplayMap(state3, displayMap3);
    for(unsigned int sy=0; sy<5; sy++) {
        unsigned int y = KeccakF::index(-1-sy-offset);
        displayPlane(fout, state1, offset, y, displayMap1);
        fout << "  |  "; 
        displayPlane(fout, state2, offset, y, displayMap2);
        fout << "  |  ";
        displayPlane(fout, state3, offset, y, displayMap3);
        fout << endl;
    }
    if (showParity1 || showParity2 || showParity3) {
        if (showParity1)
            displayParity(fout, state1, offset, displayMap1);
        else
            displayNothing(fout, state1, displayMap1);
        fout << "     ";
        if (showParity2)
            displayParity(fout, state2, offset, displayMap2);
        else
            displayNothing(fout, state2, displayMap2);
        if (showParity3) {
            fout << "     ";
            displayParity(fout, state3, offset, displayMap3);
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

RowValue getRow(const vector<LaneValue>& lanes, unsigned int y, unsigned int z)
{
    RowValue result = 0;
    for(unsigned int x=0; x<5; x++) {
        if ((lanes[KeccakF::index(x, y)] & ((LaneValue)1 << z)) != 0)
            result ^= (1 << x);
    }
    return result;
}

void setRow(vector<LaneValue>& lanes, RowValue row, unsigned int y, unsigned int z)
{
    for(unsigned int x=0; x<5; x++)
        if ((row & (1 << x)) != 0)
            lanes[KeccakF::index(x, y)] |= (LaneValue)1 << z;
        else
            lanes[KeccakF::index(x, y)] &= ~((LaneValue)1 << z);
}

RowValue getRow(const vector<SliceValue>& slices, unsigned int y, unsigned int z)
{
    return getRowFromSlice(slices[z], y);
}

void setRow(vector<SliceValue>& slices, RowValue row, unsigned int y, unsigned int z)
{
    slices[z] = (slices[z] & (~getSliceFromRow(0x1F, y))) | getSliceFromRow(row, y);
}

ColumnValue getColumn(const vector<SliceValue>& slices, unsigned int x, unsigned int z)
{
    ColumnValue column = 0;
    for(unsigned int y=0; y<5; y++)
        column |= ((slices[z] >> KeccakF::index(x, y)) & 1) << y;
    return column;
}

void setColumn(vector<SliceValue>& slices, ColumnValue column, unsigned int x, unsigned int z)
{
    const SliceValue columnMask = (1 << 0) | (1 << 5) | (1 << 10) | (1 << 15) | (1 << 20);
    slices[z] &= ~(columnMask << x);
    for(unsigned int y=0; y<5; y++)
        slices[z] |= ((column >> y) & 1) << KeccakF::index(x, y);
}

void invertColumn(vector<SliceValue>& slices, unsigned int x, unsigned int z)
{
    const SliceValue columnMask = (1 << 0) | (1 << 5) | (1 << 10) | (1 << 15) | (1 << 20);
    slices[z] ^= (columnMask << x);
}

SliceValue getSlice(const vector<LaneValue>& lanes, unsigned int z)
{
    SliceValue result = 0;
    for(unsigned int y=0; y<5; y++)
        result ^= getSliceFromRow(getRow(lanes, y, z), y);
    return result;
}

void setSlice(vector<LaneValue>& lanes, SliceValue slice, unsigned int z)
{
    for(unsigned int y=0; y<5; y++)
        setRow(lanes, getRowFromSlice(slice, y), y, z);
}

void fromLanesToSlices(const vector<LaneValue>& lanes, vector<SliceValue>& slices, unsigned int laneSize)
{
    slices.resize(laneSize);
    for(unsigned int z=0; z<laneSize; z++)
        slices[z] = getSlice(lanes, z);
}

void fromSlicesToLanes(const vector<SliceValue>& slices, vector<LaneValue>& lanes)
{
    lanes.assign(25, 0);
    for(unsigned int z=0; z<slices.size(); z++)
        setSlice(lanes, slices[z], z);
}
