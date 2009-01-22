/*
Tools for the Keccak sponge function family.
Authors: Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche

This code is hereby put in the public domain. It is given as is, 
without any guarantee.

For more information, feedback or questions, please refer to our website:
http://keccak.noekeon.org/
*/

#ifndef _KECCAKFPARTS_H_
#define _KECCAKFPARTS_H_

#include "Keccak-f.h"

using namespace std;

/** The RowValue type is one byte, containing the 5
  * bits of a row, in the least significant bits of the byte.
  */
typedef unsigned char RowValue;

/** The SliceValue type is one 32-bit word, containing the 5
  * rows of a slice, each located in 5 bits of this word. The row y is
  * in the bits corresponding to numerical value (0-31)*32^y 
  * in the word.
  * See getSliceFromRow() and getRowFromSlice() for more details.
  */
typedef UINT32 SliceValue;

/** This constant indicates the maximum value if one needs to loop through
  * all the possible slice values.
  */
const SliceValue maxSliceValue = 0x1FFFFFF;

/** This function returns a SliceValue with bits set to zero, except at 
  * row y, where the value is given by the argument @a row.
  */
inline SliceValue getSliceFromRow(const RowValue& row, const unsigned int& y)
{
    return (SliceValue)row << (5*y);
}

/** This function returns the row value at row y in the given slice value.
  */
inline RowValue getRowFromSlice(const SliceValue& slice, const unsigned int& y)
{
    return (slice >> (5*y)) & 0x1F;
}

/** This method constructs a slice value from 5 row values.
  */
SliceValue getSliceValue(RowValue row0, RowValue row1, RowValue row2, RowValue row3, RowValue row4);

/** This method returns the number of active rows in the given slice value.
  */
unsigned int getNrActiveRows(const SliceValue& slice);

/** This function translates a row value along the X axis and returns the
  * translated value. Note that 0 <= @a dx < 5 is required.
  */
inline RowValue translateRow(const RowValue& row, const unsigned int& dx)
{
    if (dx == 0) 
        return row;
    else
        return ((row << dx) | (row >> (5-dx))) & 0x1F;
}

/** Same as translateRow, but any (negative and positive) value of dx is allowed.
  */
RowValue translateRowSafely(RowValue row, int dx);

/** This function returns the Hamming weight of the given row value.
  */
unsigned int getHammingWeightRow(RowValue row);

/** This function translates a slice value along the X and Y axes and returns the
  * translated value. Note that 0 <= @a dx < 5 and 0 <= @a dy < 5 are required.
  */
SliceValue translateSlice(SliceValue slice, unsigned int dx, unsigned int dy);

/** Same as translateSlice(), but any (negative and positive) value of dx and dy is
  * allowed.
  */
SliceValue translateSliceSafely(SliceValue slice, int dx, int dy);

/** This function returns the Hamming weight of the given slice value.
  */
unsigned int getHammingWeightSlice(SliceValue slice);

/** This method outputs to @a fout the value of the slice in a human readable way.
  */
void displaySlice(ostream& fout, SliceValue slice);

/**
  * This function translates the state along the Z axis.
  */
void translateStateAlongZ(vector<SliceValue>& state, unsigned int dz);

/** This function returns the Hamming weight of the given state.
  */
unsigned int getHammingWeight(const vector<SliceValue>& state);

/** This method outputs to fout the value of the state in a human readable way.
  * The slices are displayed side by side.
  */
void displayState(ostream& fout, const vector<SliceValue>& state, bool showParity = false);

/** This method outputs to fout the value of the 2 states in a human readable way.
  */
void displayStates(ostream& fout,
                   const vector<SliceValue>& state1, bool showParity1,
                   const vector<SliceValue>& state2, bool showParity2);

/** This method outputs to fout the value of the 3 states in a human readable way.
  */
void displayStates(ostream& fout,
                   const vector<SliceValue>& state1, bool showParity1,
                   const vector<SliceValue>& state2, bool showParity2,
                   const vector<SliceValue>& state3, bool showParity3);

/** This method returns the value of a given row in a slice.
  *
  * @param  lanes   The state as a vector of lanes.
  * @param  y   The y coordinate.
  * @param  z   The z coordinate.
  */
RowValue getRow(const vector<UINT64>& lanes, unsigned int y, unsigned int z);

/** This method sets the value of a particular row in a vector of lanes.
  *
  * @param  lanes   The state as a vector of lanes.
  * @param  row     The row value.
  * @param  y   The y coordinate.
  * @param  z   The z coordinate.
  */
void setRow(vector<UINT64>& lanes, RowValue row, unsigned int y, unsigned int z);

/** This method returns the value of a given slice in a state represented as
  * a vector of lanes.
  *
  * @param  lanes   The state as a vector of lanes.
  * @param  z   The slice index (z coordinate).
  */
SliceValue getSlice(const vector<UINT64>& lanes, unsigned int z=0);

/** This method sets the value of a particular slice in a vector of lanes.
  *
  * @param  lanes   The state as a vector of lanes.
  * @param  slice   The slice value.
  * @param  z   The z coordinate.
  */
void setSlice(vector<UINT64>& lanes, SliceValue slice, unsigned int z=0);

/** This method creates the value of a state represented as a vector of slices 
  * from a state represented as a vector of lanes.
  *
  * @param  slices  The destination for the slices.
  * @param  lanes   The state as a vector of lanes.
  * @param  laneSize    The lane size, which is also the number of slices.
  */
void fromLanesToSlices(vector<SliceValue>& slices, const vector<UINT64>& lanes, unsigned int laneSize);

 /** This method creates the value of a state represented as a vector of lanes
  * from a state represented as a vector of slices.
  *
  * @param  lanes   The destination for the lanes.
  * @param  slices  The state as a vector of slices.
  */
void fromSlicesToLanes(vector<UINT64>& lanes, const vector<SliceValue>& slices);

#endif
