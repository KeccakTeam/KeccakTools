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

#ifndef _KECCAKFPARTS_H_
#define _KECCAKFPARTS_H_

#include "Keccak-f.h"
#include "Keccak-fPositions.h"

using namespace std;

// -------------------------------------------------------------
//
// Rows, columns, lanes and slices
//
// -------------------------------------------------------------

/** The number of rows and columns in Keccak-<i>f</i>.
  */
const int nrRowsAndColumns = 5;

/** The RowValue type is one byte, containing the 5
  * bits of a row, in the least significant bits of the byte.
  */
typedef unsigned char RowValue;

/** The ColumnsValue type is one byte, containing the 5
  * bits of a column, in the least significant bits of the byte.
  */
typedef unsigned char ColumnValue;

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

/** The LaneIndex type codes the x and y coordinates as the single integer x + 5y.
  */
typedef unsigned int LaneIndex;

/** This method returns the lane index from (x,y) coordinates.
  *
  * @param  x   The x coordinate, with 0 ≤ x < 5.
  * @param  y   The y coordinate, with 0 ≤ y < 5.
  */
inline LaneIndex getLaneIndex(unsigned int x, unsigned int y)
{
    return x + 5*y;
}

/** This method returns the lane index from (x,y) coordinates.
  *
  * @param  x   The x coordinate, which can be any signed integer (reduced modulo 5).
  * @param  y   The y coordinate, which can be any signed integer (reduced modulo 5).
  */
LaneIndex getLaneIndexSafely(int x, int y);

/** This method returns the value of a given bit in a state.
  *
  * @param  slices   The state as a vector of slices.
  * @param  x   The x coordinate.
  * @param  y   The y coordinate.
  * @param  z   The z coordinate.
  */
inline int getBit(const vector<SliceValue>& slices, unsigned int x, unsigned int y, unsigned int z)
{
    return (slices[z] >> (x+5*y)) & 1;
}

/** This method returns the value of a given bit in a state.
  *
  * @param  slices   The state as a vector of slices.
  * @param  p   The x, y, z coordinates.
  */
inline int getBit(const vector<SliceValue>& slices, const BitPosition& p)
{
    return getBit(slices, p.x, p.y, p.z);
}

/** This method sets to 0 a particular bit in a state.
  *
  * @param  slices   The state as a vector of slices.
  * @param  x   The x coordinate.
  * @param  y   The y coordinate.
  * @param  z   The z coordinate.
  */
inline void setBitToZero(vector<SliceValue>& slices, unsigned int x, unsigned int y, unsigned int z)
{
    slices[z] &= ~((SliceValue)1 << (x+5*y));
}

/** This method sets to 0 a particular bit in a state.
  *
  * @param  slices   The state as a vector of slices.
  * @param  p   The x, y, z coordinates.
  */
inline void setBitToZero(vector<SliceValue>& slices, const BitPosition& p)
{
    return setBitToZero(slices, p.x, p.y, p.z);
}

/** This method sets to 1 a particular bit in a state.
  *
  * @param  slices   The state as a vector of slices.
  * @param  x   The x coordinate.
  * @param  y   The y coordinate.
  * @param  z   The z coordinate.
  */
inline void setBitToOne(vector<SliceValue>& slices, unsigned int x, unsigned int y, unsigned int z)
{
    slices[z] |= (SliceValue)1 << (x+5*y);
}

/** This method sets to 1 a particular bit in a state.
  *
  * @param  slices   The state as a vector of slices.
  * @param  p   The x, y, z coordinates.
  */
inline void setBitToOne(vector<SliceValue>& slices, const BitPosition& p)
{
    setBitToOne(slices, p.x, p.y, p.z);
}

/** This method inverts a particular bit in a state.
  *
  * @param  slices   The state as a vector of slices.
  * @param  x   The x coordinate.
  * @param  y   The y coordinate.
  * @param  z   The z coordinate.
  */
inline void invertBit(vector<SliceValue>& slices, unsigned int x, unsigned int y, unsigned int z)
{
    slices[z] ^= (SliceValue)1 << (x+5*y);
}

/** This method inverts a particular bit in a state.
  *
  * @param  slices   The state as a vector of slices.
  * @param  p   The x, y, z coordinates.
  */
inline void invertBit(vector<SliceValue>& slices, const BitPosition& p)
{
    invertBit(slices, p.x, p.y, p.z);
}

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

/** This method returns the value of a given bit in a plane.
  *
  * @param  slices   The state as a vector of rows.
  * @param  x   The x coordinate.
  * @param  z   The z coordinate.
  */
inline int getBit(const vector<RowValue>& rows, unsigned int x, unsigned int z)
{
    return (rows[z] >> x) & 1;
}

/** This method sets to 0 a particular bit in a plane.
  *
  * @param  slices   The state as a vector of rows.
  * @param  x   The x coordinate.
  * @param  z   The z coordinate.
  */
inline void setBitToZero(vector<RowValue>& rows, unsigned int x, unsigned int z)
{
    rows[z] &= ~((RowValue)1 << x);
}

/** This method sets to 1 a particular bit in a plane.
  *
  * @param  slices   The state as a vector of rows.
  * @param  x   The x coordinate.
  * @param  z   The z coordinate.
  */
inline void setBitToOne(vector<RowValue>& rows, unsigned int x, unsigned int z)
{
    rows[z] |= (RowValue)1 << x;
}

/** This method returns the value of a given row in a slice.
  *
  * @param  slices   The state as a vector of slices.
  * @param  y   The y coordinate.
  * @param  z   The z coordinate.
  */
RowValue getRow(const vector<SliceValue>& slices, unsigned int y = 0, unsigned int z = 0);

/** This method returns the value of a given row in a slice.
  *
  * @param  slices   The state as a vector of slices.
  * @param  p   The y, z coordinates.
  */
inline RowValue getRow(const vector<SliceValue>& slices, const RowPosition& p)
{
    return getRow(slices, p.y, p.z);
}

/** This method sets the value of a particular row in a vector of slices.
  *
  * @param  slices   The state as a vector of slices.
  * @param  row     The row value.
  * @param  y   The y coordinate.
  * @param  z   The z coordinate.
  */
void setRow(vector<SliceValue>& slices, RowValue row, unsigned int y = 0, unsigned int z = 0);

/** This method sets the value of a particular row in a vector of slices.
  *
  * @param  slices   The state as a vector of slices.
  * @param  row     The row value.
  * @param  p   The y, z coordinates.
  */
inline void setRow(vector<SliceValue>& slices, RowValue row, const RowPosition& p)
{
    setRow(slices, row, p.y, p.z);
}

/** This method constructs a slice value from 5 row values.
  */
SliceValue getSliceValue(RowValue row0, RowValue row1, RowValue row2, RowValue row3, RowValue row4);

/** This method returns the value of a given column in a slice.
  *
  * @param  slices   The state as a vector of slices.
  * @param  x   The x coordinate.
  * @param  z   The z coordinate.
  */
ColumnValue getColumn(const vector<SliceValue>& slices, unsigned int x = 0, unsigned int z = 0);

/** This method sets the value of a particular column in a vector of slices.
  *
  * @param  slices   The state as a vector of slices.
  * @param  column   The row value.
  * @param  x   The x coordinate.
  * @param  z   The z coordinate.
  */
void setColumn(vector<SliceValue>& slices, ColumnValue column, unsigned int x = 0, unsigned int z = 0);

/** This method complements all the bits of a particular column in a vector of slices.
  *
  * @param  slices   The state as a vector of slices.
  * @param  x   The x coordinate.
  * @param  z   The z coordinate.
  */
void invertColumn(vector<SliceValue>& slices, unsigned int x = 0, unsigned int z = 0);

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

/** This function translates a slice value along the X and Y axes and returns the
  * translated value. Note that 0 <= @a dx < 5 and 0 <= @a dy < 5 are required.
  */
SliceValue translateSlice(SliceValue slice, unsigned int dx, unsigned int dy);

/** Same as translateSlice(), but any (negative and positive) value of dx and dy is
  * allowed.
  */
SliceValue translateSliceSafely(SliceValue slice, int dx, int dy);

/**
  * This function translates the state along the Z axis.
  */
void translateStateAlongZ(vector<SliceValue>& state, unsigned int dz);

/** This method returns the value of a given row in a slice.
  *
  * @param  lanes   The state as a vector of lanes.
  * @param  y   The y coordinate.
  * @param  z   The z coordinate.
  */
RowValue getRow(const vector<LaneValue>& lanes, unsigned int y = 0, unsigned int z = 0);

/** This method sets the value of a particular row in a vector of lanes.
  *
  * @param  lanes   The state as a vector of lanes.
  * @param  row     The row value.
  * @param  y   The y coordinate.
  * @param  z   The z coordinate.
  */
void setRow(vector<LaneValue>& lanes, RowValue row, unsigned int y = 0, unsigned int z = 0);

/** This method returns the value of a given slice in a state represented as
  * a vector of lanes.
  *
  * @param  lanes   The state as a vector of lanes.
  * @param  z   The slice index (z coordinate).
  */
SliceValue getSlice(const vector<LaneValue>& lanes, unsigned int z = 0);

/** This method sets the value of a particular slice in a vector of lanes.
  *
  * @param  lanes   The state as a vector of lanes.
  * @param  slice   The slice value.
  * @param  z   The z coordinate.
  */
void setSlice(vector<LaneValue>& lanes, SliceValue slice, unsigned int z = 0);

/** This method creates the value of a state represented as a vector of slices
  * from a state represented as a vector of lanes.
  *
  * @param  slices  The destination for the slices.
  * @param  lanes   The state as a vector of lanes.
  * @param  laneSize    The lane size, which is also the number of slices.
  */
void fromLanesToSlices(const vector<LaneValue>& lanes, vector<SliceValue>& slices, unsigned int laneSize);

 /** This method creates the value of a state represented as a vector of lanes
  * from a state represented as a vector of slices.
  *
  * @param  lanes   The destination for the lanes.
  * @param  slices  The state as a vector of slices.
  */
void fromSlicesToLanes(const vector<SliceValue>& slices, vector<LaneValue>& lanes);


// -------------------------------------------------------------
//
// Hamming weight and related
//
// -------------------------------------------------------------

/** This function returns the Hamming weight of the given row value.
  */
unsigned int getHammingWeightRow(RowValue row);

/** This function returns the Hamming weight of the given column value.
  */
unsigned int getHammingWeightColumn(ColumnValue column);

/** This function returns the Hamming weight of the given slice value.
  */
unsigned int getHammingWeightSlice(SliceValue slice);

/** This function returns the Hamming weight of the given state.
  */
unsigned int getHammingWeight(const vector<SliceValue>& state);

/** This function returns the Hamming weight of the given lane.
  */
unsigned int getHammingWeightLane(LaneValue lane);

/** This function returns the Hamming weight of the given state.
  */
unsigned int getHammingWeight(const vector<LaneValue>& state);

/** This method returns the number of active rows in the given slice value.
  */
unsigned int getNrActiveRows(const SliceValue& slice);

/** This method returns the number of active rows in the state given as slices.
  */
unsigned int getNrActiveRows(const vector<SliceValue>& slices);

/** This method returns the number of active rows in the state given as lanes.
  */
unsigned int getNrActiveRows(const vector<LaneValue>& lanes);


#endif
