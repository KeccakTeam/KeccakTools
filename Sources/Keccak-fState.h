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

#ifndef _KECCAKFSTATE_H_
#define _KECCAKFSTATE_H_

#include "Keccak-f.h"
#include "Keccak-fParts.h"
#include "Keccak-fPositions.h"
#include <map>

using namespace std;

inline SliceValue getSlicePoint(unsigned int x, unsigned int y)
{
    return (SliceValue)1 << (x+5*y);
}

class StateAsSlices {
public:
    vector<SliceValue> slices;
public:
    StateAsSlices() {}

    /** This method returns the value of a given slice in a state.
      *
      * @param  z   The z coordinate.
      * @return The slice value as a constant reference.
      */
    inline const SliceValue& operator[](unsigned int z) const { return slices[z]; }
    
    /** This method returns the value of a given slice in a state.
      *
      * @param  z   The z coordinate.
      * @return The slice value as a non-constant reference.
      */
    inline SliceValue& operator[](unsigned int z) { return slices[z]; }
    
    /** This method returns the value of a given bit in a state.
      *
      * @param  x   The x coordinate.
      * @param  y   The y coordinate.
      * @param  z   The z coordinate.
      * @return The bit value.
      */
    inline int getBit(unsigned int x, unsigned int y, unsigned int z) const
    {
        return (slices[z] >> (x+5*y)) & 1;
    }

    /** This method returns the value of a given bit in a state.
      *
      * @param  p   The (x,y,z) coordinates.
      * @return The bit value.
      */
    inline int getBit(const BitPosition& p) const
    {
        return getBit(p.x, p.y, p.z);
    }

    /** This method sets to 0 a particular bit in a state.
      *
      * @param  x   The x coordinate.
      * @param  y   The y coordinate.
      * @param  z   The z coordinate.
      */
    inline void setBitToZero(unsigned int x, unsigned int y, unsigned int z)
    {
        slices[z] &= ~getSlicePoint(x, y);
    }

    /** This method sets to 0 a particular bit in a state.
      *
      * @param  p   The (x,y,z) coordinates.
      */
    inline void setBitToZero(const BitPosition& p)
    {
        setBitToZero(p.x, p.y, p.z);
    }

    /** This method sets to 1 a particular bit in a state.
      *
      * @param  x   The x coordinate.
      * @param  y   The y coordinate.
      * @param  z   The z coordinate.
      */
    inline void setBitToOne(unsigned int x, unsigned int y, unsigned int z)
    {
        slices[z] |= getSlicePoint(x, y);
    }

    /** This method sets to 1 a particular bit in a state.
      *
      * @param  p   The (x,y,z) coordinates.
      */
    inline void setBitToOne(const BitPosition& p)
    {
        setBitToOne(p.x, p.y, p.z);
    }

    /** This method inverts a particular bit in a state.
      *
      * @param  x   The x coordinate.
      * @param  y   The y coordinate.
      * @param  z   The z coordinate.
      */
    inline void invertBit(unsigned int x, unsigned int y, unsigned int z)
    {
        slices[z] ^= getSlicePoint(x, y);
    }

    /** This method inverts a particular bit in a state.
      *
      * @param  p   The (x,y,z) coordinates.
      */
    inline void invertBit(const BitPosition& p)
    {
        invertBit(p.x, p.y, p.z);
    }
    
};

/** The SparseState type codes a state as map of z-coordinate, SliceValue couples.
  * It is compact for sparse states as it only stores the nonzero slices.
  * The convention used and maintained by the getSlice(), setSlice(), getBit(), 
  * setBitToZero() and invertBit() functions is that a slice with value zero 
  * does not appear in the map.
  */
class SparseStateAsSlices {
public:
    map<unsigned int, SliceValue> slices;
public:
    SparseStateAsSlices() {}

    /** This method returns the value of a given bit in a state.
      *
      * @param  x   The x coordinate.
      * @param  y   The y coordinate.
      * @param  z   The z coordinate.
      */
    inline int getBit(unsigned int x, unsigned int y, unsigned int z) const
    {
        map<unsigned int, SliceValue>::const_iterator slice = slices.find(z);
        if (slice == slices.end())
            return 0;
        else
            return ((slice->second) >> (x+5*y)) & 1;
    }

    /** This method returns the value of a given bit in a state.
      *
      * @param  p   The (x,y,z) coordinates.
      */
    inline int getBit(const BitPosition& p) const
    {
        return getBit(p.x, p.y, p.z);
    }

    /** This method sets to 0 a particular bit in a state.
      *
      * @param  x   The x coordinate.
      * @param  y   The y coordinate.
      * @param  z   The z coordinate.
      */
    inline void setBitToZero(unsigned int x, unsigned int y, unsigned int z)
    {
        map<unsigned int, SliceValue>::iterator slice = slices.find(z);
        if (slice != slices.end()) {
            slice->second &= ~getSlicePoint(x, y);
            if (slice->second == 0)
                slices.erase(slice);
        }
    }

    /** This method sets to 0 a particular bit in a state.
      *
      * @param  p   The (x,y,z) coordinates.
      */
    inline void setBitToZero(const BitPosition& p)
    {
        setBitToZero(p.x, p.y, p.z);
    }
    
    /** This method sets to 1 a particular bit in a state.
      *
      * @param  x   The x coordinate.
      * @param  y   The y coordinate.
      * @param  z   The z coordinate.
      */
    inline void setBitToOne(unsigned int x, unsigned int y, unsigned int z)
    {
        map<unsigned int, SliceValue>::iterator slice = slices.find(z);
        if (slice != slices.end())
            slice->second |= getSlicePoint(x, y);
        else
            slices[z] = getSlicePoint(x, y);
    }

    /** This method sets to 1 a particular bit in a state.
      *
      * @param  p   The (x,y,z) coordinates.
      */
    inline void setBitToOne(const BitPosition& p)
    {
        setBitToOne(p.x, p.y, p.z);
    }

    /** This method inverts a particular bit in a state.
      *
      * @param  x   The x coordinate.
      * @param  y   The y coordinate.
      * @param  z   The z coordinate.
      */
    inline void invertBit(unsigned int x, unsigned int y, unsigned int z)
    {
        map<unsigned int, SliceValue>::iterator slice = slices.find(z);
        if (slice != slices.end()) {
            slice->second ^= getSlicePoint(x, y);
            if (slice->second == 0)
                slices.erase(slice);
        }
        else
            slices[z] = getSlicePoint(x, y);
    }

    /** This method inverts a particular bit in a state.
      *
      * @param  p   The (x,y,z) coordinates.
      */
    inline void invertBit(const BitPosition& p)
    {
        invertBit(p.x, p.y, p.z);
    }

    /** This method returns the value of a given slice in a state.
      *
      * @param  z   The z coordinate.
      */
    inline SliceValue getSlice(unsigned int z) const
    {
        map<unsigned int, SliceValue>::const_iterator slice = slices.find(z);
        if (slice == slices.end())
            return 0;
        else
            return slice->second;
    }

    /** This method sets the value of a given slice in a state.
      *
      * @param  z   The z coordinate.
      * @param  value   The new slice value.
      */
    inline void setSlice(unsigned int z, SliceValue value)
    {
        map<unsigned int, SliceValue>::iterator slice = slices.find(z);
        if (slice == slices.end()) {
            if (value != 0)
                slices[z] = value;
        }
        else {
            if (value == 0)
                slices.erase(slice);
            else
                slice->second = value;
        }
    }

};

#endif
