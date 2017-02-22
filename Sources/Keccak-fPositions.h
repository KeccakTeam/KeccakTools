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

#ifndef _KECCAKFPOSITIONS_H_
#define _KECCAKFPOSITIONS_H_

#include <iostream>

using namespace std;

/** Class containing the x, y, z coordinates of a bit. */
class BitPosition {
public:
    /** The x-coordinate of the bit, 0 ≤ x < 5. */
    unsigned int x;
    /** The y-coordinate of the bit, 0 ≤ y < 5. */
    unsigned int y;
    /** The z-coordinate of the bit, 0 ≤ z < laneSize. */
    unsigned int z;
    /** The default constructor. */
    BitPosition() : x(0), y(0), z(0) {}
    /** The constructor.
      * @param  ax  The x-coordinate.
      * @param  ay  The y-coordinate.
      * @param  az  The z-coordinate.
      */
    BitPosition(unsigned int ax, unsigned int ay, unsigned int az);
    /** The copy constructor.
      * @param  aPoint  The bit position.
      */
    BitPosition(const BitPosition& aPoint);
    /** A display function, for use with the << operator.
      * @param  fout    The output stream to write to.
      * @param  point   The bit position to display.
      */
    friend ostream& operator<<(ostream& fout, const BitPosition& point);
    /** The equality operator.
      * @param  otherPoint  The bit position at the right of the operator.
      * @return True iff this object equals @a otherPoint.
      */
    bool operator==(const BitPosition& otherPoint) const;
    /** An ordering operator.
      * @param  otherPoint  The bit position at the right of the operator.
      * @return True iff this object comes before @a otherPoint.
      */
    bool operator<(const BitPosition& otherPoint) const;
    /** This method increments the coordinates while maintaining the z coordinate.
      * It first attempts to increment y. If not possible, it resets y and increments x.
      * @return It returns false if (x, y)=(4, 4).
      */
    bool nextXY();
    /** This method sets the three coordinates.
      * @param  ax  The x-coordinate.
      * @param  ay  The y-coordinate.
      * @param  az  The z-coordinate.
      */
    void set(unsigned int ax=0, unsigned int ay=0, unsigned int az=0);
    /** This method translates the position along the x axis.
      * @param  dx  The amount of translation (can be positive or negative).
      */
    void xTranslate(int dx);
    /** This method translates the position along the y axis.
      * @param  dy  The amount of translation (can be positive or negative).
      */
    void yTranslate(int dy);
    /** This method translates the position along the z axis.
      * @param  dz  The amount of translation (can be positive or negative).
      * @param  laneSize    The lane size.
      */
    void zTranslate(int dz, unsigned int aLaneSize);
};

/** Class containing the x, z coordinates of a column. */
class ColumnPosition
{
public:
    /** The x-coordinate of the column, 0 ≤ x < 5. */
    unsigned int x;
    /** The z-coordinate of the column, 0 ≤ z < laneSize. */
    unsigned int z;
public:
    /** The default constructor. */
    ColumnPosition() : x(0), z(0) {}
    /** The constructor.
      * @param  ax  The x-coordinate.
      * @param  az  The z-coordinate.
      */
    ColumnPosition(unsigned int ax, unsigned int az) : x(ax), z(az) {}
    /** A constructor taking a bit position.
      * @param  ap  The bit position.
      */
    ColumnPosition(const BitPosition& ap) : x(ap.x), z(ap.z) {}
    /** This function returns an integer between 0 and 5*laneSize.
      * @return x + 5*z
      */
    inline unsigned int getXplus5Z() const { return x+5*z; }
    /** An ordering operator, required when storing a ColumnPosition object
      * in a set or as the first member in maps.
      * @param  a   The column position at the left of the operator.
      * @param  az  The column position at the right of the operator.
      */
    friend bool operator<(const ColumnPosition& aCP, const ColumnPosition& bCP);
    /** A display function, for use with the << operator.
      * @param  fout    The output stream to write to.
      * @param  aCP The column position to display.
      */
    friend ostream& operator<<(ostream& fout, const ColumnPosition& aCP);
};

/** Class containing the y, z coordinates of a row. */
class RowPosition
{
public:
    /** The y-coordinate of the column, 0 ≤ x < 5. */
    unsigned int y;
    /** The z-coordinate of the column, 0 ≤ z < laneSize. */
    unsigned int z;
public:
    /** The default constructor. */
    RowPosition() : y(0), z(0) {}
    /** The constructor.
      * @param  ax  The y-coordinate.
      * @param  az  The z-coordinate.
      */
    RowPosition(unsigned int ay, unsigned int az) : y(ay), z(az) {}
    /** A constructor taking a bit position.
      * @param  ap  The bit position.
      */
    RowPosition(const BitPosition& ap) : y(ap.y), z(ap.z) {}
    /** This function returns an integer between 0 and 5*laneSize.
      * @return y + 5*z
      */
    inline unsigned int getYplus5Z() const { return y+5*z; }
    /** An ordering operator, required when storing a RowPosition object
      * in a set or as the first member in maps.
      * @param  aRP  The row position at the left of the operator.
      * @param  bRP  The row position at the right of the operator.
      */
    friend bool operator<(const RowPosition& aRP, const RowPosition& bRP);
    /** A display function, for use with the << operator.
      * @param  fout    The output stream to write to.
      * @param  aRP The row position to display.
      */
    friend ostream& operator<<(ostream& fout, const ColumnPosition& aCP);
};

#endif
