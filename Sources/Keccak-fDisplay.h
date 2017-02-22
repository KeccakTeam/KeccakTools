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

#ifndef _KECCAKFDISPLAY_H_
#define _KECCAKFDISPLAY_H_

#include "Keccak-f.h"
#include "Keccak-fPropagation.h"

using namespace std;

/** This method outputs to @a fout the value of the slice in a human readable way.
  */
void displaySlice(ostream& fout, SliceValue slice);

/** Class to help display state values in SVG.
 */
class KeccakDisplayInSVG
{
public:
    /** Whether the rows must be also displayed on top of individual bits. */
    bool displayRows;
    /** Whether the slices must be also displayed on top of individual bits or rows. */
    bool displaySlices;
    /** When displaying several slices, this attribute tells how many slices can be displayed together side by side. */
    int maxNumberOfHorizontalSlices;
    /** Size (in the units of the generated SVG file) of an individual bit. */
    double bitSize;
    /** The constructor. */
    int laneSize;
    KeccakDisplayInSVG() : bitSize(10.0), laneSize(8) {}
    /** Output SVG code for a state with given value.
      * Each slice is displayed and is offset by (0.4, -0.3)*bitSize units. */
    void displayState(ostream& fout, const vector<SliceValue>& state, double gOffsetX=0.0, double gOffsetY=0.0, double gScale=1.0) const;
    /** Output SVG code for a state with given value.
      * Each non-zero slice is displayed side by side,
      * up to @a maxNumberOfHorizontalSlices on the same row. */
    void displayStateSparsely(ostream& fout, const vector<SliceValue>& state, double gOffsetX=0.0, double gOffsetY=0.0, double gScale=1.0) const;
    /** Display a trail using displayStateSparsely(). */
    void displayTrail(ostream& fout, const KeccakFPropagation& DCorLC, const Trail& trail, double x=0.0, double y=0.0, double scale=1.0) const;
    /** Display the parity and parity effect of a state. */
    void displayParity(ostream& fout, const KeccakFPropagation& DCorLC, const vector<RowValue>& C, const vector<RowValue>& D, bool displayRuns=false) const;

    void displayRow(ostream& fout);
    void displayColumn(ostream& fout);
    void displayLane(ostream& fout);
    void displayPlane(ostream& fout);
    void displaySlice(ostream& fout);
    void displaySheet(ostream& fout);
    void displayState(ostream& fout);
protected:
    void displaySlice(ostream& fout, SliceValue slice, unsigned int z, double x=0.0, double y=0.0, double scale=1.0) const;
    virtual string getAdditionalBitStyles(unsigned int x, unsigned int y, unsigned int z) const;
    void getPosition(unsigned int x, unsigned int y, unsigned int z, double &px, double &py) const;
};

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

#endif
