/*
Tools for the Keccak sponge function family.
Authors: Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche

This code is hereby put in the public domain. It is given as is,
without any guarantee.

For more information, feedback or questions, please refer to our website:
http://keccak.noekeon.org/
*/

#ifndef _KECCAKFDCEQUATIONS_H_
#define _KECCAKFDCEQUATIONS_H_

#include "Keccak-fEquations.h"
#include "Keccak-fPropagation.h"

/** This class is an extension of KeccakFDCLC with additional functionality
  * to display equations related to a differential trail.
  */
class KeccakFDCEquations : public KeccakFDCLC
{
public:
    /** This constructor initializes the Keccak-<i>f</i> instance
      * with the given width, see KeccakF::KeccakF().
      */
    KeccakFDCEquations(unsigned int aWidth);
    /** This method produces a differential trail from a pair of inputs.
      * The state values are given before θ, as in the normal order of
      * the round function.
      * @param   a1     The first state of the pair, as a vector of slices.
      * @param   a2     The second state of the pair, as a vector of slices.
      * @param   trail  The output trail.
      * @see KeccakFEquations::genAbsoluteValuesBeforeChi()
      */
    void buildDCTrailFromPair(vector<SliceValue>& a1, vector<SliceValue>& a2, Trail& trail) const;
    /** This method displays the equations from round to round
      * that a pair has to satisfy to follow the given trail.
      * @note Note that for row patterns 11111, 5 equations are generated, 
      *       although the fifth is redundant with the 4 first.
      * @param   fout       The stream to display to.
      * @param   trail      The trail to follow.
      * @param   forSage    If true, the equations are displayed in a format suitable
      *                     for SAGE. In this case, the equations form a comma-separated
      *                     list of polynomials that must be equal to zero.
      *                     Otherwise, the equations are displayed in the form X=f(Y).
      */
    void genDCEquations(ostream& fout, const Trail& trail, bool forSage=false) const;
protected:
    /** This method creates the list of equations that the input of χ
      * must satisfy for the given input difference to propagate to the given
      * output difference.
      */
    void getDCEquations(const vector<SliceValue>& diffIn, const vector<SliceValue>& diffOut,
        const vector<SymbolicLane>& input, vector<SymbolicBit>& inputRelations) const;
    /** This method creates the list of equations that the input of χ (for one row)
      * must satisfy for the given input difference to propagate to the given
      * output difference.
      */
    void getDCEquations(RowValue diffIn, RowValue diffOut,
        const vector<SymbolicBit>& inputVariables, vector<SymbolicBit>& inputRelations) const;
    /** This method produces the display of the equations. */
    void displayEquations(ostream& fout, const vector<SymbolicLane>& state, const string& prefixOutput, bool forSage=false) const;
};

#endif
