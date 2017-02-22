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
      * @param  startRoundIndex The index of the first round to perform.
      * @param  nrRounds    The number of rounds to perform.
      * @see KeccakFEquations::genAbsoluteValuesBeforeChi()
      */
    void buildDCTrailFromPair(const vector<SliceValue>& a1, const vector<SliceValue>& a2, Trail& trail, int startRoundIndex, unsigned nrRounds) const;
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
    /** This method checks whether a given pair follows a given trail.
      * The state value is given before θ, as in the normal order of
      * the round function.
      * The other state value of the pair is determined by the given trail's
      * first round difference.
      * @param   a1     The first state of the pair, as a vector of slices.
      * @param   givenTrail     The trail that the given pair should follow.
      * @param   actualTrail    The output trail that the given pair follows.
      * @param  startRoundIndex The index of the first round to perform.
      * @return A Boolean telling whether the pair follows the given trail.
      */
    bool checkPairGivenDCTrail(const vector<SliceValue>& a1, const Trail& givenTrail, Trail& actualTrail, int startRoundIndex) const;
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
