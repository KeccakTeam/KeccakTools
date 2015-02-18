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

#ifndef _KECCAKFPROPAGATION_H_
#define _KECCAKFPROPAGATION_H_

#include <string>
#include "Keccak-fAffineBases.h"
#include "Keccak-fDCLC.h"
#include "Keccak-fParts.h"
using namespace std;

class ReverseStateIterator;

/** This class provides the necessary tools to compute the propagation of
  * either differences or linear patterns through the rounds of Keccak-<i>f</i>.
  * To provide methods that work similarly for linear (LC) and differential cryptanalysis (DC),
  * an instance of this class is specialized in either DC or LC.
  * The convention of the direction of propagation is as described in the Keccak main document:
  * - "direct" means in the direction that allows to describe the output patterns in an affine space
  *   (same direction as the rounds for DC, inverse rounds for LC);
  * - "reverse" means the opposite of "direct" 
  *   (hence same direction as the inverse rounds for DC, rounds for LC).
  *
  * In this context, the words "before" and "after" refer to the "direct" direction.
  */
class KeccakFPropagation {
public:
    /** The output row patterns:
      * - for DC, same as KeccakFDCLC::diffChi;
      * - for LC, same as KeccakFDCLC::corrInvChi.
      */
    vector<ListOfRowPatterns> directRowOutputListPerInput;
    /** The output row patterns in the reverse direction:
      * - for DC, same as KeccakFDCLC::diffInvChi;
      * - for LC, same as KeccakFDCLC::corrChi.
      */
    vector<ListOfRowPatterns> reverseRowOutputListPerInput;
    /** This attribute contains the same as directRowOutputListPerInput
      * but in the form of an affine space representation.
      */
    vector<AffineSpaceOfRows> affinePerInput;
    /** This is a link to the 'parent' KeccakFDCLC class.
      */
    const KeccakFDCLC& parent;
    /** This attribute contains the lane size (a copy of parent.laneSize).
      */
    unsigned int laneSize;
    /** This attribute contains a string to help build appropriate file names: "DC" or "LC".
      */
    const string name;
protected:
    /** The λ mode contains the lambdaMode attribute to pass to KeccakFDCLC
      * to compute the linear part in the "direct" direction.
      */
    KeccakFDCLC::LambdaMode lambdaMode;
    /** The λ mode contains the lambdaMode attribute to pass to KeccakFDCLC
      * to compute the linear part in the "reverse" direction.
      */
    KeccakFDCLC::LambdaMode reverseLambdaMode;
private:
    /** This attribute contains the propagation weight of every possible slice value.
      */
    vector<unsigned char> weightPerSlice;
    /** This attribute contains the minimum reverse weight of every possible slice value.
      */
    vector<unsigned char> minReverseWeightPerSlice;
    /** This vector tells whether a pattern x at the input of χ is compatible
      * with output y. This can be found in chiCompatibilityTable[x+32*y].
      * See also isChiCompatible().
      */
    vector<bool> chiCompatibilityTable;
public:
    /** This type allows one to specify the type of propagation: differential (DC) or linear (LC). */
    enum DCorLC { DC = 0, LC };
    /** This constructor initializes the different attributes as a function of the
      * Keccak-<i>f</i> instance referenced by @a aParent and
      * whether the instance handles DC or LC (@a aDCorLC).
      * @param   aParent    A reference to the Keccak-<i>f</i> instance as a KeccakFDCLC object.
      * @param   aDCorLC    The propagation type.
      */
    KeccakFPropagation(const KeccakFDCLC& aParent, KeccakFPropagation::DCorLC aDCorLC);
    /** This method returns the propagation type (DC or LC) handled by the instance.
      * @return The propagation type as a DCorLC value.
      */
    DCorLC getPropagationType() const;
    /** This function displays the possible patterns and their weights. 
      * @param  out The stream to display to.
      */
    void display(ostream& out) const;
    /** This method returns the propagation weight of a slice.
      * @param   slice  The value of a slice.
      * @return The propagation weight of the given slice.
      */
    inline unsigned int getWeight(const SliceValue& slice) const { return weightPerSlice[slice]; }
    /** This method returns the propagation weight of a row.
      * @param   row    The value of a slice.
      * @return The propagation weight of the given row.
      */
    inline unsigned int getWeightRow(const RowValue& row) const { return getWeight(getSliceFromRow(row, 0)); }
    /** This method returns the propagation weight of a state.
      * @param   state  The value of a state given as a vector of slices.
      * @return The propagation weight of the given state.
      */
    unsigned int getWeight(const vector<SliceValue>& state) const;
    /** This method returns the minimum reverse weight of a slice.
      * @param   slice  The value of a slice.
      * @return The minimum weight of the given slice.
      */
    inline unsigned int getMinReverseWeight(const SliceValue& slice) const { return minReverseWeightPerSlice[slice]; }
    /** This method returns the minimum reverse weight of a slice.
      * @param   slice  The value of a slice.
      * @return The minimum weight of the given slice.
      */
    inline unsigned int getMinReverseWeightRow(const RowValue& row) const { return getMinReverseWeight(getSliceFromRow(row, 0)); }
    /** This method returns the minimum reverse weight of a state.
      * @param   state  The value of a state given as a vector of slices.
      * @return The minimum reverse weight of the given state.
      */
    unsigned int getMinReverseWeight(const vector<SliceValue>& state) const;
    /** This method returns the minimum reverse weight of a state, to which
      * the reverse λ is first applied. 
      * This allows to give a state value before χ (so after λ),
      * which is then converted to a state value after the χ of the previous round
      * (so before λ).
      * @param   state  The value of a state given as a vector of slices.
      * @return The minimum reverse weight.
      */
    unsigned int getMinReverseWeightAfterLambda(const vector<SliceValue>& state) const;
    /** This method multiplies the vector (dx, dy)<sup>T</sup> by the matrix π or π<sup>-1</sup>
      * to the left:
      * - for DC, π is used;
      * - for LC, π<sup>-1</sup> is used.
      * @param  dx  The x coordinate to update.
      * @param  dy  The y coordinate to update.
      */
    void directPi(unsigned int& dx, unsigned int& dy) const;
    /** This method multiplies the vector (dx, dy)<sup>T</sup> by the matrix π or π<sup>-1</sup> 
      * to the left:
      * - for DC, π<sup>-1</sup> is used;
      * - for LC, π is used.
      * @param  dx  The x coordinate to update.
      * @param  dy  The y coordinate to update.
      */
    void reversePi(unsigned int& dx, unsigned int& dy) const;
    /** This method moves the given bit position through ρ and π in the direct direction:
      * - for DC, this computes π(ρ(x, y, z));
      * - for LC, this computes ρ<sup>-1</sup>(π<sup>-1</sup>(x, y, z)).
      * @param  point   The coordinates (x, y, z) to update.
      */
    void directRhoPi(BitPosition& point) const;
    /** This method moves the given bit position through ρ and π in the reverse direction:
      * - for DC, this computes ρ<sup>-1</sup>(π<sup>-1</sup>(x, y, z));
      * - for LC, this computes π(ρ(x, y, z)).
      * @param  point   The coordinates (x, y, z) to update.
      */
    void reverseRhoPi(BitPosition& point) const;
    /** This method moves the given bit position through the operations after θ in the direct direction:
      * - for DC, this computes π(ρ(x, y, z));
      * - for LC, this leaves the given bit position unchanged.
      * @param  point   The coordinates (x, y, z) to update.
      */
    void directRhoPiAfterTheta(BitPosition& point) const;
    /** This method moves the given bit position through the operations before θ in the reverse direction:
      * - for DC, this leaves the given bit position unchanged;
      * - for LC, this computes π(ρ(x, y, z)).
      * @param  point   The coordinates (x, y, z) to update.
      */
    void reverseRhoPiBeforeTheta(BitPosition& point) const;
    /** This method computes a lower bound on the propagation weight 
      * for any state having the given the Hamming weight.
      * The formula is given in Section 3.1 of "The Keccak reference".
      * @param  hammingWeight   The Hamming weight.
      * @return The computed lower bound.
      */
    unsigned int getLowerBoundOnWeightGivenHammingWeight(unsigned int hammingWeight) const;
    /** This method computes a lower bound on the propagation weight 
      * for any state with given lower bounds on its Hamming weight and number of active rows.
      * The formulas are as follows. Let <i>l</i> be the resulting lower bound,
      * <i>n</i> be the number of active rows and <i>h</i> be the Hamming weight.
      * - For DC:
      *     - if (<i>n</i> ≥ <i>h</i>): <i>l</i> = 2<i>n</i>;
      *     - else if (5<i>n</i> ≤ <i>h</i>): refer to getLowerBoundOnWeightGivenHammingWeight();
      *     - else <i>l</i> = ceil((3<i>n</i> + <i>h</i>)/2).
      * - For LC:
      *     - if (2<i>n</i> ≥ <i>h</i>): <i>l</i> = 2<i>n</i>;
      *     - else if (5<i>n</i> ≤ <i>h</i>): refer to getLowerBoundOnWeightGivenHammingWeight();
      *     - else <i>l</i> = 2×ceil((<i>n</i> + <i>h</i>)/3).
      * @param  hammingWeight   Lower bound on the Hamming weight.
      * @param  nrActiveRows    Lower bound on the number of active rows.
      * @return The computed lower bound.
      */
    unsigned int getLowerBoundOnWeightGivenHammingWeightAndNrActiveRows(unsigned int hammingWeight, unsigned int nrOfActiveRows) const;
    /** This method computes a lower bound on the minimum reverse weight 
      * for any state having the given Hamming weight.
      * The formula is given in Section 3.1 of "The Keccak reference".
      * @param  hammingWeight   The Hamming weight.
      * @return The computed lower bound.
      */
    unsigned int getLowerBoundOnReverseWeightGivenHammingWeight(unsigned int hammingWeight) const;
    /** This method computes a lower bound on the minimum reverse weight 
      * for any state with given lower bounds on its Hamming weight and number of active rows.
      * The formulas are as follows. Let <i>l</i> be the resulting lower bound,
      * <i>n</i> be the number of active rows and <i>h</i> be the Hamming weight.
      * - For DC:
      *     - if (3<i>n</i> ≥ <i>h</i>): <i>l</i> = 2<i>n</i>;
      *     - else if (5<i>n</i> ≤ <i>h</i>): refer to getLowerBoundOnReverseWeightGivenHammingWeight();
      *     - else <i>l</i> = ceil((<i>n</i> + <i>h</i>)/2).
      * - For LC:
      *     - if (4<i>n</i> ≥ <i>h</i>): <i>l</i> = 2<i>n</i>;
      *     - else refer to getLowerBoundOnReverseWeightGivenHammingWeight().
      * @param  hammingWeight   Lower bound on the Hamming weight.
      * @param  nrActiveRows    Lower bound on the number of active rows.
      * @return The computed lower bound.
      */
    unsigned int getLowerBoundOnReverseWeightGivenHammingWeightAndNrActiveRows(unsigned int hammingWeight, unsigned int nrActiveRows) const;
    /** This method returns a slice value that has the bits set that are active in all slice values that satisfy the two following conditions:
      * - they are compatible through χ with the slice value in parameter;
      * - they are in the kernel.
      * @note   Actually, the LC variant is not really useful here.
      * @param   sliceBeforeChi      Slice value before χ.
      */
    SliceValue getMinimumInKernelSliceAfterChi(const SliceValue& sliceBeforeChi) const;
    /** This method returns a slice value that has the bits set that are active in all slice values that satisfy the two following conditions:
      * - they are compatible through χ<sup>-1</sup> with the slice value in parameter;
      * - they are in the kernel.
      * @note   Actually, the DC variant is not really useful here.
      * @param   sliceAfterChi      Slice value before χ.
      */
    SliceValue getMinimumInKernelSliceBeforeChi(const SliceValue& sliceAfterChi) const;
    /** This method builds an affine set of slices from a given slice pattern.
      * The produced base defines the slices just after χ (so before λ).
      * The parities considered are also those of the slices after χ.
      * @param  slice   The slice before χ to propagate.
      * @return The affine space as a AffineSpaceOfSlices object.
      */
    AffineSpaceOfSlices buildSliceBase(SliceValue slice) const;
    /** This method builds an affine set of states corresponding to the propagation
      * of a given input state through χ and λ.
      * The affine space produced thus covers the propagation through a whole round.
      * The parities in the AffineSpaceOfStates object are those before θ.
      * @param   state  The state before χ to propagate, given as a vector of slices.
      * @param   packedIfPossible If true, the produced object will have AffineSpaceOfStates::packed
      *                           set to true, unless the parities do not fit in the PackedParity type.
      *                           If false, the produced object will have AffineSpaceOfStates::packed
      *                           set to false.
      * @return The affine space as a AffineSpaceOfStates object.
      */
    AffineSpaceOfStates buildStateBase(const vector<SliceValue>& state, bool packedIfPossible = false) const;
    /** This method builds an iterator over the possible states propagating through χ 
      * in the "reverse" direction. The iterator can be restricted to run through the states
      * only up to a given maximum propagation weight.
      * @param   stateAfterChi  The state just after χ given as a vector of slices.
      * @param   maxWeight      The maximum propagation weight considered by the iterator.
      *                         If 0, the iterator runs through all the possible states.
      * @return The iterator as a ReverseStateIterator object.
      */
    ReverseStateIterator getReverseStateIterator(const vector<SliceValue>& stateAfterChi, unsigned int maxWeight = 0) const;
    /** This method returns true iff the input row pattern is compatible with the output row pattern.
      * @param   beforeChi  The row value at the input of χ.
      * @param   afterChi   The row value at the output of χ.
      * @return It returns true iff the given values are compatible through χ.
      */
    inline bool isChiCompatible(const RowValue& beforeChi, const RowValue& afterChi) const
    {
        return chiCompatibilityTable[beforeChi+32*afterChi];
    }
    /** This method returns true iff the given state before χ is compatible with the given state after χ.
      * @param   beforeChi  The state value at the input of χ.
      * @param   afterChi   The state value at the output of χ.
      * @return It returns true iff the given values are compatible through χ.
      */
    bool isChiCompatible(const vector<SliceValue>& beforeChi, const vector<SliceValue>& afterChi) const;
    /** This method returns true iff two trails can be chained, i.e., if last state of the first trail
      * is compatible through χ and λ with the first state of the second trail.
      * @param   first  The first trail.
      * @param   second The second trail.
      * @return It returns true iff the given trails can be chained.
      */
    bool isRoundCompatible(const Trail& first, const Trail& second) const;
    /** This method returns true iff θ (or θ<sup>T</sup>) is the first step of the linear step between two χ's.
      * @return It returns true iff θ (or θ<sup>T</sup>) is the first step of the linear step between two χ's.
      */
    bool isThetaJustAfterChi() const;
    /** This method applies λ in the "direct" direction:
      * - DC: θ then ρ then π;
      * - LC: π<sup>-1</sup> then ρ<sup>-1</sup> then θ<sup>T</sup>.
      * @param   in     The input state value given as a vector of slices.
      * @param   out    The output state value returned as a vector of slices.
      */
    void directLambda(const vector<SliceValue>& in, vector<SliceValue>& out) const;
    /** This method applies λ in the "reverse" direction:
      * - DC: π<sup>-1</sup> then ρ<sup>-1</sup> then θ<sup>-1</sup>;
      * - LC: θ<sup>-1T</sup> then ρ then π.
      * @param   in     The input state value given as a vector of slices.
      * @param   out    The output state value returned as a vector of slices.
      */
    void reverseLambda(const vector<SliceValue>& in, vector<SliceValue>& out) const;
    /** This method applies the part of λ before θ in the "direct" direction:
      * - DC: identity;
      * - LC: π<sup>-1</sup> then ρ<sup>-1</sup>.
      * @param   in     The input state value given as a vector of slices.
      * @param   out    The output state value returned as a vector of slices.
      */
    void directLambdaBeforeTheta(const vector<SliceValue>& in, vector<SliceValue>& out) const;
    /** This method applies the part of λ before θ in the "reverse" direction:
      * - DC: identity;
      * - LC: ρ then π.
      * @param   in     The input state value given as a vector of slices.
      * @param   out    The output state value returned as a vector of slices.
      */
    void reverseLambdaBeforeTheta(const vector<SliceValue>& in, vector<SliceValue>& out) const;
    /** This method applies θ in the "direct" direction:
      * - DC: θ;
      * - LC: θ<sup>T</sup>.
      * @param   in     The input state value given as a vector of slices.
      * @param   out    The output state value returned as a vector of slices.
      */
    void directTheta(const vector<SliceValue>& in, vector<SliceValue>& out) const;
    /** This method applies θ in the "reverse" direction:
      * - DC: θ<sup>-1</sup>;
      * - LC: θ<sup>-1T</sup>.
      * @param   in     The input state value given as a vector of slices.
      * @param   out    The output state value returned as a vector of slices.
      */
    void reverseTheta(const vector<SliceValue>& in, vector<SliceValue>& out) const;
    /** This method applies the part of λ after θ in the "direct" direction:
      * - DC: ρ then π;
      * - LC: identity.
      * @param   in     The input state value given as a vector of slices.
      * @param   out    The output state value returned as a vector of slices.
      */
    void directLambdaAfterTheta(const vector<SliceValue>& in, vector<SliceValue>& out) const;
    /** This method applies the part of λ after θ in the "reverse" direction:
      * - DC: π<sup>-1</sup> then ρ<sup>-1</sup>;
      * - LC: identity.
      * @param   in     The input state value given as a vector of slices.
      * @param   out    The output state value returned as a vector of slices.
      */
    void reverseLambdaAfterTheta(const vector<SliceValue>& in, vector<SliceValue>& out) const;
    /** This function computes the θ-effect from the parity, in the "direct" direction:
      * - DC: the θ-effect;
      * - LC: the θ<sup>T</sup>-effect.
      * @param  C   The parity as a vector of lanes.
      * @param  D   The resulting θ-effect.
      */
    void directThetaEffectFromParities(const vector<LaneValue>& C, vector<LaneValue>& D) const;
    /** This function computes the θ-effect from the parity, in the "direct" direction:
      * - DC: the θ-effect;
      * - LC: the θ<sup>T</sup>-effect.
      * @param  C   The parity as a vector of rows.
      * @param  D   The resulting θ-effect.
      */
    void directThetaEffectFromParities(const vector<RowValue>& C, vector<RowValue>& D) const;
    /** This function converts the t coordinate into (x,z) coordinates.
      * - DC: (x,z)=(-2t,t)
      * - LC: (x,z)=(2t,-t)
      * Note that the t coordinate is considered modulo 5w, so the input t can be ≥ 5w.
      * @param  t   The t coordinate.
      * @param  x   The resulting x coordinate.
      * @param  z   The resulting z coordinate.
      */
    void getXandZfromT(unsigned int t, unsigned int& x, unsigned int& z) const;
    /** This function translates a point expressed in the t coordinate along the x axis.
      * - DC: x goes to x+1 (e.g., t←t+192 if lane size is 64)
      * - LC: x goes to x-1 (e.g., t←t+192 if lane size is 64)
      * @param  t   The t coordinate of the point to translate.
      * @return The t coordinate after translation.
      */
    unsigned int translateAlongXinT(unsigned int t) const;
    /** This methods reads all the trails in a file, checks their consistency
      * and then produces a report in the given output stream.
      * See also Trail::produceHumanReadableFile().
      * @param   fileNameIn The name of the file containing the trails.
      * @param   fout       The output stream to send the report to.
      * @param   maxWeight  The maximum weight to display trails. 
      *                     If 0, the maximum weight of trails to display is
      *                     computed automatically so that a reasonable number
      *                     of trails are displayed in the report.
      * @return The number of trails read and checked.
      */
    UINT64 displayTrailsAndCheck(const string& fileNameIn, ostream& fout, unsigned int maxWeight = 0) const;
    /** Displays the parity pattern and its effect on an ostream.
     *  @param  fout    The stream to display to.
     *  @param C    The parity as a vector of rows.
     */
    void displayParity(ostream& fout, const vector<RowValue>& C) const;
    /** Displays the parity pattern and its effect on an ostream.
     *  @param  fout    The stream to display to.
     *  @param p    The parity, packed.
     */
    void displayParity(ostream& fout, PackedParity p) const;
    /** This method builds a file name by prepending "DC" or "LC" as a prefix 
      * and appending a given suffix to the name produced by 
      * KeccakFDCLC::getName().
      * @param   suffix The given suffix.
      * @return The constructed file name.
      */
    string buildFileName(const string& suffix) const;
    /** This method builds a file name by prepending "DC" or "LC" and a given prefix 
      * and appending a given suffix to the name produced by 
      * KeccakFDCLC::getName().
      * @param   prefix The given prefix.
      * @param   suffix The given suffix.
      * @return The constructed file name.
      */
    string buildFileName(const string& prefix, const string& suffix) const;
    /** This method transforms a trail core into a trail or a trail prefix by
      * choosing a first state arbitrarily.
      * The choice is not entirely arbitrary in the sense that the weight of
      * the chosen first state matches the minimum reverse weight of the state
      * it propagates to.
      * @param  trail   The trail to transform.
      */
    void specifyFirstStateArbitrarily(Trail& trail) const;
    /** This method transforms a trail prefix into a fully specified trail by
      * choosing a state after the last χ arbitrarily.
      * @param  trail   The trail to transform.
      */
    void specifyStateAfterLastChiArbitrarily(Trail& trail) const;
private:
    /** This method initializes affinePerInput.
      */
    void initializeAffine();
    /** This method initializes weightPerSlice.
      */
    void initializeWeight();
    /** This method initializes minReverseWeightPerSlice.
      */
    void initializeMinReverseWeight();
    /** This method initializes chiCompatibilityTable.
      */
    void initializeChiCompatibilityTable();
    unsigned int weightOfSlice(SliceValue slice) const;
};

/** This class implements an iterator over the possible state values
  * before χ given a state after χ.
  */
class ReverseStateIterator
{
private:
    vector<ListOfRowPatterns> patterns;
    vector<unsigned int> Ys, Zs;
    vector<unsigned int> indexes;
    unsigned int size;
    unsigned int minWeight, maxWeight;
    vector<SliceValue> current;
    unsigned int currentWeight;
    UINT64 index;
    bool end;
public:
    /** This constructor initializes the iterator based on a state value after χ,
      * the KeccakFPropagation instance, which determines the compatible states.
      * With this constructor, the iterator runs through all the possible states
      * regardless of their weight.
      * @param   stateAfterChi  The state value after χ as a vector of slices.
      * @param   DCorLC         A reference to the KeccakFPropagation instance that
      *                         determines the type of propagation.
      */
    ReverseStateIterator(const vector<SliceValue>& stateAfterChi, const KeccakFPropagation& DCorLC);
    /** This constructor initializes the iterator based on a state value after χ,
      * the KeccakFPropagation instance, which determines the compatible states,
      * and a maximum of the propagation weight.
      * @param   stateAfterChi  The state value after χ as a vector of slices.
      * @param   DCorLC         A reference to the KeccakFPropagation instance that
      *                         determines the type of propagation.
      * @param   aMaxWeight     The iterator will run through the states whose propagation
      *                         weight is not higher than this parameter. 
      */
    ReverseStateIterator(const vector<SliceValue>& stateAfterChi, const KeccakFPropagation& DCorLC, unsigned int aMaxWeight);
    /** This method tells whether the iterator has reached the end of the possible states.
      * @return It returns true iff there are no more states to run through.
      */
    bool isEnd() const;
    /** This method tells wether the set of states to run through is empty.
      * @return It returns true iff there are no states to run through.
      */
    bool isEmpty() const;
    /** This method moves the iterator to the next state. */
    void operator++();
    /** This method returns a constant reference to the current state.
      * @return A constant reference to the current state as a vector of slices.
      */
    const vector<SliceValue>& operator*() const;
    /** This method returns the propagation weight of the current state.
      * @return The weight of the current state.
      */
    unsigned int getCurrentWeight() const;
private:
    void initialize(const vector<SliceValue>& stateAfterChi, const KeccakFPropagation& DCorLC);
    void next();
};

#endif
