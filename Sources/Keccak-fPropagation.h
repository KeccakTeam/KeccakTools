/*
Tools for the Keccak sponge function family.
Authors: Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche

This code is hereby put in the public domain. It is given as is,
without any guarantee.

For more information, feedback or questions, please refer to our website:
http://keccak.noekeon.org/
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
    /** This method multiplies the vector (dx, dy)<sup>t</sup> by the matrix π or π<sup>-1</sup>
      * to the left:
      * - for DC, π is used; 
      * - for LC, π<sup>-1</sup> is used.
      */
    void directPi(unsigned int& dx, unsigned int& dy) const;
    /** This method multiplies the vector (dx, dy)<sup>t</sup> by the matrix π or π<sup>-1</sup> 
      * to the left:
      * - for DC, π<sup>-1</sup> is used; 
      * - for LC, π is used.
      */
    void reversePi(unsigned int& dx, unsigned int& dy) const;
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
    /** This method returns true iff θ (or θ<sup>t</sup>) is the first step of the linear step between two χ's.
      * @return It returns true iff θ (or θ<sup>t</sup>) is the first step of the linear step between two χ's.
      */
    bool isThetaJustAfterChi() const;
    /** This method applies λ in the "direct" direction:
      * - DC: θ then ρ then π;
      * - LC: π<sup>-1</sup> then ρ<sup>-1</sup> then θ<sup>t</sup>.
      * @param   in     The input state value given as a vector of slices.
      * @param   out    The output state value returned as a vector of slices.
      */
    void directLambda(const vector<SliceValue>& in, vector<SliceValue>& out) const;
    /** This method applies λ in the "reverse" direction:
      * - DC: π<sup>-1</sup> then ρ<sup>-1</sup> then θ<sup>-1</sup>;
      * - LC: θ<sup>-1t</sup> then ρ then π.
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
      * the KeccakFPropagation instance, which determines the compatible states,
      * and an optional maximum of the propagation weight.
      * @param   stateAfterChi  The state value after χ as a vector of slices.
      * @param   DCorLC         A reference to the KeccakFPropagation instance that
      *                         determines the type of propagation.
      * @param   aMaxWeight     The iterator will run through the states whose propagation
      *                         weight is not higher than this parameter. 
      *                         If 0, the iterator runs through all the possible states.
      */
    ReverseStateIterator(const vector<SliceValue>& stateAfterChi, const KeccakFPropagation& DCorLC, unsigned int aMaxWeight = 0);
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
    void next();
};

#endif
