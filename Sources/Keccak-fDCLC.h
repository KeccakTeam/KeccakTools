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

#ifndef _KECCAKFDCLC_H_
#define _KECCAKFDCLC_H_

#include <string>
#include "Keccak-fParts.h"
#include "Keccak-fTrails.h"

class KeccakFPropagation;

/** As a utility class for KeccakFDCLC, this class lists the output row patterns
  * (differences or linear masks) compatible with a given row input pattern
  * through either χ or χ<sup>-1</sup>.
  * The main attribute is @a values, which lists all the corresponding output row patterns.
  * Associated with each output value (same index in the vector),
  * its propagation weight is in the vector @a weights.
  * Note that the values are sorted in increasing weight order,
  * so one can start looking for the outputs with lowest weight first.
  * The maxWeight and minWeight attributes indicate the maximum and minimum
  * propagation weights present, respectively.
  */
class ListOfRowPatterns {
public:
    /** The list of output row patterns (differences or selection vectors)
      * compatible with a given row input pattern.
      */
    vector<RowValue> values;
    /** The list of propagation weights, for each row pattern, i.e.,
      * weights[i] is the weight of values[i].
      */
    vector<int> weights;
    /** True iff maxWeight and minWeight are meaningful. */
    bool minMaxInitialized;
    /** The maximum propagation weight in weights.*/
    int maxWeight;
    /** The minimum propagation weight in weights.*/
    int minWeight;
public:
    /** This constructor initializes the output patterns to the empty set.
      */
    ListOfRowPatterns() : minMaxInitialized(false) {}
    /** This function adds a possible output pattern,
      * along with an associated propagation weight,
      * while inserting it at the right
      * place to ensure that they are listed with increasing propagation weights.
      * @param  aValue  An output pattern.
      * @param  aWeight Its associated propagation weight.
      */
    void add(RowValue aValue, int aWeight);
    /** This function displays the possible patterns and their weights.
      * @param  fout    The stream to display to.
      */
    void display(ostream& fout) const;
};

/** This class is an extension of KeccakF with additional functionality
  * aimed at differential and linear cryptanalysis.
  */
class KeccakFDCLC : public KeccakF {
public:
    /** This attribute contains the output difference patterns for the non-linear function χ.
      * The index of the vector corresponds to the input difference pattern.
      */
    vector<ListOfRowPatterns> diffChi;
    /** This attribute contains the output difference patterns for the non-linear function χ<sup>-1</sup>.
      * The index of the vector corresponds to the input difference pattern.
      */
    vector<ListOfRowPatterns> diffInvChi;
    /** This attribute contains the output linear masks for the non-linear function χ.
      * The index of the vector corresponds to the input mask.
      */
    vector<ListOfRowPatterns> corrChi;
    /** This attribute contains the output linear masks for the non-linear function χ<sup>-1</sup>.
      * The index of the vector corresponds to the input mask.
      */
    vector<ListOfRowPatterns> corrInvChi;
    /** This attribute indicates, for each λ mode (see LambdaMode), if θ is just after χ.
      * When true, this means that lambdaBeforeTheta() is the identity.
      */
    vector<bool> thetaJustAfterChi;
    /** This attribute indicates, for each λ mode (see LambdaMode), if θ is just before χ.
      * When true, this means that lambdaAfterTheta() is the identity.
      */
    vector<bool> thetaJustBeforeChi;
private:
    /** For λ mode (see LambdaMode), this attribute contains
      * the contribution of an input row to an output slice via the linear
      * function λ.
      *
      * The indexes are as follows, from outer to inner:
      * - the λ mode (λ, its inverse, its transpose, the transpose of its inverse)
      * - the output slice index (oz)
      * - the input slice index (iz)
      * - the input row index (iy)
      * - the input row value input[iy,iz]
      *
      * So, given a row value B located at row iy in slice 0&lt;=iz&lt;laneSize,
      * its linear contribution to the output of the λ function in Straight mode is given
      * by lambdaRowToSlice[(int)Straight][oz][iz][iy][B], where oz is
      * varied from 0 to laneSize-1 to see all the output slices.
      */
    vector<vector<vector<vector<vector<SliceValue> > > > > lambdaRowToSlice;
    /** Same as lambdaRowToSlice, but only for the linear part before θ.
      */
    vector<vector<vector<vector<vector<SliceValue> > > > > lambdaBeforeThetaRowToSlice;
    /** Same as lambdaRowToSlice, but only for the linear part after θ.
      */
    vector<vector<vector<vector<vector<SliceValue> > > > > lambdaAfterThetaRowToSlice;
public:
    /** In this context, λ represents the linear operations in Keccak-<i>f</i>
      * between two applications of χ.
      * The λ mode indicates whether to perform these operations, their transpose,
      * their inverse or the transpose of their inverse.
      * - Straight: λ(a) means π(ρ(θ(a)));
      * - Inverse: λ(a) means θ<sup>-1</sup>(ρ<sup>-1</sup>(π<sup>-1</sup>(a)));
      * - Transpose: λ(a) means θ<sup>t</sup>(ρ<sup>t</sup>(π<sup>t</sup>(a)))=θ<sup>t</sup>(ρ<sup>-1</sup>(π<sup>-1</sup>(a)));
      * - Dual: λ(a) means π(ρ(θ<sup>-1t</sup>(a))).
      */
    enum LambdaMode {
        Straight = 0,
        Inverse,
        Transpose,
        Dual,
        EndOfLambdaModes
    };
    /** This constructor has the same parameters as KeccakF::KeccakF.
      */
    KeccakFDCLC(unsigned int aWidth);
    /** This method retuns a string describing the instance. */
    string getDescription() const;
    /** This method displays all the possible output patterns
      * (differences and selection vectors) for each input pattern.
      * If the paramters DC and LC point to actual KeccakFPropagation
      * instances, the affine descriptions are also given.
      * @param  fout    The stream to display to.
      * @param  DC      A pointer to the KeccakFPropagation instance specialized in differential cryptanalysis (can be null).
      * @param  LC      A pointer to the KeccakFPropagation instance specialized in linear cryptanalysis (can be null).
      */
    void displayAll(ostream& fout, KeccakFPropagation *DC=0, KeccakFPropagation *LC=0) const;
    /**
      * Template method that applies θ<sup>t</sup> (θ transposed).
      *
      * @param  A       The given state is organized as a vector of lanes.
      */
    template<class Lane> void thetaTransposed(vector<Lane>& A) const;
    /** This method inverts the order of bits in lanes and of lanes in planes.
      * @param  state   The state to process as a vector of lanes.
      */
    void thetaTransEnvelope(vector<LaneValue>& state) const;
    /** This method applies the linear transformation (see LambdaMode) between two χ's.
      * @param   state  The state to process as a vector of lanes.
      * @param   mode   The λ mode.
      */
    template<class Lane> void lambda(vector<Lane>& state, LambdaMode mode) const;
    /** This method applies the λ function (see LambdaMode)
      * using the lambdaRowToSlice table.
      * The input is a vector of laneSize slice values, and so is the output.
      * The mode argument gives the λ mode (i.e., inverse/transpose) to use.
      * @param   in     The input state as a vector of slices.
      * @param   out    The output state as a vector of slices.
      * @param   mode   The λ mode.
      * @pre This assumes that @a in has size equal to @a laneSize.
      */
    void lambda(const vector<SliceValue>& in, vector<SliceValue>& out, LambdaMode mode) const;
    /** Among the linear transformation steps (see LambdaMode)
      * between two χ's, this method applies the linear steps that come before θ.
      * @param   state  The state to process as a vector of lanes.
      * @param   mode   The λ mode.
      */
    template<class Lane> void lambdaBeforeTheta(vector<Lane>& state, LambdaMode mode) const;
    /** This method is the same as lambdaBeforeTheta() but works on states represented as slices
      * and internally uses the lambdaBeforeThetaRowToSlice table.
      * @param   in     The input state as a vector of slices.
      * @param   out    The output state as a vector of slices.
      * @param   mode   The λ mode.
      * @pre This assumes that @a in has size equal to @a laneSize.
      */
    void lambdaBeforeTheta(const vector<SliceValue>& in, vector<SliceValue>& out, LambdaMode mode) const;
    /** Among the linear transformation steps (or its inverse and/or transpose)
      * between two χ's, this method applies the linear steps that come after θ.
      * @param   state  The state to process as a vector of lanes.
      * @param   mode   The λ mode.
      */
    template<class Lane> void lambdaAfterTheta(vector<Lane>& state, LambdaMode mode) const;
    /** This method is the same as lambdaAfterTheta() but works on states represented as slices
      * and internally uses the lambdaAfterThetaRowToSlice table.
      * @param   in     The input state as a vector of slices.
      * @param   out    The output state as a vector of slices.
      * @param   mode   The λ mode.
      * @pre This assumes that @a in has size equal to @a laneSize.
      */
    void lambdaAfterTheta(const vector<SliceValue>& in, vector<SliceValue>& out, LambdaMode mode) const;
    /** Apply the χ function on a single row value.
      * @param   a      The input row value.
      * @return The output row value.
      */
    RowValue chiOnRow(RowValue a) const;
    /** Apply the χ<sup>-1</sup> function on a single row value.
      * @param   a      The input row value.
      * @return The output row value.
      */
    RowValue inverseChiOnRow(RowValue a) const;
    /** This method creates the value of a state represented as a vector of slices
      * from a state represented as a vector of lanes.
      * @param  lanes   The state as a vector of lanes.
      * @param  slices  The output state as a vector of slices.
      */
    inline void fromLanesToSlices(const vector<LaneValue>& lanes, vector<SliceValue>& slices) const
    {
        ::fromLanesToSlices(lanes, slices, laneSize);
    }
    /** This method checks the consistency of a DC trail given as parameter.
      * If an inconsistency is found, a KeccakException is thrown and details
      * about the inconsistency are displayed in the error console (cerr).
      * The aspects tested are:
      * - the propagation weights declared in the trail match the propagation weights of the specified differences;
      * - between two rounds, the specified differences are compatible.
      * .
      * @param  trail   The trail to test the consistence of.
      * @param  DC      A pointer to the KeccakFPropagation instance specialized in differential cryptanalysis.
      *                 This pointer can be zero. It is only used for display purposes, in case
      *                 an inconsistency is detected.
      */
    void checkDCTrail(const Trail& trail, KeccakFPropagation *DC=0) const;
    /** This method checks the consistency of a LC trail given as parameter.
      * If an inconsistency is found, a KeccakException is thrown and details
      * about the inconsistency are displayed in the error console (cerr).
      * The aspects tested are:
      * - the propagation weights declared in the trail match the propagation weights of the specified masks;
      * - between two rounds, the specified masks are compatible.
      * .
      * @param  trail   The trail to test the consistence of.
      * @param  LC      A pointer to the KeccakFPropagation instance specialized in linear cryptanalysis.
      *                 This pointer can be zero. It is only used for display purposes, in case
      *                 an inconsistency is detected.
      */
    void checkLCTrail(const Trail& trail, KeccakFPropagation *LC=0) const;
    /** This method computes the θ-gap of the state given as input.
      * (See the Keccak main document for a definition of the θ-gap.)
      * @param   state  The state of which to compute the θ-gap as a vector of lanes.
      * @return The θ-gap value.
      */
    unsigned int getThetaGap(const vector<LaneValue>& state) const;
    /** This method computes the θ-gap from the parity of a state given as input.
      * (See the Keccak main document for a definition of the θ-gap.)
      * @param   parity The parity of which to compute the θ-gap as a vector of 5 lanes.
      * @return The θ-gap value.
      */
    unsigned int getThetaGapFromParity(const vector<LaneValue>& parity) const;
    /** This method computes the θ-effect of a state given as input.
      * (See the Keccak main document for a definition of the θ-effect.)
      * @param   C  The parity as a vector of 5 lanes.
      * @param   D  The θ-effect as a vector of 5 lanes.
      */
    template<class Lane> void getThetaEffectFromParity(const vector<Lane>& C, vector<Lane>& D) const;
    /** This method computes the θ<sup>t</sup>-effect of a state given as input.
      * (The θ<sup>t</sup>-effect is defined like the θ-effect but on the transposed θ.)
      * @param   C  The parity as a vector of 5 lanes.
      * @param   D  The θ<sup>t</sup>-effect as a vector of 5 lanes.
      */
    template<class Lane> void getThetaTransposedEffectFromParity(const vector<Lane>& C, vector<Lane>& D) const;
    /**
      * Method that returns a short string that uniquely identifies the
      * Keccak-<i>f</i> instance.
      * Unlike the function in KeccakF, this function does not take the number
      * of rounds into account in the description.
      * @return The name of the instance.
      */
    virtual string getName() const;
private:
    void initializeAll();
    void initializeLambdaLookupTables();
};

template<class Lane>
void KeccakFDCLC::thetaTransposed(vector<Lane>& A) const
{
    vector<Lane> C(5);
    for(unsigned int x=0; x<5; x++) {
        C[x] = A[index(x,0)];
        for(unsigned int y=1; y<5; y++)
            C[x] ^= A[index(x,y)];
    }
    vector<Lane> D(5);
    for(unsigned int x=0; x<5; x++) {
        Lane temp = C[index(x-1)];
        ROL(temp, -1);
        D[x] = temp ^ C[index(x+1)];
    }
    for(int x=0; x<5; x++)
        for(unsigned int y=0; y<5; y++)
            A[index(x,y)] ^= D[x];
}

template<class Lane>
void KeccakFDCLC::lambda(vector<Lane>& state, LambdaMode mode) const
{
    if (mode == Straight) {
        theta(state);
        rho(state);
        pi(state);
    }
    else if (mode == Inverse) {
        inversePi(state);
        inverseRho(state);
        inverseTheta(state);
    }
    else if (mode == Transpose) {
        inversePi(state);
        inverseRho(state);
        thetaTransposed(state);
    }
    else if (mode == Dual) {
        thetaTransEnvelope(state);
        inverseTheta(state);
        thetaTransEnvelope(state);
        rho(state);
        pi(state);
    }
}

template<class Lane>
void KeccakFDCLC::lambdaBeforeTheta(vector<Lane>& state, LambdaMode mode) const
{
    if ((mode == Transpose) || (mode == Inverse)) {
        inversePi(state);
        inverseRho(state);
    }
}

template<class Lane>
void KeccakFDCLC::lambdaAfterTheta(vector<Lane>& state, LambdaMode mode) const
{
    if ((mode == Straight) || (mode == Dual)) {
        rho(state);
        pi(state);
    }
}

template<class Lane>
void KeccakFDCLC::getThetaEffectFromParity(const vector<Lane>& C, vector<Lane>& D) const
{
    if (D.size() != 5)
        D.resize(5);
    for(unsigned int x=0; x<5; x++) {
        Lane temp = C[index(x+1)];
        ROL(temp, 1);
        D[x] = temp ^ C[index(x-1)];
    }
}

template<class Lane>
void KeccakFDCLC::getThetaTransposedEffectFromParity(const vector<Lane>& C, vector<Lane>& D) const
{
    if (D.size() != 5)
        D.resize(5);
    for(unsigned int x=0; x<5; x++) {
        Lane temp = C[index(x-1)];
        ROL(temp, -1);
        D[x] = temp ^ C[index(x+1)];
    }
}

#endif
