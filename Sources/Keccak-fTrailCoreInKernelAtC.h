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

#ifndef _KECCAKFTRAILCOREINKERNELATC_H_
#define _KECCAKFTRAILCOREINKERNELATC_H_

#include "Keccak-fPropagation.h"
#include "Keccak-fTrailCore3Rounds.h"
#include "Keccak-fState.h"
#include "Keccak-fDCLC.h"
#include <stack>
#include <set>


/** Class that implements an iterator over states leading to low weight three round trail cores with C in the kernel
  */
class TrailCoreInKernelAtC : public TrailCore3Rounds
{
protected:

    /** State at D that has as active bits the orbital points at B only.
      * - If DC: partialStateAtD is λ(stateAtB, restriced to orbital points).
      * - If LC: partialStateAtD is λ<sup>-1</sup>(stateAtB, restriced to orbital points).
      */
    vector<SliceValue> partialStateAtD;

    /** This attribute keeps track of the propagation weight of stateAtB.
      * It is redundant and there for optimization.
      */
    unsigned int weightAtB;

    /** This attribute keeps track of the minimum number of knot points
      * that must be added before the state at B is tame.
      */
    unsigned int knotPointDeficit;

    /** This attribute keeps track of the increase in weight due to the knot points that
      * must be added before the state at B is tame.
      */
    unsigned int knotWeightAtBDeficit;

    /** This attribute keeps track of the Hamming weight of partialStateAtD.
      * It is redundant and there for optimization.
      */
    unsigned int partialHammingWeightAtD;

    /** This attribute keeps track of the number of active rows in partialStateAtD.
      * It is redundant and there for optimization.
      */
    unsigned int partialNrActiveRowsAtD;

    /** This struct groups data elements useful when adding vortices.
      */
    struct CoreInfo {
       /** Hamming weight of the state at A.
         */
        unsigned int hammingWeightAtA;
       /** Number of active rows of the state at A.
         */
        unsigned int nrActiveRowsAtA;
       /** The state at B.
         */
        vector<SliceValue> stateAtB;
       /** The propagation weight of the state at B.
         */
        unsigned int weightAtB;
       /** A state containing all points at C that will certainly be present in a state at C in the kernel and compatible with the state at B.
         */
        vector<SliceValue> partialStateAtC;
       /** Hamming weight at D based on the partialStateAtC.
         */
        unsigned int hammingWeightAtD;
       /** Number of active rows at D based on the partialStateAtC.
         */
        unsigned int nrActiveRowsAtD;
       /** Weight of the trail core based on the state at B and partialStateAtC.
         */
        unsigned int partialWeight;
       /** Length of the vortex to be added.
         */
        unsigned int vortexLength;
       /** Index of the vortex to be added in vortexBase[vortexLength/2].
         */
        unsigned int vortexIndex;
       /** Z offset by which the vortex to be added must be translated.
         */
        unsigned int vortexZOffset;
    };

    /** This attribute keeps track of the iteration, in the light of adding vortices.
      * When a vortex is added, a new entry is pushed back.
      * The first element contains a state with knots or an empty state.
      */
    vector<CoreInfo> outCore;


protected:

    /** See TrailCore3Rounds::isStateAtBWellFormed().
      */
    bool isStateAtBWellFormed() const;

    /** See TrailCore3Rounds::canAffordGeneric().
      * Here, @a minDeltaNrKnotsOrRuns indicates the number of knots that will be added.
      */
    bool canAffordGeneric(  unsigned int deltaNrKnotPointsWorkingChain,
                            unsigned int minDeltaNrKnotsOrRuns,
                            unsigned int nrOrbitalPointsPerDeltaChain,
                            unsigned int deltaNrOrbitalPointsWorkingChain) const;

    /** See TrailCore3Rounds::addPoint().
      */
    void addPoint(const BitPosition& pB, bool toKnotSlice, bool isBackgroundPoint = false);

    /** See TrailCore3Rounds::removePoint().
      */
    void removePoint(bool fromKnotSlice);

    /** See TrailCore3Rounds::convertKnotPointToOrbitalPoint().
      */
    void convertKnotPointToOrbitalPoint();

    /** See TrailCore3Rounds::mayBeEndPoint().
      */
    bool mayBeEndPoint(const BitPosition& pB);

    /** See TrailCore3Rounds::maybeStartPointSliceAndGoThere().
      */
    bool mayBeStartPointSliceAndGoThere(unsigned int& z, bool zIsInitialized);

    /** This method returns a lower bound for the weight of the state resulting from adding to the state at the back of outCore the
      * vortex specified at the back of outCore. The lower bound is based on the propagation weight at B and the Hamming weights and
      * number of active rows at A and D.
      */
    unsigned int computeLowerWeightAssumingVortexIsAdded();

public:

    /** Constructor that initializes the attributes.
      * An initial call to next() is necessary for having the first valid state.
      * @param   backgroundAtA  The background at A.
      * @param   aTabooAtB      State indicating the bits at B where no active points may be put.
      * @param   aMaxWeight     The maximum propagation weight up to which one must generate all states.
      * @param   aP             Reference to KeccakFTamePropagation instance defining laneSize and LC/DC.
      */
    TrailCoreInKernelAtC(const vector<SliceValue>& backgroundAtA,
                    const vector<SliceValue>& aForbiddenAreaAtB,
                    unsigned int aMaxWeight,
                    const KeccakFDCLC& aParent,
                    KeccakFPropagation::DCorLC aDCorLC);

    /** This method generates the next three-round trail core inside the kernel at C with weight not above maxWeight.
      * @return  Whether a trail core was found.
      */
    bool next();

    /** This method returns a constant reference to the current vortex being iterated.
     */
    const CoreInfo& getTopCoreInfo() const;

    /** This function displays the attributes of the TrailCoreInKernelAtC object.
      */
    friend ostream& operator<<(ostream& fout, const TrailCoreInKernelAtC& aL);
};

#endif
