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

#ifndef _KECCAKFTRAILCORE3ROUNDS_H_
#define _KECCAKFTRAILCORE3ROUNDS_H_

#include "Keccak-fPositions.h"
#include "Keccak-fPropagation.h"
#include "Keccak-fState.h"
#include "Keccak-fDCLC.h"
#include <stack>
#include <set>


/** Struct that groups information about knots, useful in generating light 3-round trails. 
  */
struct KnotInformation {
   /** Whether the knot is an orbital. 
     */
    bool isOrbital;                      
   /** Number of active rows in the knot. 
     */
    unsigned int nrActiveRows;           
   /** Mininum number of knot points to add to make this a tame knot. 
     */
    unsigned int knotPointDeficit;       
   /** Minimum increase in weight due to the knot points that must be added before this knot is tame.
     */
    unsigned int knotWeightAtBDeficit; 
};

/** Class that implements an iterator over trail seeds for light 3-round trails, given a background. 
  */
class TrailCore3Rounds  : public KeccakFPropagation
{
protected:

    /** Vector with element i containing information about the slice with value i, relevant when it is a knot.
      */
    vector<UINT8> knotInfoLUT;

    /** Maximum lower weight of 3-round trail cores to be generated.
      */
    unsigned int maxWeight;

    /** Keeps track of the knot slices by their z coordinate and relevant information related to them.
      */
    map<unsigned int,KnotInformation> knots;

    /** Contains the z coordinates of the knots that are there due to the background. 
      * They qualify as tame knots if they contain an orbital while knots without background don't.
      */
    set<unsigned int> knotsWithBackground;
    
    /** Vector of chains, where each chain is actually a vector of bit positions at B.
      * All chains are between knots, except the last chain that may be incomplete.
      * If the last chain is incomplete, it is called the working chain.
      * The chains respect an order: 
      * - later chains are never shorter than earlier chains, with the exception of the working chain
      * - when chains have the same length, the starting bit position of the earlier chain is smaller than that of the later.
      */
    vector<vector<BitPosition> > chains;

    /** y-offset accompanying the chains, required in the iteration process.
      * yOffset[i][j] gives information on how to compute chains[i][j+1] from chains[i][j]. 
      * If j is even, this is the y offset at A, otherwise it is the y offset at B.
      * For the last point of working chain, yOffset[i].back() can be used as the current index when extending it. 
      */
    vector<vector<unsigned int> > yOffsets;

    /** This stack keeps track of whether a knot point added a knot (true) or not (false).
      * Adding a chain may add a knot if it starts in a slice that is not a knot yet.
      * Ending a chain in a slice that is not a knot may also add a knot.
      * Chain start and end points are called knot points.
      */
    stack<bool> knotPointAddedKnot;

    /** The state at A.
      * - If DC: stateAtA is λ<sup>-1</sup>(stateAtB).
      * - If LC: stateAtA is λ(stateAtB).
      */
    vector<SliceValue> stateAtA;

    /** The state at B.
      * - If DC: stateAtB is λ(stateAtA).
      * - If LC: stateAtB is λ<sup>-1</sup>(stateAtA).
      */
    vector<SliceValue> stateAtB;

    /** This attribute is a mask indicating bit positions at B where no active bits may exist.
      * This is related to the background at A: an active bit in one of these positions would modify the background.
      */
    vector<SliceValue> tabooAtB;

    /** This attribute keeps track of the minimum length of the last chain.
      * It is redundant and there for optimization.
      */
    unsigned int minimumWorkingChainLength;
    
    /** This attribute keeps track of whether the start point of the working chain is free.
      * It is redundant and there for optimization.
      */
    bool startPointWorkingChainIsFree;

    /** This attribute keeps track of the number of active rows in stateAtA.
      * It is redundant and there for optimization.
      */
    unsigned int nrActiveRowsAtA; 

    /** This attribute keeps track of the Hamming weight of stateAtA.
      * It is redundant and there for optimization.
      */
    unsigned int hammingWeightAtA;

    /** Struct that groups information about a vortex, useful when adding vortices to tame tame states. 
      */
    struct VortexInfo {
        /** State at B. 
          */
        SparseStateAsSlices stateAtB;
        /** Number of active rows at A. 
          */
        unsigned int nrActiveRowsAtA;
        /** Number of active rows at D, assuming C is in the kernel.
          */
        unsigned int nrActiveRowsAtD;
    };

    /** Database of all vortices that lead to 3-round trail cores with weight not above maxWeight.
      * Vortex[u][v] contains a vortex of length 2*u, v is just an index and has no special significance.
      */
    vector<vector<VortexInfo> > vortexBase;

protected:
    /** This method builds the table knotInfoLUT 
      */
    void initializeKnotInfoLUT();

    /** This method codes the knot information in 8 bits for storage.
      * @return A byte coding the given parameters.
      */
    UINT8 packKnotInfo(unsigned int knotPointDeficit, unsigned int knotWeightAtBDeficit, unsigned int nrActiveRows, bool isOrbital) const;

    /** This method returns a Boolean indicating whether the current weight allows adding a point 
      * to the working chain.
      * @return        Whether a point may be added to the working chain.
      */
    bool canAffordExtendingChain() const;

    /** This method returns a Boolean indicating whether the current weight allows adding a chain (true) or not (false). 
      * @return        Whether a chain may be added.
      */
    bool canAffordAddingChain() const;

    /** This method initializes the state with the background in the parameters.
      * After executing this function the states at A and B are consistent through λ.
      * @param   backgroundAtA  The background at A to be put into the state.
      */
    void populateStatesWithBackground(const vector<SliceValue>& backgroundAtA);

    /** This method adds a point at the end of a working chain that wants to become a vortex.
      * @param   pB         Bit position to be added.
      * @param   chainAtB   Working chain with Bit Positions at B.
      * @param   yOffset    Y offsets of the chain.
      * @param   rowsAtA    Keeps track of the active row positions at A.
      * @param   rowsAtD    Keeps track of the active row positions at D.
      */
     void addVortexPoint(const BitPosition& pB,
                        vector<BitPosition>& chainAtB,
                        vector<unsigned int>& yOffset,
                        map<RowPosition,unsigned int>& rowsAtA,
                        map<RowPosition,unsigned int>& rowsAtD,
                        map<unsigned int,unsigned int>& slicesAtB) const;

    /** This method removes a point from the back of a working chain that wants to become a vortex.
      * @param   chainAtB   Working chain with Bit Positions at B .
      * @param   yOffset    Y offsets of the chain.
      * @param   rowsAtA    Keeps track of the active row positions at A.
      * @param   rowsAtD    Keeps track of the active row positions at D.
      */
    void removeVortexPoint(vector<BitPosition >& pointsAtB,
                           vector<unsigned int>& yOffset,
                           map<RowPosition,unsigned int>& rowsAtA,
                           map<RowPosition,unsigned int>& rowsAtD,
                           map<unsigned int,unsigned int>& slicesAtB) const;

    /** This method adds a vortex to vortexBase, provided the corresponding state at B is minimal in z
      * @param   chainAtB  Chain containing the coordinates of the vortex points at B.
      * @param   nrActiveRowsAtA  self-explanatory.
      * @param   nrActiveRowsAtD  self-explanatory.
      */
    void addVortexToBaseIfMinimal(const vector<BitPosition>& pointsAtB,
                                  unsigned int nrActiveRowsAtA,
                                  unsigned int nrActiveRowsAtD);

    /** This method assures the vortexBase contains all vortices up to and including maxWeight.
      */
    void initializeVortexBase();

    /** This virtual method returns whether the state is well formed corresponding to the subclass.
      * @return        Whether the state is well formed.
      */ 
    virtual bool isStateAtBWellFormed() const = 0;

    /** This virtual method indicates whether iteration can be continued given some assumptions that are specified 
      * in the parameters.
      * @param   deltaNrKnotPointsWorkingChain      Number of knot points that will be added by the working chain.
      * @param   minDeltaNrKnotsOrRuns              Minimum number of knots (in kernel at C) or runs (out kernel at C) that will be added.
      * @param   nrOrbitalPointsPerDeltaChain       Minimum number of orbital points that each chain to be added must contain.
      * @param   deltaNrOrbitalPointsWorkingChain   Minimum number of orbital points that must be added to the working chain.
      * @return       Whether the iteration can be continued.
      */
    virtual bool canAffordGeneric(  unsigned int deltaNrKnotPointsWorkingChain,
                                    unsigned int minDeltaNrKnotsOrRuns,
                                    unsigned int nrOrbitalPointsPerDeltaChain,
                                    unsigned int deltaNrOrbitalPointsWorkingChain) const = 0;

    /** This virtual method adds the point with bit position at B given in parameter.
      * It also adapts all the weight stacks. 
      * As the processing is different between knot points and orbital points, this is indicated in parameter.
      * @param   pB             Bit position at B of the point to be added.
      * @param   toKnotSlice    Whether it is a knot point (true) or an orbital point (false).
      * @param   isBackgroundPoint  Whether the point added is a background point (true) or not (false).
      */
    virtual void addPoint(const BitPosition& pB, bool toKnotSlice, bool isBackgroundPoint = false) = 0;

    /** This virtual method removes the last point of the working chain.
      * It also adapts all the weight stacks.
      * As the processing is different between knot points and orbital points, this is indicated in parameter.
      * @param   fromKnotSlice  Whether it is a knot point (true) or an orbital point (false).
      */
    virtual void removePoint(bool fromKnotSlice) = 0;

    /** This virtual method converts the end point (a knot point) of the working chain to an orbital point.
      * It also adapts all the weight stacks.
      */
    virtual void convertKnotPointToOrbitalPoint() = 0;

    /** This virtual method indicates whether the bit position in parameter may be 
      * an end point given the current weights and maxWeight.
      * @return       Whether the point in parameter may be an endpoint
      */
    virtual bool mayBeEndPoint(const BitPosition& pB) = 0;

    /** This virtual method attempts to find the z coordinate larger than that of in parameter (or the first one) 
      * from which the working chain may start.
      * If successful, the z coordinate in parameter will have the correct z value and the method returns true.
      * If unsuccessful, it returns false.
      * @param   z      The z coordinate where the working chain may start.
      *                             If @a zIsInitialized is true, this parameter contains 
      *                             the current bit position from where the search 
      *                             for a new start slice starts.
      * @param  zIsInitialized  Whether z is already initialized. 
      *                         If set to false, it looks for the first z coordinate.
      * @return        Whether a next chain may start.
      */
    virtual bool mayBeStartPointSliceAndGoThere(unsigned int& z, bool zIsInitialized) = 0;

    /** This method attempts to complete the working chain. When calling it, the last point of the working chain is not an end point, so it is not in a knot.
      * if unsuccessful, the working chain has a single point only: its start point. This start point will have the same value as when this method was called.
      * @return        Whether completing the chain was successful.
      */
    bool completeChain();

    /** This method iterates the working chain to the next one compatible with maxWeight.
      * It successful, it returns true. 
      * if unsuccessful, it returns false and the working chain is removed.
      * @return        Whether finding a next chain was successful.
      */
    bool nextChain();

    /** This method updates the minimum required length of the working chain.
      */
    void updateMinimumWorkingChainLength();

    /** This method attempts to find the next 3-round trail core state with at least one knot, respecting the maximum weight.
      * if unsuccessful, the number of chains is zero.
      * @return        Whether finding a next 3-round trail core state was successful.
      */ 
    bool nextWithKnots();

    /** This method populates the knot information instance in parameter based on the slice value in parameter.
      * It takes into account whether it has a background knot point or not.
      * @param   aKnotInformation   Knot information object to be populated.
      * @param   aSliceValue        Slice value of the knot to be reported on in the knot information.
      * @param   knotHasSinglePoint Indicates whether the knot has a single point.
      * @param   hasBackground      Indicates whether the knot has a background knot point or not.
      */ 
    void populateKnotInfo(KnotInformation& aKnotInformation,const SliceValue& aSliceValue, bool knotHasSinglePoint, bool hasBackground) const;

public:

    /** Constructor that initializes the attributes. 
      * An initial call to next() is necessary for having the first valid state.
      * @param   backgroundAtA  The background at A.
      * @param   aTabooAtB      State indicating the bits at B where no active points may be put.
      * @param   aMaxWeight     The maximum propagation weight up to which one must generate all trail cores.
      * @param   aP             Reference to KeccakFTamePropagation instance defining laneSize and LC/DC.
      */
    TrailCore3Rounds(const vector<SliceValue>& backgroundAtA,
                    const vector<SliceValue>& aTabooAtB,
                    unsigned int aMaxWeight,
                    const KeccakFDCLC& aParent,
                    KeccakFPropagation::DCorLC aDCorLC);
};

#endif
