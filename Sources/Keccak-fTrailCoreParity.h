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

#ifndef _KECCAKFTRAILCOREPARITY_H_
#define _KECCAKFTRAILCOREPARITY_H_

#include <stack>
#include <vector>
#include "Keccak-fPositions.h"
#include "Keccak-fPropagation.h"
#include "Keccak-fTrails.h"
#include "Keccak-fState.h"
#include "progress.h"

typedef vector<SliceValue> StateAsVectorOfSlices;

/** Class containing the column position and y-coordinates of the two bits in an orbital. */
class OrbitalPosition : public ColumnPosition
{
public:
    /** The y-coordinate of the first bit in the orbital. */
    unsigned int y0;
    /** The y-coordinate of the second bit in the orbital. */
    unsigned int y1;
public:
    /** The default constructor. */
    OrbitalPosition() : ColumnPosition(), y0(0), y1(1) {}
    /** This method sets the position to the first available one with y-coordinates
      * at least as specified by the parameter yMin.
      * @param  yMin    The value in yMin[x+5*z] indicates the minimum y-coordinate of the bits in the orbital.
      * @param  laneSize    The lane size.
      */
    bool first(const vector<unsigned int>& yMin, unsigned int laneSize);
    /** This method sets the position to the next available one with y-coordinates
      * at least as specified by the parameter yMin.
      * @param  yMin    The value in yMin[x+5*z] indicates the minimum y-coordinate of the bits in the orbital.
      * @param  laneSize    The lane size.
      */
    bool next(const vector<unsigned int>& yMin, unsigned int laneSize);
    /** This method sets the position to the next available orbital after @a other with y-coordinates
      * at least as specified by the parameter yMin and with y-coordinates higher than those
      * of the specified orbital @a other if in the same column.
      * @param  other   The other orbital position.
      * @param  yMin    The value in yMin[x+5*z] indicates the minimum y-coordinate of the bits in the orbital.
      * @param  laneSize    The lane size.
      */
    bool successorOf(const OrbitalPosition& other, const vector<unsigned int>& yMin, unsigned int laneSize);
};

/** This abstract class iterates on all 2-round trail cores with a given parity, 
  * limited by some budget defined in a derived class.
  */
class KeccakFTrailWithGivenParityIterator : public TrailIterator
{
protected:
    unsigned int laneSize;
    bool orbitals;
    vector<RowValue> C, D;
    vector<ColumnPosition> UOcolumns, Acolumns;
    bool initialized;
    bool end, empty;
    UINT64 index;
    Trail trail;

    void initialize();

    // Common part
    /** This abstract method is called when the iteration adds active bits
      * in an affected column.
      * It returns true iff the budget allows this addition.
      * If this method returns true, the pop() method will later be called to remove these active bits.
      * If this method returns false, the caller understands that nothing changed.
      * @param  columnBeforeTheta   The column position before θ.
      * @param  value   The value of the column before θ.
      * @return True iff the active bits could and were added.
      */
    virtual bool pushValueInAffectedColumn(const ColumnPosition& columnBeforeTheta, ColumnValue value) = 0;
    /** This abstract method is called when the iteration adds an active bit
      * in an unaffected odd column.
      * It returns true iff the budget allows this addition.
      * If this method returns true, the pop() method will later be called to remove the active bit.
      * If this method returns false, the caller understands that nothing changed.
      * @param  columnBeforeTheta   The column position before θ.
      * @param  y   The y coordinate of the bit before θ.
      * @return True iff the active bit could and was added.
      */
    virtual bool pushBitInUnaffectedOddColumn(const ColumnPosition& columnBeforeTheta, unsigned int y) = 0;
    /** This abstract method is called when the iteration adds 
      * an active pair of bits, called an orbital.
      * It returns true iff the budget allows this addition.
      * If this method returns true, the pop() method will later be called to remove the active bits.
      * If this method returns false, the caller understands that nothing changed.
      * @param  orbital The orbital position before θ.
      * @return True iff the active bits could and were added.
      */
    virtual bool pushOrbitalInUnaffectedColumn(const OrbitalPosition& orbital) = 0;
    /** This abstract method is called when the iteration removes active bits
      * added, corresponding to a call either 
      * to pushValueInAffectedColumn(), pushBitInUnaffectedOddColumn() or 
      * pushOrbitalInUnaffectedColumn(), as on a stack.
      */
    virtual void pop() = 0;
    /** This abstract method is called to populate the @a trail attribute.
     */
    virtual void getTrail() = 0;
    bool first();
    bool next();
    
    // Stack 1 for affected columns
    unsigned int S1_height;
    stack<unsigned int> S1_valueIndex;
    static const ColumnValue evenValues[16], oddValues[16];
    bool S1_push(unsigned int valueIndex);
    unsigned int S1_pop();
    bool S1_firstTop();
    bool S1_nextTop();
    bool S1_first();
    bool S1_next();
    
    // Stack 2 for unaffected odd columns
    unsigned int S2_height;
    stack<unsigned int> S2_y;
    bool S2_push(unsigned int y);
    unsigned int S2_pop();
    bool S2_firstTop();
    bool S2_nextTop();
    bool S2_first();
    bool S2_next();
    
    // Stack 3 for orbitals
    unsigned int S3_height;
    vector<unsigned int> S3_yMin;
    stack<OrbitalPosition> S3_position;
    bool S3_push(const OrbitalPosition& orbital);
    bool S3_addNewOrbital();
    bool S3_nextTop();
    bool S3_next();

public:
    /** The constructor.
      * @param   DCorLC The propagation context of the trail, 
      *                 as a reference to a KeccakFPropagation object.
      * @param  aParity The parity as a vector of row values.
      * @param  aOrbitals   If false, only the states with minimum total 
      *     Hamming weight (before and after θ) are generated, i.e.,
      *     no orbitals are added.
      */
    KeccakFTrailWithGivenParityIterator(const KeccakFPropagation& aDCorLC, 
        const vector<RowValue>& aParity, bool aOrbitals = true);

    /** See TrailIterator::isEnd(). */
    bool isEnd();
    /** See TrailIterator::isEmpty(). */
    bool isEmpty();
    /** See TrailIterator::isBounded(). */
    bool isBounded();
    /** See TrailIterator::getIndex(). */
    UINT64 getIndex();
    /** See TrailIterator::getCount(). */
    UINT64 getCount();
    /** See TrailIterator::operator++(). */
    void operator++();
    /** See TrailIterator::operator*(). */
    const Trail& operator*();
};

/** This class iterates on all 2-round trail cores with a given parity, 
 * such that the weight of the core is not greater than a given value.
 */
class KeccakFTwoRoundTrailCoreWithGivenParityIterator : public KeccakFTrailWithGivenParityIterator 
{
protected:
    int maxWeight;
    stack<StateAsVectorOfSlices> stack_stateAtA, stack_stateAtB;
    stack<unsigned int> stack_weight;
public: 
    /** The constructor.
      * @param   DCorLC The propagation context of the trail, 
      *                 as a reference to a KeccakFPropagation object.
      * @param  aParity The parity as a vector of row values.
      * @param  aOrbitals   See KeccakFTrailWithGivenParityIterator::KeccakFTrailWithGivenParityIterator().
      * @param  aMaxWeight The maximum trail core weight w(a)+w<sup>rev</sup>(λ<sup>-1</sup>(a)).
      */
    KeccakFTwoRoundTrailCoreWithGivenParityIterator(const KeccakFPropagation& aDCorLC,
        const vector<RowValue>& aParity, int aMaxWeight, bool aOrbitals = true);
private:
    int setBitToOneAndGetDeltaWeight(vector<SliceValue>& state, const BitPosition& p) const;
    int setBitToOneAndGetDeltaMinReverseWeight(vector<SliceValue>& state, const BitPosition& p) const;
    int setValueInAffectedColumnAndGetDeltaTotalWeight(vector<SliceValue>& stateAtA, vector<SliceValue>& stateAtB, const ColumnPosition& columnBeforeTheta, ColumnValue valueBeforeTheta) const;
    int setBitInUnaffectedColumnAndGetDeltaTotalWeight(vector<SliceValue>& stateAtA, vector<SliceValue>& stateAtB, unsigned int x, unsigned int y, unsigned int z) const;
    int setBitInUnaffectedColumnAndGetDeltaTotalWeight(vector<SliceValue>& stateAtA, vector<SliceValue>& stateAtB, const ColumnPosition& columnBeforeTheta, unsigned int y) const;
    int setOrbitalInUnaffectedColumnAndGetDeltaTotalWeight(vector<SliceValue>& stateAtA, vector<SliceValue>& stateAtB, const OrbitalPosition& orbital) const;
protected:
    bool pushValueInAffectedColumn(const ColumnPosition& columnBeforeTheta, ColumnValue valueBeforeTheta);
    bool pushBitInUnaffectedOddColumn(const ColumnPosition& columnBeforeTheta, unsigned int y);
    bool pushOrbitalInUnaffectedColumn(const OrbitalPosition& orbital);
    void pop();
    void getTrail();
};

#endif
