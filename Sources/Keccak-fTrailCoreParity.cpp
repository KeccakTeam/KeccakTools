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

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stack>
#include <time.h>
#include "Keccak-fTrailCoreParity.h"

bool OrbitalPosition::first(const vector<unsigned int>& yMin, unsigned int laneSize)
{
    x = 0;
    z = 0;
    y0 = yMin[getXplus5Z()];
    while(y0 >= 4) {
        if (x < 4)
            x++;
        else {
            if (z < laneSize-1) {
                z++;
                x = 0;
            }
            else
                return false;
        }
        y0 = yMin[getXplus5Z()];
    }
    y1 = y0+1;
    return true;
}

bool OrbitalPosition::next(const vector<unsigned int>& yMin, unsigned int laneSize)
{
    if (y1 < 4)
        y1++;
    else {
        if (y0 < 3) {
            y0++;
            y1 = y0+1;
        }
        else {
            do {
                if (x < 4)
                    x++;
                else {
                    if (z < laneSize-1) {
                        z++;
                        x = 0;
                    }
                    else
                        return false;
                }
                y0 = yMin[getXplus5Z()];
            } while(y0 >= 4);
            y1 = y0+1;
        }
    }
    return true;
}

bool OrbitalPosition::successorOf(const OrbitalPosition& other, const vector<unsigned int>& yMin, unsigned int laneSize)
{
    x = other.x;
    z = other.z;
    y0 = other.y1+1;
    while(y0 >= 4) {
        if (x < 4)
            x++;
        else {
            if (z < laneSize-1) {
                z++;
                x = 0;
            }
            else
                return false;
        }
        y0 = yMin[getXplus5Z()];
    }
    y1 = y0+1;
    return true;
}

int KeccakFTwoRoundTrailCoreWithGivenParityIterator::setValueInAffectedColumnAndGetDeltaTotalWeight(vector<SliceValue>& stateAtA, vector<SliceValue>& stateAtB, const ColumnPosition& columnBeforeTheta, ColumnValue valueBeforeTheta) const
{
    int delta = 0;
    for(unsigned int y=0; y<5; y++) {
        bool bitBeforeTheta = ((valueBeforeTheta >> y) & 1) != 0;
        BitPosition p(columnBeforeTheta.x, y, columnBeforeTheta.z);
        if (bitBeforeTheta) {
            DCorLC.reverseRhoPiBeforeTheta(p);
            delta += setBitToOneAndGetDeltaMinReverseWeight(stateAtA, p);
        }
        else {
            DCorLC.directRhoPiAfterTheta(p);
            delta += setBitToOneAndGetDeltaWeight(stateAtB, p);
        }
    }
    return delta;
}

int KeccakFTwoRoundTrailCoreWithGivenParityIterator::setBitInUnaffectedColumnAndGetDeltaTotalWeight(vector<SliceValue>& stateAtA, vector<SliceValue>& stateAtB, unsigned int x, unsigned int y, unsigned int z) const
{
    int delta = 0;
    {
        BitPosition p(x, y, z);
        DCorLC.reverseRhoPiBeforeTheta(p);
        delta += setBitToOneAndGetDeltaMinReverseWeight(stateAtA, p);
    }
    {
        BitPosition p(x, y, z);
        DCorLC.directRhoPiAfterTheta(p);
        delta += setBitToOneAndGetDeltaWeight(stateAtB, p);
    }
    return delta;
}

int KeccakFTwoRoundTrailCoreWithGivenParityIterator::setBitInUnaffectedColumnAndGetDeltaTotalWeight(vector<SliceValue>& stateAtA, vector<SliceValue>& stateAtB, const ColumnPosition& columnBeforeTheta, unsigned int y) const
{
    return setBitInUnaffectedColumnAndGetDeltaTotalWeight(stateAtA, stateAtB, columnBeforeTheta.x, y, columnBeforeTheta.z);
}

int KeccakFTwoRoundTrailCoreWithGivenParityIterator::setOrbitalInUnaffectedColumnAndGetDeltaTotalWeight(vector<SliceValue>& stateAtA, vector<SliceValue>& stateAtB, const OrbitalPosition& orbital) const
{
    return 
        setBitInUnaffectedColumnAndGetDeltaTotalWeight(stateAtA, stateAtB, 
            orbital.x, orbital.y0, orbital.z)
    +   setBitInUnaffectedColumnAndGetDeltaTotalWeight(stateAtA, stateAtB, 
            orbital.x, orbital.y1, orbital.z);
}

const ColumnValue KeccakFTrailWithGivenParityIterator::evenValues[16] = {
    0x00, 0x03, 0x05, 0x06, 0x09, 0x0A, 0x0C, 0x0F,
    0x11, 0x12, 0x14, 0x17, 0x18, 0x1B, 0x1D, 0x1E };
const ColumnValue KeccakFTrailWithGivenParityIterator::oddValues[16] = {
    0x01, 0x02, 0x04, 0x07, 0x08, 0x0B, 0x0D, 0x0E,
    0x10, 0x13, 0x15, 0x16, 0x19, 0x1A, 0x1C, 0x1F };

KeccakFTrailWithGivenParityIterator::KeccakFTrailWithGivenParityIterator(const KeccakFPropagation& aDCorLC,
        const vector<RowValue>& aParity, bool aOrbitals)
    : TrailIterator(aDCorLC), laneSize(DCorLC.laneSize), C(aParity), initialized(false), orbitals(aOrbitals)
{
    DCorLC.directThetaEffectFromParities(C, D);
    S3_yMin.assign(5*laneSize, 0);
    for(unsigned int x=0; x<5; x++)
    for(unsigned int z=0; z<laneSize; z++) {
        bool odd = (getBit(C, x, z) != 0);
        bool affected = (getBit(D, x, z) != 0);
        if (affected) {
            ColumnPosition column(x, z);
            Acolumns.push_back(column);
            S3_yMin[column.getXplus5Z()] = 5; // no orbitals here
        }
        else
            if (odd)
                UOcolumns.push_back(ColumnPosition(x, z));
    }
    S1_height = 0;
    S2_height = 0;
    S3_height = 0;
}
   
void KeccakFTrailWithGivenParityIterator::initialize()
{
    index = 0;
    if (first()) {
        getTrail();
        end = false;
        empty = false;
    }
    else {
        end = true;
        empty = true;
    }
    initialized = true;
}

bool KeccakFTrailWithGivenParityIterator::S1_push(unsigned int valueIndex)
{
    bool odd = (getBit(C, Acolumns[S1_height].x, Acolumns[S1_height].z) != 0);
    bool canAfford;
    if (odd)
        canAfford = pushValueInAffectedColumn(Acolumns[S1_height], oddValues[valueIndex]);
    else
        canAfford = pushValueInAffectedColumn(Acolumns[S1_height], evenValues[valueIndex]);
    if (canAfford) {
        S1_valueIndex.push(valueIndex);
        S1_height++;
        return true;
    }
    else
        return false;
}

unsigned int KeccakFTrailWithGivenParityIterator::S1_pop()
{
    pop();
    S1_height--;
    unsigned int valueIndex = S1_valueIndex.top();
    S1_valueIndex.pop();
    return valueIndex;
}

bool KeccakFTrailWithGivenParityIterator::S1_firstTop()
{
    unsigned int valueIndex = 0;
    bool success;
    do {
        success = S1_push(valueIndex);
        valueIndex++;
    } while((!success) && (valueIndex < 16));
    return success;
}

bool KeccakFTrailWithGivenParityIterator::S1_nextTop()
{
    unsigned int valueIndex = S1_pop();
    valueIndex++;
    while(valueIndex < 16) {
        if (S1_push(valueIndex))
            return true;
        valueIndex++;
    }
    return false;
}

bool KeccakFTrailWithGivenParityIterator::S1_first()
{
    while(S1_height < Acolumns.size()) {
        bool success = true;
        while((S1_height < Acolumns.size()) && success) {
            success = S1_firstTop();
        }
        while((S1_height > 0) && (!success)) {
            success = S1_nextTop();
        }
        if (S1_height == 0)
            return false;
    }
    return true;
}

bool KeccakFTrailWithGivenParityIterator::S1_next()
{
    bool success = true;
    do {
        success = S1_nextTop();
    } while((S1_height > 0) && (!success));
    if (S1_height == 0)
        return false;
    else
        return S1_first();
}

bool KeccakFTrailWithGivenParityIterator::S2_push(unsigned int y)
{
    bool canAfford = pushBitInUnaffectedOddColumn(UOcolumns[S2_height], y);
    if (canAfford) {
        S3_yMin[UOcolumns[S2_height].getXplus5Z()] = y+1; // orbitals start after bit set by stack 2
        S2_y.push(y);
        S2_height++;
        return true;
    }
    else
        return false;
}

unsigned int KeccakFTrailWithGivenParityIterator::S2_pop()
{
    pop();
    S2_height--;
    unsigned int y = S2_y.top();
    S2_y.pop();
    return y;
}

bool KeccakFTrailWithGivenParityIterator::S2_firstTop()
{
    unsigned int y = 0;
    bool success;
    do {
        success = S2_push(y);
        y++;
    } while((!success) && (y < 5));
    return success;
}

bool KeccakFTrailWithGivenParityIterator::S2_nextTop()
{
    unsigned int y = S2_pop();
    y++;
    while(y < 5) {
        if (S2_push(y))
            return true;
        y++;
    }
    return false;
}

bool KeccakFTrailWithGivenParityIterator::S2_first()
{
    while(S2_height < UOcolumns.size()) {
        bool success = true;
        while((S2_height < UOcolumns.size()) && success) {
            success = S2_firstTop();
        }
        while((S2_height > 0) && (!success)) {
            success = S2_nextTop();
        }
        if (S2_height == 0)
            return false;
    }
    return true;
}

bool KeccakFTrailWithGivenParityIterator::S2_next()
{
    bool success = true;
    do {
        success = S2_nextTop();
    } while((S2_height > 0) && (!success));
    if (S2_height == 0)
        return false;
    else
        return S2_first();
}

bool KeccakFTrailWithGivenParityIterator::S3_push(const OrbitalPosition& orbital)
{
    bool canAfford = pushOrbitalInUnaffectedColumn(orbital);
    if (canAfford) {
        S3_position.push(orbital);
        S3_height++;
        return true;
    }
    else
        return false;
}

bool KeccakFTrailWithGivenParityIterator::S3_addNewOrbital()
{
    OrbitalPosition orbital;
    if (S3_height == 0) {
        if (!orbital.first(S3_yMin, laneSize))
            return false;
    }
    else {
        if (!orbital.successorOf(S3_position.top(), S3_yMin, laneSize))
            return false;
    }
    bool success, notLast = true;
    do {
        success = S3_push(orbital);
        if (!success)
            notLast = orbital.next(S3_yMin, laneSize);
    } while((!success) && notLast);
    return success;
}

bool KeccakFTrailWithGivenParityIterator::S3_nextTop()
{
    OrbitalPosition orbital = S3_position.top();
    pop();
    S3_position.pop();
    S3_height--;
    while(orbital.next(S3_yMin, laneSize)) {
        if (S3_push(orbital))
            return true;
    }
    return false;
}

bool KeccakFTrailWithGivenParityIterator::S3_next()
{
    if (S3_addNewOrbital())
        return true;
    else {
        bool success = false;
         while((!success) && (S3_height > 0))
            success = S3_nextTop();
        return success;
    }
}

bool KeccakFTrailWithGivenParityIterator::first()
{
    if (Acolumns.size() > 0) {
        if (!S1_first())
            return false;
        if (UOcolumns.size() > 0) {
            while (!S2_first()) {
                if (!S1_next())
                    return false;
            }
        }
    }
    return true;
}

bool KeccakFTrailWithGivenParityIterator::next()
{
    if ((!orbitals) || (!S3_next())) {
        if (UOcolumns.size() > 0) {
            if (S2_next())
                return true;
            else {
                if (Acolumns.size() > 0) {
                    bool success;
                    do {
                        if (!S1_next())
                            return false;
                        success = S2_first();
                    } while (!success);
                    return true;
                }
                else
                    return false;
            }
        }
        else {
            if (Acolumns.size() > 0)
                return S1_next();
            else
                return false;
        }
    }
    else
        return true;
}

bool KeccakFTrailWithGivenParityIterator::isEnd()
{
    if (!initialized) initialize();
    return end;
}

bool KeccakFTrailWithGivenParityIterator::isEmpty()
{
    if (!initialized) initialize();
    return empty;
}

bool KeccakFTrailWithGivenParityIterator::isBounded()
{
    if (!initialized) initialize();
    return false;
}

UINT64 KeccakFTrailWithGivenParityIterator::getIndex()
{
    if (!initialized) initialize();
    return index;
}

UINT64 KeccakFTrailWithGivenParityIterator::getCount()
{
    if (!initialized) initialize();
    return (end ? index : 0);
}

void KeccakFTrailWithGivenParityIterator::operator++()
{
    if (!initialized) initialize();
    if (!end) {
        ++index;
        if (next())
            getTrail();
        else
            end = true;
    }
}

const Trail& KeccakFTrailWithGivenParityIterator::operator*()
{
    if (!initialized) initialize();
    return trail;
}

KeccakFTwoRoundTrailCoreWithGivenParityIterator::KeccakFTwoRoundTrailCoreWithGivenParityIterator(const KeccakFPropagation& aDCorLC,
        const vector<RowValue>& aParity, int aMaxWeight, bool aOrbitals)
    : KeccakFTrailWithGivenParityIterator(aDCorLC,aParity, aOrbitals), maxWeight(aMaxWeight)
{
    stack_stateAtA.push(StateAsVectorOfSlices(laneSize, 0));
    stack_stateAtB.push(StateAsVectorOfSlices(laneSize, 0));
    stack_weight.push(0);
}

int KeccakFTwoRoundTrailCoreWithGivenParityIterator::setBitToOneAndGetDeltaWeight(vector<SliceValue>& state, const BitPosition& p) const
{
    int weightBefore = DCorLC.getWeight(state[p.z]);
    setBitToOne(state, p.x, p.y, p.z);
    return (int)DCorLC.getWeight(state[p.z]) - weightBefore;
}

int KeccakFTwoRoundTrailCoreWithGivenParityIterator::setBitToOneAndGetDeltaMinReverseWeight(vector<SliceValue>& state, const BitPosition& p) const
{
    int weightBefore = DCorLC.getMinReverseWeight(state[p.z]);
    setBitToOne(state, p.x, p.y, p.z);
    return (int)DCorLC.getMinReverseWeight(state[p.z]) - weightBefore;
}

bool KeccakFTwoRoundTrailCoreWithGivenParityIterator::pushValueInAffectedColumn(const ColumnPosition& columnBeforeTheta, ColumnValue valueBeforeTheta)
{
    stack_stateAtA.push(stack_stateAtA.top());
    stack_stateAtB.push(stack_stateAtB.top());
    int deltaWeight = setValueInAffectedColumnAndGetDeltaTotalWeight(stack_stateAtA.top(), stack_stateAtB.top(), columnBeforeTheta, valueBeforeTheta);
    int newWeight = deltaWeight + stack_weight.top();
    if (newWeight <= maxWeight) {
        stack_weight.push(newWeight);
        return true;
    }
    else {
        stack_stateAtA.pop();
        stack_stateAtB.pop();
        return false;
    }
}

bool KeccakFTwoRoundTrailCoreWithGivenParityIterator::pushBitInUnaffectedOddColumn(const ColumnPosition& columnBeforeTheta, unsigned int y)
{
    stack_stateAtA.push(stack_stateAtA.top());
    stack_stateAtB.push(stack_stateAtB.top());
    int deltaWeight = setBitInUnaffectedColumnAndGetDeltaTotalWeight(stack_stateAtA.top(), stack_stateAtB.top(), columnBeforeTheta, y);
    int newWeight = deltaWeight + stack_weight.top();
    if (newWeight <= maxWeight) {
        stack_weight.push(newWeight);
        return true;
    }
    else {
        stack_stateAtA.pop();
        stack_stateAtB.pop();
        return false;
    }
}

bool KeccakFTwoRoundTrailCoreWithGivenParityIterator::pushOrbitalInUnaffectedColumn(const OrbitalPosition& orbital)
{
    stack_stateAtA.push(stack_stateAtA.top());
    stack_stateAtB.push(stack_stateAtB.top());
    int deltaWeight = setOrbitalInUnaffectedColumnAndGetDeltaTotalWeight(stack_stateAtA.top(), stack_stateAtB.top(), orbital);
    int newWeight = deltaWeight + stack_weight.top();
    if (newWeight <= maxWeight) {
        stack_weight.push(newWeight);
        return true;
    }
    else {
        stack_stateAtA.pop();
        stack_stateAtB.pop();
        return false;
    }
}

void KeccakFTwoRoundTrailCoreWithGivenParityIterator::pop()
{
    stack_stateAtA.pop();
    stack_stateAtB.pop();
    stack_weight.pop();
}

void KeccakFTwoRoundTrailCoreWithGivenParityIterator::getTrail()
{
    trail.clear();
    trail.setFirstStateReverseMinimumWeight(DCorLC.getMinReverseWeight(stack_stateAtA.top()));
    trail.append(stack_stateAtB.top(), DCorLC.getWeight(stack_stateAtB.top()));
}
