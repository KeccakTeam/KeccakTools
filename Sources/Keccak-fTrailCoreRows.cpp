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

#include <sstream>
#include "Keccak-fTrailCoreRows.h"

KeccakFTrailCoreRows::KeccakFTrailCoreRows(const KeccakFDCLC& aParent, KeccakFPropagation::DCorLC aDCorLC)
    : KeccakFPropagation(aParent, aDCorLC)
{
}

void KeccakFTrailCoreRows::generateTrailCoresBasedOnRows(TrailFetcher& trailsOut, int maxNrRowsAtA, int maxNrRowsAtB, int maxWeight)
{
    if ((maxNrRowsAtA > 3) && (maxNrRowsAtB > 3))
        throw KeccakException("This method generates up to 3 rows only.");
    generateTrailCoresBasedOnRows(trailsOut, (maxNrRowsAtA < maxNrRowsAtB), maxNrRowsAtA, maxNrRowsAtB, maxWeight);
}

void KeccakFTrailCoreRows::generateTrailCoresUpToGivenWeight(TrailFetcher& trailsOut, int maxMinRevWeightAtA, int maxWeightAtB, int maxWeight)
{
    if ((maxMinRevWeightAtA > 7) && (maxWeightAtB > 7))
        throw KeccakException("This method generates up to 3 rows only.");
    generateTrailCoresUpToGivenWeight(trailsOut, (maxMinRevWeightAtA < maxWeightAtB), maxMinRevWeightAtA, maxWeightAtB, maxWeight);
}

void KeccakFTrailCoreRows::generateTrailCoresBasedOnRows(TrailFetcher& trailsOut, bool startingFromA, int maxNrRowsAtA, int maxNrRowsAtB, int maxWeight)
{
    int maxNrRows = (startingFromA ? maxNrRowsAtA : maxNrRowsAtB);
    if (maxNrRows >= 1) {
        const unsigned int z = 0;
        for(unsigned int y=0; y<5; y++)
        for(RowValue a=1; a<32; a++) {
            vector<SliceValue> in(laneSize, 0);
            in[z] = getSliceFromRow(a, y);
            filterGeneratedTrailCores(trailsOut, in, startingFromA, maxNrRowsAtA, maxNrRowsAtB, maxWeight);
        }
    }
    if (maxNrRows >= 2) {
        progress.stack("Generating 2 rows");
        const unsigned int z1 = 0;
        for(unsigned int y1=0; y1<5; y1++)
        for(unsigned int z2=0; z2<=laneSize/2; z2++)
        for(unsigned int y2=0; y2<5; y2++)
        if ((z1 != z2) || (y1 < y2)) {
            for(RowValue a1=1; a1<32; a1++)
            for(RowValue a2=1; a2<32; a2++) {
                vector<SliceValue> in(laneSize, 0);
                in[z1] ^= getSliceFromRow(a1, y1);
                in[z2] ^= getSliceFromRow(a2, y2);
                filterGeneratedTrailCores(trailsOut, in, startingFromA, maxNrRowsAtA, maxNrRowsAtB, maxWeight);
            }
            ++progress;
        }
        progress.unstack();
    }
    if (maxNrRows >= 3) {
        progress.stack("Generating 3 rows");
        const unsigned int z1 = 0;
        for(unsigned int z2=0; z2<laneSize; z2++)
        for(unsigned int z3=z2; z3<laneSize; z3++) {
            LaneValue test1 = (LaneValue)1 | ((LaneValue)1 << z2) | ((LaneValue)1 << z3);
            LaneValue test2 = test1;
            parent.ROL(test2, -(int)z2);
            LaneValue test3 = test1;
            parent.ROL(test3, -(int)z3);
            if ((test1 <= test2) && (test1 <= test3)) {
                {
                    stringstream str;
                    str << "Rows in slices " << dec << z1 << ", " << z2 << ", " << z3;
                    progress.stack(str.str());
                }
                for(unsigned int y1=0; y1<5; y1++)
                for(unsigned int y2=0; y2<5; y2++)
                for(unsigned int y3=0; y3<5; y3++)
                if (((z1 != z2) || (y1 < y2)) && ((z2 != z3) || (y2 < y3))) {
                    for(RowValue a1=1; a1<32; a1++)
                    for(RowValue a2=1; a2<32; a2++)
                    for(RowValue a3=1; a3<32; a3++) {
                        vector<SliceValue> in(laneSize, 0);
                        in[z1] ^= getSliceFromRow(a1, y1);
                        in[z2] ^= getSliceFromRow(a2, y2);
                        in[z3] ^= getSliceFromRow(a3, y3);
                        filterGeneratedTrailCores(trailsOut, in, startingFromA, maxNrRowsAtA, maxNrRowsAtB, maxWeight);
                    }
                    ++progress;
                }
                progress.unstack();
                ++progress;
            }
        }
        progress.unstack();
    }
}

void KeccakFTrailCoreRows::generateTrailCoresUpToGivenWeight(TrailFetcher& trailsOut, bool startingFromA, int maxMinRevWeightAtA, int maxWeightAtB, int maxWeight)
{
    int maxMinRevWeightAtAorB = (startingFromA ? maxMinRevWeightAtA : maxWeightAtB);
    int maxNrRows = maxMinRevWeightAtAorB/2;
    if (maxNrRows >= 1) {
        const unsigned int z = 0;
        for(RowValue a=1; a<32; a++) {
            int weight = (startingFromA ? getMinReverseWeightRow(a) : getWeightRow(a));
            if (weight <= maxMinRevWeightAtAorB) for(unsigned int y=0; y<5; y++) {
                vector<SliceValue> in(laneSize, 0);
                in[z] = getSliceFromRow(a, y);
                filterGeneratedTrailCoresUpToGivenWeight(trailsOut, in, startingFromA, maxMinRevWeightAtA, maxWeightAtB, maxWeight);
            }
        }
    }
    if (maxNrRows >= 2) {
        progress.stack("Generating 2 rows");
        for(RowValue a1=1; a1<32; a1++)
        for(RowValue a2=1; a2<32; a2++) {
            int weight = (startingFromA ? 
                getMinReverseWeightRow(a1)+getMinReverseWeightRow(a2) 
              : getWeightRow(a1)+getWeightRow(a2));
            if (weight <= maxMinRevWeightAtAorB) {
                {
                    stringstream str;
                    str << "Row values " << hex << (int)a1 << "," << (int)a2 
                        << " of " << (startingFromA ? "min. rev. " : "") << "weight " 
                        << dec << weight;
                    progress.stack(str.str());
                }
                const unsigned int z1 = 0;
                for(unsigned int y1=0; y1<5; y1++)
                for(unsigned int z2=0; z2<=laneSize/2; z2++)
                for(unsigned int y2=0; y2<5; y2++)
                if ((z1 != z2) || (y1 < y2)) {
                    vector<SliceValue> in(laneSize, 0);
                    in[z1] ^= getSliceFromRow(a1, y1);
                    in[z2] ^= getSliceFromRow(a2, y2);
                    filterGeneratedTrailCoresUpToGivenWeight(trailsOut, in, startingFromA, maxMinRevWeightAtA, maxWeightAtB, maxWeight);
                    ++progress;
                }
                progress.unstack();
            }
            ++progress;
        }
        progress.unstack();
    }
    if (maxNrRows >= 3) {
        progress.stack("Generating 3 rows");
        for(RowValue a1=1; a1<32; a1++)
        for(RowValue a2=1; a2<32; a2++)
        for(RowValue a3=1; a3<32; a3++) {
            int weight = (startingFromA ? 
                getMinReverseWeightRow(a1)+getMinReverseWeightRow(a2)+getMinReverseWeightRow(a3)  
              : getWeightRow(a1)+getWeightRow(a2)+getWeightRow(a3));
            if (weight <= maxMinRevWeightAtAorB) {
                {
                    stringstream str;
                    str << "Row values " << hex << (int)a1 << "," << (int)a2 << "," << (int)a3
                        << " of " << (startingFromA ? "min. rev. " : "") << "weight " 
                        << dec << weight;
                    progress.stack(str.str());
                }
                const unsigned int z1 = 0;
                for(unsigned int z2=0; z2<laneSize; z2++)
                for(unsigned int z3=z2; z3<laneSize; z3++) {
                    LaneValue test1 = (LaneValue)1 | ((LaneValue)1 << z2) | ((LaneValue)1 << z3);
                    LaneValue test2 = test1;
                    parent.ROL(test2, -(int)z2);
                    LaneValue test3 = test1;
                    parent.ROL(test3, -(int)z3);
                    if ((test1 <= test2) && (test1 <= test3)) {
                        for(unsigned int y1=0; y1<5; y1++)
                        for(unsigned int y2=0; y2<5; y2++)
                        for(unsigned int y3=0; y3<5; y3++)
                        if (((z1 != z2) || (y1 < y2)) && ((z2 != z3) || (y2 < y3))) {
                            vector<SliceValue> in(laneSize, 0);
                            in[z1] ^= getSliceFromRow(a1, y1);
                            in[z2] ^= getSliceFromRow(a2, y2);
                            in[z3] ^= getSliceFromRow(a3, y3);
                            filterGeneratedTrailCoresUpToGivenWeight(trailsOut, in, startingFromA, maxMinRevWeightAtA, maxWeightAtB, maxWeight);
                            ++progress;
                        }
                    }
                }
                progress.unstack();
            }
            ++progress;
        }
        progress.unstack();
    }
}

void KeccakFTrailCoreRows::filterGeneratedTrailCores(TrailFetcher& trailsOut, const vector<SliceValue>& stateAtAorB, bool stateAtA, int maxNrRowsAtA, int maxNrRowsAtB, int maxWeight)
{
    vector<SliceValue> stateAtBorA;
    if (stateAtA) {
        directLambda(stateAtAorB, stateAtBorA);
        if (getNrActiveRows(stateAtBorA) > maxNrRowsAtB)
            return;
        if ((getMinReverseWeight(stateAtAorB) + getWeight(stateAtBorA)) > maxWeight)
            return;
    }
    else {
        reverseLambda(stateAtAorB, stateAtBorA);
        if (getNrActiveRows(stateAtBorA) > maxNrRowsAtA)
            return;
        if ((getMinReverseWeight(stateAtBorA) + getWeight(stateAtAorB)) > maxWeight)
            return;
    }
    Trail trail;
    if (stateAtA) {
        trail.setFirstStateReverseMinimumWeight(getMinReverseWeight(stateAtAorB));
        trail.append(stateAtBorA, getWeight(stateAtBorA));
    }
    else {
        trail.setFirstStateReverseMinimumWeight(getMinReverseWeight(stateAtBorA));
        trail.append(stateAtAorB, getWeight(stateAtAorB));
    }
    trailsOut.fetchTrail(trail);
}

void KeccakFTrailCoreRows::filterGeneratedTrailCoresUpToGivenWeight(TrailFetcher& trailsOut, const vector<SliceValue>& stateAtAorB, bool stateAtA, int maxMinRevWeightAtA, int maxWeightAtB, int maxWeight)
{
    vector<SliceValue> stateAtBorA;
    if (stateAtA) {
        directLambda(stateAtAorB, stateAtBorA);
        if (getWeight(stateAtBorA) > maxWeightAtB)
            return;
        if ((getMinReverseWeight(stateAtAorB) + getWeight(stateAtBorA)) > maxWeight)
            return;
    }
    else {
        reverseLambda(stateAtAorB, stateAtBorA);
        if (getMinReverseWeight(stateAtBorA) > maxMinRevWeightAtA)
            return;
        if ((getMinReverseWeight(stateAtBorA) + getWeight(stateAtAorB)) > maxWeight)
            return;
    }
    Trail trail;
    if (stateAtA) {
        trail.setFirstStateReverseMinimumWeight(getMinReverseWeight(stateAtAorB));
        trail.append(stateAtBorA, getWeight(stateAtBorA));
    }
    else {
        trail.setFirstStateReverseMinimumWeight(getMinReverseWeight(stateAtBorA));
        trail.append(stateAtAorB, getWeight(stateAtAorB));
    }
    trailsOut.fetchTrail(trail);
}
