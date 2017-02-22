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

#include "Keccak-fDisplay.h"
#include "Keccak-fTrailCoreInKernelAtC.h"
#include "translationsymmetry.h"

using namespace std;

TrailCoreInKernelAtC::TrailCoreInKernelAtC(const vector<SliceValue>& backgroundAtA,
                                           const vector<SliceValue>& aTabooAtB,
                                           unsigned int aMaxWeight,
                                           const KeccakFDCLC& aParent,
                                           KeccakFPropagation::DCorLC aDCorLC) :
    TrailCore3Rounds(backgroundAtA, aTabooAtB, aMaxWeight, aParent, aDCorLC)
{
    partialStateAtD.assign(laneSize,0);

    weightAtB = 0;

    knotPointDeficit = 0;
    knotWeightAtBDeficit = 0;

    partialHammingWeightAtD = 0;
    partialNrActiveRowsAtD = 0;

    populateStatesWithBackground(backgroundAtA);

    nrActiveRowsAtA  = getNrActiveRows(stateAtA);
    hammingWeightAtA = getHammingWeight(stateAtA);

    if (weightAtB == 0) {
        CoreInfo workCoreInfo;

        workCoreInfo.hammingWeightAtA = 0;
        workCoreInfo.hammingWeightAtA = 0;
        workCoreInfo.nrActiveRowsAtA = 0;

        workCoreInfo.stateAtB.assign(laneSize,0);
        workCoreInfo.weightAtB = 0;

        workCoreInfo.partialStateAtC.assign(laneSize,0);
        workCoreInfo.hammingWeightAtD = 0;
        workCoreInfo.nrActiveRowsAtD = 0;
        workCoreInfo.partialWeight = 0;

        workCoreInfo.vortexLength = 0;
        workCoreInfo.vortexIndex = 0;
        workCoreInfo.vortexZOffset = 0;

        outCore.push_back(workCoreInfo);
    }
}

bool TrailCoreInKernelAtC::isStateAtBWellFormed() const
{
    return (0 == knotPointDeficit);
}

bool TrailCoreInKernelAtC::canAffordGeneric(unsigned int deltaNrKnotPointsWorkingChain,
                                            unsigned int additionalKnotOrRun,
                                            unsigned int nrOrbitalPointsPerDeltaChain,
                                            unsigned int deltaNrOrbitalPointsWorkingChain) const
{
    // Each knot needs to have at least 2 (resp. 3) points in the case of a knot with (resp. without) background. If not the case yet, this implies that chains will have to be added. As each chain contributes 2 knot points, the deficit in chains is ceil(deficit in knot points / 2).
    // The projected deficit in knot points is the current knot point deficit plus 2 per projected additional knot, minus the knot points brought by the working chain (deltaNrKnotPointsWorkingChain).
    int chainDeficit = knotPointDeficit + 2*additionalKnotOrRun; // as a new knot or run has at least knotPointDeficit 2

    chainDeficit = max(0, (int)chainDeficit-(int)deltaNrKnotPointsWorkingChain);
    chainDeficit = (chainDeficit + 1)/2;

    // The number of orbital points that still have to be added has two contributions. The first one is the orbitals that have yet to be added to the working chain. And the second one is the number of orbitals for each chain that still has to be added.
    unsigned int orbitalPointDeficit = nrOrbitalPointsPerDeltaChain*chainDeficit + deltaNrOrbitalPointsWorkingChain;

    // The projected Hamming weight at A is augmented with the number of orbital points and knot points still to be added.
    unsigned int projectedHammingWeightAtA = hammingWeightAtA + orbitalPointDeficit + knotPointDeficit;

    // The projected weight at B is augmented with the orbital points (each contributing a weight 2 because in different rows) and the knot points still to be added.
    unsigned int projectedWeightAtB = weightAtB + (2*orbitalPointDeficit) + knotWeightAtBDeficit;

    // The projected Hamming weight at A is augmented with the number of orbital points still to be added. The knot points do not necessarily survive through χ, hence they are not counted.
    unsigned int projectedPartialHammingWeightAtD = partialHammingWeightAtD + orbitalPointDeficit;

    unsigned int lowerWeight = getLowerBoundOnReverseWeightGivenHammingWeightAndNrActiveRows(projectedHammingWeightAtA, nrActiveRowsAtA)
                                + projectedWeightAtB
                                + getLowerBoundOnWeightGivenHammingWeightAndNrActiveRows(projectedPartialHammingWeightAtD, partialNrActiveRowsAtD);

    return (lowerWeight <= maxWeight);
}

void TrailCoreInKernelAtC::addPoint(const BitPosition& pB, bool toKnotSlice, bool isBackgroundPoint)
{
    if (!isBackgroundPoint){
        yOffsets.back().push_back(0);
        chains.back().push_back(pB);

        // Dealing with the impact at A
        BitPosition pA(pB);
        reverseRhoPi(pA);
        hammingWeightAtA += 1;
        if (getRow(stateAtA,RowPosition(pA)) == 0) nrActiveRowsAtA += 1;
        setBitToOne(stateAtA, pA);
    }
    else knotsWithBackground.insert(pB.z);
    // Note: a background point cannot be removed with removePoint().

    // Dealing with the impact at B and D
    if (toKnotSlice) {
        weightAtB -= getWeight(stateAtB[pB.z]);
        setBitToOne(stateAtB, pB);
        weightAtB += getWeight(stateAtB[pB.z]);

        KnotInformation newKnotInfo;
        bool hasBackground = (knotsWithBackground.count(pB.z) != 0);
        if (knots.count(pB.z) == 0) {
            knotPointAddedKnot.push(true);
            populateKnotInfo(newKnotInfo, stateAtB[pB.z], true, hasBackground);
            if (hasBackground) {
                knotPointDeficit        += 1;
                knotWeightAtBDeficit    += 2;
                partialHammingWeightAtD += 1;
            }
            else {
                knotPointDeficit        += 2;
                knotWeightAtBDeficit    += 3;
                partialHammingWeightAtD += 1;
            }
        }
        else {
            knotPointAddedKnot.push(false);
            KnotInformation oldKnotInfo = knots[pB.z];
            populateKnotInfo(newKnotInfo,stateAtB[pB.z], false, hasBackground);
            knotPointDeficit        += newKnotInfo.knotPointDeficit     - oldKnotInfo.knotPointDeficit;
            knotWeightAtBDeficit    += newKnotInfo.knotWeightAtBDeficit - oldKnotInfo.knotWeightAtBDeficit;
            partialHammingWeightAtD += newKnotInfo.nrActiveRows         - oldKnotInfo.nrActiveRows;
        }
        // partialStateAtD is not updated because no bit is for sure active at C even after adding more points to the knot.
        // Same remark for nrActiveRows.
        // partialHammingWeightAtD is for sure at least one bit per active row.
        knots[pB.z] = newKnotInfo;
    }
    else {
        // An orbital point is alone in a row at B, hence has weight 2.
        weightAtB += 2;
        setBitToOne(stateAtB, pB);

        // Dealing with the impact at D
        BitPosition pD(pB);
        directRhoPi(pD);
        partialHammingWeightAtD += 1;
        if (getRow(partialStateAtD, RowPosition(pD)) == 0) partialNrActiveRowsAtD += 1;
        setBitToOne(partialStateAtD, pD);
    }
}

void TrailCoreInKernelAtC::removePoint(bool fromKnotSlice)
{
    BitPosition pB(chains.back().back());
    yOffsets.back().pop_back();
    chains.back().pop_back();

    BitPosition pA(pB);
    reverseRhoPi(pA);
    setBitToZero(stateAtA, pA);
    hammingWeightAtA -= 1;
    if (getRow(stateAtA, RowPosition(pA)) == 0) nrActiveRowsAtA -= 1;

    if (fromKnotSlice) {
        if (knotPointAddedKnot.top()) {
            setBitToZero(stateAtB, pB);
            weightAtB               -= 2;
            knotPointDeficit        -= 2;
            knotWeightAtBDeficit    -= 3;
            partialHammingWeightAtD -= 1;
            knots.erase(pB.z);
        }
        else {
            weightAtB -= getWeight(stateAtB[pB.z]);
            setBitToZero(stateAtB, pB);
            weightAtB += getWeight(stateAtB[pB.z]);

            KnotInformation oldKnotInfo = knots[pB.z];
            KnotInformation updatedKnotInfo;
            bool hasBackground = (0 != knotsWithBackground.count(pB.z));
            bool knotHasSinglePoint = (getHammingWeightSlice(stateAtB[pB.z]) == 1);
            populateKnotInfo(updatedKnotInfo,stateAtB[pB.z],knotHasSinglePoint,hasBackground);
            knotPointDeficit        += updatedKnotInfo.knotPointDeficit     - oldKnotInfo.knotPointDeficit;
            knotWeightAtBDeficit    += updatedKnotInfo.knotWeightAtBDeficit - oldKnotInfo.knotWeightAtBDeficit;
            partialHammingWeightAtD += updatedKnotInfo.nrActiveRows         - oldKnotInfo.nrActiveRows;
            knots[pB.z] = updatedKnotInfo;
        }
        knotPointAddedKnot.pop();
    }
    else {
        setBitToZero(stateAtB, pB);
        weightAtB    -= 2;
        BitPosition pD(pB);
        directRhoPi(pD);
        setBitToZero(partialStateAtD, pD);
        partialHammingWeightAtD -= 1;
        if (getRow(partialStateAtD, RowPosition(pD)) == 0) partialNrActiveRowsAtD -= 1;
    }
}

void TrailCoreInKernelAtC::convertKnotPointToOrbitalPoint()
{
    // The point pB is assumed to be the only point in its knot.
    BitPosition pB = chains.back().back();

    knotPointAddedKnot.pop(); // Popped because addPoint(toKnotSlice=false) does not push onto knotPointAddedKnot.
    knots.erase(pB.z);
    knotPointDeficit     -= 2;
    knotWeightAtBDeficit -= 3;

    BitPosition pD(pB);
    directRhoPi(pD);

    // A knot point does not count for partialStateAtD or partialNrActiveRowsAtD, but an orbital point does.
    if (getRow(partialStateAtD, RowPosition(pD)) == 0) partialNrActiveRowsAtD += 1;
    setBitToOne(partialStateAtD, pD);
}

bool TrailCoreInKernelAtC::mayBeEndPoint(const BitPosition& pB)
{
    (void)pB;
    // This method is called when a chain arrives in a slice that is not yet a knot. So the question to be answered is:
    // does the current weight allow adding a knot?
    // As the chain is not yet complete, there is still one knot point to be added. Assuming the knot point introduces a new knot, it will create a new knot.
    return canAffordGeneric(1, 1, chains.back().size() - 1, 0);
}

bool TrailCoreInKernelAtC::mayBeStartPointSliceAndGoThere(unsigned int& z, bool zIsInitialized)
{
    map<unsigned int,KnotInformation>::const_iterator it;
    if (zIsInitialized)
        it = knots.upper_bound(z);
    else
        it = knots.begin();
    if (it == knots.end()) return false;
    z = it->first;
    return true;
}

unsigned int TrailCoreInKernelAtC::computeLowerWeightAssumingVortexIsAdded()
{
    const VortexInfo& v = vortexBase[outCore.back().vortexLength/2][outCore.back().vortexIndex];

    unsigned int localNrActiveRowsAtA = max(outCore.back().nrActiveRowsAtA, v.nrActiveRowsAtA);
    unsigned int localLowerWeight = getLowerBoundOnReverseWeightGivenHammingWeightAndNrActiveRows(
        outCore.back().hammingWeightAtA + outCore.back().vortexLength,
        localNrActiveRowsAtA);

    localLowerWeight += outCore.back().weightAtB + 2*outCore.back().vortexLength;

    unsigned int localNrActiveRowsAtD = max(outCore.back().nrActiveRowsAtD,v.nrActiveRowsAtD);
    localLowerWeight += getLowerBoundOnWeightGivenHammingWeightAndNrActiveRows(
        outCore.back().hammingWeightAtD + outCore.back().vortexLength,
        localNrActiveRowsAtD);

    return localLowerWeight;
}

bool TrailCoreInKernelAtC::next()
{
    do {
        if (outCore.empty()){
            if (!nextWithKnots()) return false;
            CoreInfo workCoreInfo;

            workCoreInfo.hammingWeightAtA = hammingWeightAtA;
            workCoreInfo.hammingWeightAtA = hammingWeightAtA;
            workCoreInfo.nrActiveRowsAtA = nrActiveRowsAtA;

            workCoreInfo.stateAtB = stateAtB;
            workCoreInfo.weightAtB = weightAtB;

            workCoreInfo.partialStateAtC = stateAtB;
            map<unsigned int,KnotInformation>::const_iterator it = knots.begin();
            while(it != knots.end()) {
                workCoreInfo.partialStateAtC[it->first] = getMinimumInKernelSliceAfterChi(workCoreInfo.partialStateAtC[it->first]);
                it++;
            }
            vector<SliceValue> localStateAtD(laneSize);
            directLambdaAfterTheta(workCoreInfo.partialStateAtC,localStateAtD);
            workCoreInfo.hammingWeightAtD = getHammingWeight(localStateAtD);
            workCoreInfo.nrActiveRowsAtD = getNrActiveRows(localStateAtD);
            workCoreInfo.partialWeight = getMinReverseWeight(stateAtA) + workCoreInfo.weightAtB + getWeight(localStateAtD);

            workCoreInfo.vortexLength = 0;
            workCoreInfo.vortexIndex = 0;
            workCoreInfo.vortexZOffset = 0;
            outCore.push_back(workCoreInfo);
            if (outCore.back().partialWeight <= maxWeight)
                return true;
            else
                outCore.pop_back();
        }
        else {
            bool foundGoodVortexToAdd = true;
            if ((!knots.empty() || (outCore.size() > 1)) && // In absence of knots, the first vortex has a fixed position: vortexZOffset = 0
                    (outCore.back().vortexIndex < vortexBase[outCore.back().vortexLength/2].size()) &&  // vortexIndex must point to an existing entry
                    (outCore.back().vortexZOffset < laneSize-1)) {
                outCore.back().vortexZOffset += 1;
                if ((outCore.size() > 1) &&
                        (outCore[0].vortexLength == outCore.back().vortexLength) &&
                        (outCore[0].vortexIndex  == outCore.back().vortexIndex ))   {
                    vector<unsigned int> zPattern(laneSize, 0);
                    for (unsigned int i=0 ; i<outCore.size() ; i++)
                        zPattern[outCore[i].vortexZOffset] = 1;
                    foundGoodVortexToAdd = isMinimalSymmetrically(zPattern);
                }
            }
            else if ((int)outCore.back().vortexIndex < (int)vortexBase[outCore.back().vortexLength/2].size()-1) {
                outCore.back().vortexZOffset = 0;
                outCore.back().vortexIndex += 1;
            }
            else {
                outCore.back().vortexZOffset = 0;
                outCore.back().vortexIndex = 0;
                outCore.back().vortexLength += 2;
                if (outCore.back().vortexLength/2 >= vortexBase.size()) {
                    outCore.pop_back();
                    foundGoodVortexToAdd = false;
                }
                else if (vortexBase[outCore.back().vortexLength/2].empty())
                    foundGoodVortexToAdd = false;
                else if (outCore.back().partialWeight + 2*outCore.back().vortexLength > maxWeight) {
                    outCore.pop_back();
                    foundGoodVortexToAdd = false;
                }
            }

            if (foundGoodVortexToAdd) {
                foundGoodVortexToAdd = foundGoodVortexToAdd && (computeLowerWeightAssumingVortexIsAdded() <= maxWeight);
                const VortexInfo& v = vortexBase[outCore.back().vortexLength/2][outCore.back().vortexIndex];
                if (foundGoodVortexToAdd) { // Now test the vortex to add for overlap with the state up to now and its tabooAtB
                    map<unsigned int,SliceValue>::const_iterator it = v.stateAtB.slices.begin();
                    while (foundGoodVortexToAdd && (it != v.stateAtB.slices.end())){
                        unsigned int localZ = ((it->first)+outCore.back().vortexZOffset)%laneSize;
                        foundGoodVortexToAdd = (outCore.back().stateAtB[localZ] == 0);
                        foundGoodVortexToAdd = foundGoodVortexToAdd  && (((tabooAtB[localZ])&(it->second)) == 0);
                        it++;
                    }
                }
                if (foundGoodVortexToAdd) { // Now really adding the vortex
                    outCore.push_back(outCore.back());
                    map<unsigned int,SliceValue>::const_iterator it = v.stateAtB.slices.begin();
                    while (it != v.stateAtB.slices.end()){
                        unsigned int localZ = ((it->first)+outCore.back().vortexZOffset)%laneSize;
                        outCore.back().stateAtB[localZ] = it->second;
                        outCore.back().partialStateAtC[localZ] = it->second;
                        it++;
                    }
                    outCore.back().weightAtB += 2*outCore.back().vortexLength;

                    vector<SliceValue> localStateAtA(laneSize);
                    reverseLambda(outCore.back().stateAtB,localStateAtA);
                    outCore.back().hammingWeightAtA = getHammingWeight(localStateAtA);
                    outCore.back().nrActiveRowsAtA = getNrActiveRows(localStateAtA);

                    vector<SliceValue> localStateAtD(laneSize);
                    directLambdaAfterTheta(outCore.back().partialStateAtC,localStateAtD);

                    outCore.back().hammingWeightAtD = getHammingWeight(localStateAtD);
                    outCore.back().nrActiveRowsAtD = getNrActiveRows(localStateAtD);
                    outCore.back().partialWeight = getMinReverseWeight(localStateAtA) + outCore.back().weightAtB + getWeight(localStateAtD);
                    if (outCore.back().partialWeight > maxWeight)
                        outCore.pop_back();
                    else
                        return true;
                }
            }
        }
    }
    while(true);
}

const TrailCoreInKernelAtC::CoreInfo& TrailCoreInKernelAtC::getTopCoreInfo() const
{
    return outCore.back();
}

ostream& operator<<(ostream& fout, const TrailCoreInKernelAtC& aL)
{
    for(unsigned int chainNr = 0 ; chainNr < aL.chains.size() ; chainNr++) {
        fout << "chain " << chainNr << ":  ";
        for(unsigned int pointNr = 0 ; pointNr < aL.chains[chainNr].size() ; pointNr++)
            fout << aL.chains[chainNr][pointNr] << aL.yOffsets[chainNr][pointNr] << " ";
        fout << endl;
    }
    for (unsigned int vortexNr = 0 ; vortexNr < aL.outCore.size()-1 ; vortexNr++) {
        fout << "vortex " << vortexNr << ": ";
        fout << "length " << aL.outCore[vortexNr].vortexLength;
        fout << " index " << aL.outCore[vortexNr].vortexIndex;
        fout << " offset " << aL.outCore[vortexNr].vortexZOffset;
        fout << endl;
    }
    fout << "state at B" << endl;
    displayState(fout,aL.outCore.back().stateAtB, false);
    fout << "rowsAtA " << aL.nrActiveRowsAtA;
    fout << " HWAtA " << aL.hammingWeightAtA;
    fout << " dfctWAtB " << aL.knotWeightAtBDeficit;
    fout << " dfctKn " << aL.knotPointDeficit;
    fout << " rowsAtD " << aL.partialNrActiveRowsAtD;
    fout << " HWAtD " << aL.partialHammingWeightAtD;
    fout << " nrKnots " << aL.knots.size();
    fout << " nrChains " << aL.chains.size();
    fout << endl;
    return fout;
}
