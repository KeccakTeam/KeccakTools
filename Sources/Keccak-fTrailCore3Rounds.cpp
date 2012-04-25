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
#include "Keccak-fTrailCore3Rounds.h"
#include "translationsymmetry.h"

using namespace std; 

TrailCore3Rounds::TrailCore3Rounds(const vector<SliceValue>& backgroundAtA,
                                   const vector<SliceValue>& aTabooAtB,
                                   unsigned int aMaxWeight,
                                   const KeccakFDCLC& aParent,
                                   KeccakFPropagation::DCorLC aDCorLC) :
    KeccakFPropagation(aParent, aDCorLC), 
    maxWeight(aMaxWeight), 
    tabooAtB(aTabooAtB)
{
    initializeKnotInfoLUT();

    initializeVortexBase();

    stateAtA.assign(laneSize,0);
    stateAtB.assign(laneSize,0);

    minimumWorkingChainLength = 2;
}

void TrailCore3Rounds::populateStatesWithBackground(const vector<SliceValue>& backgroundAtA)
{
    vector<SliceValue> backgroundAtB(laneSize);
    directLambda(backgroundAtA,backgroundAtB);
    for(unsigned int z=0 ; z<laneSize ; z++) {
        for(unsigned int y=0 ; y<5 ; y++) {
            for(unsigned int x=0 ; x<5 ; x++) {
                if (getBit(backgroundAtB,x,y,z) != 0) {
                    BitPosition pB(x, y, z);
                    addPoint(pB, true, true); // these are all background knot points  
                }
            }
        }
    }
    stateAtA = backgroundAtA;
}

void TrailCore3Rounds::initializeKnotInfoLUT()
{
    if (getPropagationType() == KeccakFPropagation::DC){
        string fileName = "KnotInfoDC" + name + ".cache";
        ifstream fin(fileName.c_str(), ios::binary);
        if (!fin) {
            vector<bool> isTameKnot;
            for(SliceValue s=0; s<=maxSliceValue ; s++) {
                if ( (0  == (s%0x40000)) ) cout << "phase 1 " << s << " of " << maxSliceValue+1 << endl;
                vector<SliceValue> state(laneSize,0);
                AffineSpaceOfSlices base = buildSliceBase(s);
                SliceValue dummy;
                bool isTame = base.getOffsetWithGivenParity(0,dummy);
                isTameKnot.push_back(isTame);
                if ((getHammingWeightSlice(s) == 2) && (getParity(s) == 0)) isTameKnot.back() = false;
            }
            for(SliceValue s=0; s<=maxSliceValue ; s++) {
                if ( (0  == (s%0x40000)) ) cout << "phase 2 " << s << " of " << maxSliceValue+1 << endl;
                unsigned int knotPointDeficit;
                unsigned int knotWeightAtBDeficit;
                unsigned int nrActiveRows = getNrActiveRows(s);
                bool isOrbital;
                if (isTameKnot[s]) {
                    knotPointDeficit = 0;
                    knotWeightAtBDeficit = 0;
                    isOrbital = false;
                }
                else {
                    unsigned int HW = getHammingWeightSlice(s);
                    if (HW == 0){
                        knotPointDeficit = 0;
                        knotWeightAtBDeficit = 0;
                        isOrbital = false;
                    }
                    else if (HW == 1){
                        knotPointDeficit = 2;
                        knotWeightAtBDeficit = 3;
                        isOrbital = false;
                    }
                    else if ((HW == 2) && (getParity(s) == 0)){
                        knotPointDeficit = 1;
                        knotWeightAtBDeficit = 1;
                        isOrbital = true;
                    }
                    else {
                        knotPointDeficit = 1;
                        isOrbital = false;
                        knotWeightAtBDeficit = 4;
                        for (unsigned int x=0 ; x<5 ; x++){
                            for (unsigned int y=0 ; y<5 ; y++){
                                SliceValue strayBit = getSlicePoint(x,y);
                                strayBit |= s;
                                if (isTameKnot[strayBit]){
                                    int tmp = getWeight(strayBit) - getWeight(s);
                                    if (tmp < (int) knotWeightAtBDeficit) knotWeightAtBDeficit =  (unsigned int) tmp;
                                }
                            }
                        }
                    }
                }
                knotInfoLUT.push_back(packKnotInfo(knotPointDeficit,knotWeightAtBDeficit,nrActiveRows,isOrbital));
            }
            ofstream fout(fileName.c_str(), ios::binary);
            for(SliceValue s=0 ; s<=maxSliceValue ; s++ ) {
                static unsigned char tmp[1];
                tmp[0] =  knotInfoLUT[s];
                fout.write((char *)tmp, 1);
            }
        }
        else {
            for(SliceValue s=0 ; s<=maxSliceValue ; s++ ) {
                unsigned char tmp[1];
                fin.read((char *)tmp,1);
                knotInfoLUT.push_back(tmp[0]);
            }
        }
    }
    else throw KeccakException("implementation of KnotInfoLUT for LC is under construction");
}

UINT8 TrailCore3Rounds::packKnotInfo(unsigned int knotPointDeficit,unsigned int knotWeightAtBDeficit, unsigned int nrActiveRows, bool isOrbital) const
{
    UINT8 tmp = nrActiveRows*2 + knotPointDeficit*16 + knotWeightAtBDeficit*64;
    if (isOrbital) return (1+tmp);
    else return tmp;
}

void TrailCore3Rounds::populateKnotInfo(KnotInformation& aKnotInfo,const SliceValue& aSliceValue, bool knotHasSinglePoint, bool hasBackground) const
{
    unsigned int tmp            = knotInfoLUT[aSliceValue];
    aKnotInfo.isOrbital         = (1 == (tmp&1));
    if (hasBackground && knotHasSinglePoint){
        aKnotInfo.nrActiveRows         = 1;
        aKnotInfo.knotPointDeficit     = 1;
        aKnotInfo.knotWeightAtBDeficit = 2;
    }
    else if (hasBackground && aKnotInfo.isOrbital) {
        aKnotInfo.nrActiveRows         = 2;
        aKnotInfo.knotPointDeficit     = 0;
        aKnotInfo.knotWeightAtBDeficit = 0;
    }
    else {
        tmp >>= 1; aKnotInfo.nrActiveRows         = tmp&0x7;
        tmp >>= 3; aKnotInfo.knotPointDeficit     = tmp&0x3;
        tmp >>= 2; aKnotInfo.knotWeightAtBDeficit = tmp&0x3;
    }
}

bool TrailCore3Rounds::canAffordExtendingChain() const
{
    unsigned int nrOrbitalPointsPerDeltaChain;
    unsigned int deltaNrOrbitalPointsWorkingChain;

    if (chains.back().size() >= minimumWorkingChainLength) {
        int chainOdd = chains.back().size()%2; // the length of a chain, and also its number of orbital points, is even
        deltaNrOrbitalPointsWorkingChain = 1 - chainOdd; // If the working chain has an odd number of points, the next point can be a knot point. Only if the working chain has an even number of points, the next point is for sure an orbital point.
        nrOrbitalPointsPerDeltaChain = chains.back().size() - chainOdd; // If the working chain has an odd number of points, it contains (its size minus 1) orbital points. If the working chain has an even number of points, it currently contains (its size minus 2) orbital points, which becomes (its size) after extension.
    }
    else {
        deltaNrOrbitalPointsWorkingChain = minimumWorkingChainLength - chains.back().size() - 1; // After extension, the working chain will contain at least (its current size minus 1) orbital points. The working chain has to reach at least (minimumWorkingChainLength - 2) orbital points.
        nrOrbitalPointsPerDeltaChain = minimumWorkingChainLength - 2; // Any new chain has to reach at least (minimumWorkingChainLength - 2) orbital points.
    }
    return canAffordGeneric(1, 0, nrOrbitalPointsPerDeltaChain, (unsigned int) deltaNrOrbitalPointsWorkingChain); // By extending the chain, at least one point will be added, but not necessarily a knot or a run.
}

bool TrailCore3Rounds::canAffordAddingChain() const
{
    // Adding a chain necessarily adds two knot points. The new chain, as well as future ones, must have at least (minimumWorkingChainLength - 2) orbital points.
    return canAffordGeneric(2, 0, minimumWorkingChainLength - 2, minimumWorkingChainLength - 2);
}

// completeChain assumes the working chain has at least one point: the start point.
// It will not modify this point and if unsuccessful, the working chain contains this point.
bool TrailCore3Rounds::completeChain()
{
    bool canChainBeExtended = canAffordExtendingChain();
    do {
        if (!canChainBeExtended) {
            if (chains.back().size() == 1) return false;
            else removePoint(false);
        }
        while((chains.back().size() > 1) && (yOffsets.back().back() == 4)) 
            removePoint(false); // yOffset equal to 4 does not have a successor
        if (yOffsets.back().back() == 4) 
            return false;

        if (((chains.back().size())%2) == 0) {   // In this case, we have to add the second point in an orbital slice.
            bool success;
            BitPosition pB;
            do {
                yOffsets.back().back() += 1; // this is the y offset at B
                pB = chains.back().back();
                pB.yTranslate(yOffsets.back().back());
                success = (0 == getBit(tabooAtB, pB)); // pB shall not lie in tabooAtB
            }
            while((!success) && (yOffsets.back().back() < 4));
            if (success) {
                addPoint(pB, false);
                canChainBeExtended = canAffordExtendingChain();
            }
        }
        else {  // In this case, we have to add the first point in an empty slice (to become an orbital slice) or the end point of the chain (in a knot).
            bool success;
            BitPosition pB;
            bool isPotentialEndPoint;
            do {
                yOffsets.back().back() += 1;        // this is the offset at A
                pB = chains.back().back();
                reverseRhoPi(pB);
                pB.yTranslate(yOffsets.back().back());
                directRhoPi(pB);
                isPotentialEndPoint = 
                    (chains.back().size()+1 >= minimumWorkingChainLength )  // chain is long enough
                    && (chains.back()[0] < pB);                             // end point is larger than start point
                success = (0 == getBit(tabooAtB, pB));                      // it does not arrive in the taboo zone
                if (success) {
                    if (knots.count(pB.z) != 0) {                           // if it arrives in a knot slice
                        success = 
                            (getBit(stateAtB, pB) == 0) //  it does not collide with an already active point
                            && isPotentialEndPoint;     //  it satisfies all conditions of an end point
                    }
                    else {
                        success = (0 == stateAtB[pB.z]);                    // it arrives in an empty slice
                        isPotentialEndPoint = isPotentialEndPoint && success && mayBeEndPoint(pB);  // can we afford adding a knot (in kernel at C) or a run (out kernel at C)?
                    }
                }
            }
            while((!success) && (yOffsets.back().back() < 4));

            if (success) {
                addPoint(pB,isPotentialEndPoint);
                if (isPotentialEndPoint)
                    return true;
                else
                    canChainBeExtended = canAffordExtendingChain();
            }
        }
    }
    while(true);
}

void TrailCore3Rounds::updateMinimumWorkingChainLength()
{
    if (chains.size() <= 1)
        minimumWorkingChainLength = 2; // if there is only one chain or none, the minimum working chain length is 2
    else {
        // If there are more than one chain, then they must be in non-decreasing order of size, so the working chain must reach at least the size of the last one.
        minimumWorkingChainLength = chains[chains.size()-2].size();
        // Within a group of chains of equal sizes, their start point must be in increasing order. Hence, if the working chain starts at a position less than the start point of the last chain, it necessarily has to belong to a group of larger size.
        if (chains.back()[0] < chains[chains.size()-2][0])
            minimumWorkingChainLength += 2;
    }
}

bool TrailCore3Rounds::nextChain()
{
    if (chains.back().size() > 1) {
        // If the last added point introduced a new knot, it has to be first converted to an orbital point, as the completeChain() function assumes the working chain ends with an orbital point.
        if (knotPointAddedKnot.top()) 
            convertKnotPointToOrbitalPoint(); // end point becomes intermediate point
        else
            removePoint(true);
    }
    bool success = true;
    do {
        // First, let's see if we can find a chain starting with the current start point.
        if ((chains.back().size() > 0) && (completeChain()))
            return true;
        else {  // If not, we must iterate over the start point.
            BitPosition pB;
            bool pBisInitialized = (chains.back().size() > 0);
            if (pBisInitialized) {
                pB = chains.back()[0];
                removePoint(true);
            }
            do {
                if ((!pBisInitialized) || (!pB.nextXY())) {
                    if (pBisInitialized && (knots.size() == 0)) return false; // implying this is the start point of first chain and the background is zero in which case the first chain shall start in slice 0
                    pB.x = 0;
                    pB.y = 0;
                    if (startPointWorkingChainIsFree) {
                        // If the start point can be freely chosen, we simply run through the different z coordinates.
                        if (!pBisInitialized) {
                            pB.z = 0;
                            pBisInitialized = true;
                        }
                        else
                            pB.z = pB.z + 1;
                        if (pB.z == laneSize) return false;
                    }
                    else {
                        // Otherwise, we rely on mayBeStartPointSliceAndGoThere().
                        if (!mayBeStartPointSliceAndGoThere(pB.z, pBisInitialized)) return false;
                        pBisInitialized = true;
                    }
                }
            }
            while((getBit(stateAtB, pB) != 0)  ||   // shall not collide with existing point
                  (getBit(tabooAtB, pB) != 0)  ||   // shall not be in the taboo zone
                  ( (stateAtB[pB.z] != 0) && (knots.count(pB.z) == 0) ) );  // shall not be in an orbital slice
            addPoint(pB, true);
            updateMinimumWorkingChainLength(); 
        }
    }
    while(true);
}

bool TrailCore3Rounds::nextWithKnots()
{
    if ((!knots.empty()) && chains.empty() && isStateAtBWellFormed()) return true; // this deals with non-zero backgrounds that result in a well-formed state at B
    do {
        if (knots.empty() || (canAffordAddingChain())) {
            chains.push_back(vector<BitPosition>());
            yOffsets.push_back(vector<unsigned int>());
            // Having a "free" starting point implies a knot or a run can be added. When this function is called, the working chain is empty, hence it necessarily adds two knot points. The new chain, as well as future ones, must have at least (minimumWorkingChainLength - 2) orbital points.
            startPointWorkingChainIsFree = canAffordGeneric(2, 1, minimumWorkingChainLength - 2, minimumWorkingChainLength - 2);
        }
        while((!chains.empty()) && (!nextChain())) {
            chains.pop_back();
            yOffsets.pop_back();
            updateMinimumWorkingChainLength();
        }
        if (chains.empty()) return false;
    }
    while (!isStateAtBWellFormed()); 
    return true;
}

void TrailCore3Rounds::addVortexPoint(const BitPosition& pB,
                                      vector<BitPosition >& chainAtB,
                                      vector<unsigned int>& yOffset,
                                      map<RowPosition,unsigned int>& rowsAtA, 
                                      map<RowPosition,unsigned int>& rowsAtD,
                                      map<unsigned int,unsigned int>& slicesAtB) const
{
    chainAtB.push_back(pB);
    yOffset.push_back(0);
    
    slicesAtB[pB.z] += 1;

    BitPosition pA(pB); 
    reverseRhoPi(pA);
    rowsAtA[RowPosition(pA)] += 1;
    
    BitPosition pD(pB); 
    directRhoPi(pD); // For an orbital point at B, the corresponding bit position at C is active. If C is in the kernel, then this active bit will appear at D after applying ρ and π.
    rowsAtD[RowPosition(pD)] += 1;
}

void TrailCore3Rounds::removeVortexPoint(vector<BitPosition >& chainAtB,
                                        vector<unsigned int>& yOffset,
                                        map<RowPosition,unsigned int>& rowsAtA, 
                                        map<RowPosition,unsigned int>& rowsAtD,
                                        map<unsigned int,unsigned int>& slicesAtB) const
{
    BitPosition pB(chainAtB.back());
    if (slicesAtB[pB.z] == 1) slicesAtB.erase(pB.z);
    else slicesAtB[pB.z] -= 1;

    BitPosition pA(pB); 
    reverseRhoPi(pA); 
    RowPosition rA(pA); 
    if (rowsAtA[rA] == 1) rowsAtA.erase(rA); 
    else rowsAtA[rA] -= 1;

    BitPosition pD(pB); 
    directRhoPi(pD);        
    RowPosition rD(pD); 
    if (rowsAtD[rD] == 1) rowsAtD.erase(rD); 
    else rowsAtD[rD] -= 1;

    chainAtB.pop_back();
    yOffset.pop_back();
}

void TrailCore3Rounds::addVortexToBaseIfMinimal(const vector<BitPosition >& chainAtB, unsigned int theNrActiveRowsAtA, unsigned int theNrActiveRowsAtD)
{
    vector<SliceValue> workState(laneSize,0);
    for (unsigned int i=0 ; i<chainAtB.size() ; i++) setBitToOne(workState,chainAtB[i]);

    if ((chainAtB[0] < chainAtB.back()) && (isMinimalSymmetrically(workState))) {
        while(2*vortexBase.size() <= chainAtB.size()) vortexBase.push_back(vector<VortexInfo>(0));
        VortexInfo workVortexInfo;
        workVortexInfo.nrActiveRowsAtA = theNrActiveRowsAtA;
        workVortexInfo.nrActiveRowsAtD = theNrActiveRowsAtD;
        for(unsigned int i=0 ; i<chainAtB.size() ; i++) workVortexInfo.stateAtB.setBitToOne(chainAtB[i]);
        vortexBase[chainAtB.size()/2].push_back(workVortexInfo);
    }
}

void TrailCore3Rounds::initializeVortexBase()
{
    vector<BitPosition> chainAtB;
    vector<unsigned int> yOffset;
    map<RowPosition,unsigned int> rowsAtA;
    map<RowPosition,unsigned int> rowsAtD;
    map<unsigned int,unsigned int> slicesAtB;

    BitPosition pB(0,0,0);
    addVortexPoint(pB,chainAtB,yOffset,rowsAtA,rowsAtD,slicesAtB);

    do {
        if ((chainAtB.size() == 1) && (yOffset.back() == 4)) { // time to iterate the first point of the working chain
            pB = chainAtB.back();
            removeVortexPoint(chainAtB,yOffset,rowsAtA,rowsAtD,slicesAtB);
            if (!pB.nextXY()) return;
            addVortexPoint(pB,chainAtB,yOffset,rowsAtA,rowsAtD,slicesAtB);
        }

        if ((chainAtB.size()%2) == 0) { // The next point to add has to be chained to the last one (i.e., in the same column at B).
            yOffset.back() += 1;
            pB = chainAtB.back();
            pB.yTranslate(yOffset.back());
            addVortexPoint(pB,chainAtB,yOffset,rowsAtA,rowsAtD,slicesAtB); 
            if (2*chainAtB.size() + 2*rowsAtA.size() + 2*rowsAtD.size() + 2 > maxWeight) removeVortexPoint(chainAtB,yOffset,rowsAtA,rowsAtD,slicesAtB); // if the maximum weight does not allow adding one more point afterwards, remove this point again
        }
        else { // The next point to add has to be peer to the last one (i.e., in the same column at A).
            yOffset.back() += 1;
            pB = chainAtB.back();
            reverseRhoPi(pB);
            pB.yTranslate(yOffset.back());
            directRhoPi(pB);
            if ( slicesAtB.count(pB.z) == 0 ) {  // if the point arrives in an empty slice, add it
                addVortexPoint(pB,chainAtB,yOffset,rowsAtA,rowsAtD,slicesAtB);
                if (2*chainAtB.size() + 2*rowsAtA.size() + 2*rowsAtD.size() + 4 > maxWeight) removeVortexPoint(chainAtB,yOffset,rowsAtA,rowsAtD,slicesAtB); // if the maximum weight does not allow adding two more points afterwards, remove this point again
            }
            else if ((pB.z == chainAtB[0].z) && (pB.x == chainAtB[0].x) && (pB.y != chainAtB[0].y)) { // we have hit a vortex
                addVortexPoint(pB,chainAtB,yOffset,rowsAtA,rowsAtD,slicesAtB);
                if (2*chainAtB.size() + 2*rowsAtA.size() + 2*rowsAtD.size() <= maxWeight) {
                    addVortexToBaseIfMinimal(chainAtB,rowsAtA.size(),rowsAtD.size());
                }
                removeVortexPoint(chainAtB,yOffset,rowsAtA,rowsAtD,slicesAtB); // after storing the vortex, remove the last point of the working chain
            }
        }
        while((chainAtB.size() > 1) && (yOffset.back() == 4))
            removeVortexPoint(chainAtB,yOffset,rowsAtA,rowsAtD,slicesAtB); // purging the trailing part of the working chain
    }
    while(true);
}

