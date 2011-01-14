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

#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include "Keccak-fDCLC.h"
#include "Keccak-fPropagation.h"
#include "Keccak-fTrails.h"

const RowValue maskRowValue = 0x1F;

KeccakFPropagation::KeccakFPropagation(const KeccakFDCLC& aParent, KeccakFPropagation::DCorLC aDCorLC)
    : name((aDCorLC == DC) ? "DC" : "LC"), parent(aParent), laneSize(parent.getWidth()/25),
        lambdaMode((aDCorLC == DC) ? KeccakFDCLC::Straight : KeccakFDCLC::Transpose),
        reverseLambdaMode((aDCorLC == DC) ? KeccakFDCLC::Inverse : KeccakFDCLC::Dual),
        directRowOutputListPerInput((aDCorLC == DC) ? aParent.diffChi : aParent.corrInvChi),
        reverseRowOutputListPerInput((aDCorLC == DC) ? aParent.diffInvChi : aParent.corrChi)
{
    initializeAffine();
    initializeWeight();
	initializeMinReverseWeight();
	initializeChiCompatibilityTable();
}

void KeccakFPropagation::initializeAffine()
{
    if (getPropagationType() == KeccakFPropagation::DC) {
        for(RowValue row=0; row<(1<<nrRowsAndColumns); row++) {
            AffineSpaceOfRows a;
            if (row == ((1<<nrRowsAndColumns)-1)) {
                for(unsigned int i=0; i<nrRowsAndColumns-1; i++)
                    a.addGenerator(translateRowSafely(0x3, i));
            }
            else {
                for(int i=0; i<nrRowsAndColumns; i++) {
                    RowValue t = translateRowSafely(row, i);
                    // in:   X--
                    // out: 1000
                    if ((t & 0xE) == 0x2)
                        a.addGenerator(translateRowSafely(0x1, -i));
                    // in:   X-X
                    // out: 1100
                    if ((t & 0xE) == 0xA)
                        a.addGenerator(translateRowSafely(0x3, -i));
                    // in:   XX
                    // out: 100
                    if ((t & 0x6) == 0x6)
                        a.addGenerator(translateRowSafely(0x1, -i));
                    // in:  --X
                    // out: 100
                    if ((t & 0x7) == 0x4)
                        a.addGenerator(translateRowSafely(0x1, -i));
                }
            }
            a.setOffset(parent.chiOnRow(0)^parent.chiOnRow(row));
            affinePerInput.push_back(a);
        }
    }
    else {
        for(RowValue row=0; row<(1<<nrRowsAndColumns); row++) {
            AffineSpaceOfRows a;
            RowValue offs = 0;  
            if (row == ((1<<nrRowsAndColumns)-1)) {
                for(unsigned int i=0; i<nrRowsAndColumns-1; i++)
                    a.addGenerator(translateRowSafely(0x5, i));
                offs = 1;    
            }
            else if (row != 0){
                RowValue t = row;  
                int strt = 0;
                while ( (t&0x1) != 0){
                      strt++ ;
                      t = translateRowSafely(row,-strt);
                }
                int i=0;      
                while (i<nrRowsAndColumns){
                      if ( (t&0x3) == 0x0 ){
                        t = translateRowSafely(t,-1);
                        i +=1;
                      }
                      else if ( (t&0x3) == 0x2 ){
                        offs ^= translateRowSafely(0x2,i+strt);
                        t = translateRowSafely(t,-1);
                        i +=1;
                      }
                      else if ( (t&0x3) == 0x1 ){
                        a.addGenerator(translateRowSafely(0x2,i+strt));
                        a.addGenerator(translateRowSafely(0x4,i+strt));
                        t = translateRowSafely(t,-1);
                        i +=1;
                      }
                      else if ( (t&0x3) == 0x3 ){
                        a.addGenerator(translateRowSafely(0xA,i+strt));
                        a.addGenerator(translateRowSafely(0x4,i+strt));
                        t = translateRowSafely(t,-2);
                        i +=2;
                      }
                }
            }
            a.setOffset(offs);
            affinePerInput.push_back(a);
        }
    }
}

void KeccakFPropagation::initializeWeight()
{
    for(SliceValue slice=0; slice<=maxSliceValue; slice++) {
        unsigned int weightOfThisSlice = weightOfSlice(slice);
        weightPerSlice.push_back(weightOfThisSlice);
    }
}

void KeccakFPropagation::initializeMinReverseWeight()
{
    for(SliceValue slice=0; slice<=maxSliceValue; slice++) {
        unsigned int minReverseWeight = 0;
        for(unsigned int y=0; y<nrRowsAndColumns; y++) {
            RowValue row = getRowFromSlice(slice, y);
            minReverseWeight += reverseRowOutputListPerInput[row].minWeight;
        }
        minReverseWeightPerSlice.push_back(minReverseWeight);
    }
}

KeccakFPropagation::DCorLC KeccakFPropagation::getPropagationType() const
{
    if (lambdaMode == KeccakFDCLC::Straight)
        return DC;
    else if (lambdaMode == KeccakFDCLC::Transpose)
        return LC;
    else
        throw KeccakException("The lambda mode does not match either DC or LC propagation.");
}

unsigned int KeccakFPropagation::weightOfSlice(SliceValue slice) const
{
    unsigned int weight = 0;
    for(unsigned int y=0; y<nrRowsAndColumns; y++) {
        RowValue row = getRowFromSlice(slice, y);
        weight += affinePerInput[row].getWeight();
    }
    return weight;
}

unsigned int KeccakFPropagation::getWeight(const vector<SliceValue>& state) const
{
    unsigned int weight = 0;
    for(unsigned int i=0; i<state.size(); i++)
        weight += getWeight(state[i]);
    return weight;
}

unsigned int KeccakFPropagation::getMinReverseWeight(const vector<SliceValue>& state) const
{
    unsigned int weight = 0;
    for(unsigned int i=0; i<state.size(); i++)
        weight += getMinReverseWeight(state[i]);
    return weight;
}

unsigned int KeccakFPropagation::getMinReverseWeightAfterLambda(const vector<SliceValue>& state) const
{
    vector<SliceValue> stateBeforeLambda;
    reverseLambda(state, stateBeforeLambda);
    return getMinReverseWeight(stateBeforeLambda);
}

void KeccakFPropagation::directPi(unsigned int& dx, unsigned int& dy) const
{
    if ((lambdaMode == KeccakFDCLC::Straight) || (lambdaMode == KeccakFDCLC::Dual))
        KeccakF::pi(dx, dy, dx, dy);
    else
        KeccakF::inversePi(dx, dy, dx, dy);
}

void KeccakFPropagation::reversePi(unsigned int& dx, unsigned int& dy) const
{
    if ((lambdaMode == KeccakFDCLC::Straight) || (lambdaMode == KeccakFDCLC::Dual))
        KeccakF::inversePi(dx, dy, dx, dy);
    else
        KeccakF::pi(dx, dy, dx, dy);
}

void KeccakFPropagation::display(ostream& out) const
{
    if (getPropagationType() == KeccakFPropagation::DC)
        out << "DC analysis tables; patterns are differences." << endl;
    else
        out << "LC analysis tables; patterns are linear masks." << endl;

    vector<vector<RowValue> > rowValuesPerWeight;
    unsigned int maxWeight = 0;
    for(unsigned int i=0; i<affinePerInput.size(); i++)
        if (affinePerInput[i].getWeight() > maxWeight)
            maxWeight = affinePerInput[i].getWeight();

    rowValuesPerWeight.resize(maxWeight+1);
    for(unsigned int i=0; i<affinePerInput.size(); i++)
        rowValuesPerWeight[affinePerInput[i].getWeight()].push_back(i);

    for(unsigned int i=0; i<rowValuesPerWeight.size(); i++) {
        if (rowValuesPerWeight[i].size() > 0) {
            out << "Weight " << dec << i << ": ";
            for(unsigned int j=0; j<rowValuesPerWeight[i].size(); j++) {
                if (j > 0)
                    out << ", ";
                out << hex << (int)rowValuesPerWeight[i][j];
            }
            out << endl;
        }
    }
}

AffineSpaceOfStates KeccakFPropagation::buildStateBase(const vector<SliceValue>& state, bool packedIfPossible) const
{
    static const bool debug = false;
    bool packed = packedIfPossible && ((laneSize*nrRowsAndColumns) <= (sizeof(PackedParity)*8));

    vector<vector<SliceValue> > genValues;   // to store the generator values while processing
    vector<PackedParity> genParitiesPacked;
    vector<vector<RowValue> > genParities;
    vector<SliceValue> offset(laneSize, 0);
    for(unsigned int z=0; z<laneSize; z++) {
        vector<SliceValue> v(laneSize, 0);
        for(unsigned int y=0; y<nrRowsAndColumns; y++) {
            RowValue row = getRowFromSlice(state[z], y);
            offset[z] ^= getSliceFromRow(affinePerInput[row].offset, y);
            for(unsigned int i=0; i<affinePerInput[row].generators.size(); i++) {
                RowValue b = affinePerInput[row].generators[i];
                v[z] = getSliceFromRow(b, y);
                vector<SliceValue> stateAfterLambda;
                parent.lambda(v, stateAfterLambda, lambdaMode);
                genValues.push_back(stateAfterLambda);
                vector<SliceValue> stateBeforeTheta;
                parent.lambdaBeforeTheta(v, stateBeforeTheta, lambdaMode);
                if (packed)
                    genParitiesPacked.push_back(getParity(stateBeforeTheta));
                else {
                    vector<RowValue> parities;
                    getParity(stateBeforeTheta, parities);
                    genParities.push_back(parities);
                }
                if (debug) {
                    cout << "Generator: " << endl;
                    displayStates(cout, v, false, stateBeforeTheta, true, stateAfterLambda, false);
                }
            }
        }
    }
    vector<SliceValue> offsetAfterLambda;
    parent.lambda(offset, offsetAfterLambda, lambdaMode);
    vector<SliceValue> offsetBeforeTheta;
    parent.lambdaBeforeTheta(offset, offsetBeforeTheta, lambdaMode);
    if (debug) {
        cout << "Offset: " << endl;
        displayStates(cout, offset, false, offsetBeforeTheta, true, offsetAfterLambda, false);
    }
    if (packed) {
        PackedParity offsetParitiesPacked = getParity(offsetBeforeTheta);
        AffineSpaceOfStates a(laneSize, genValues, genParitiesPacked, offsetAfterLambda, offsetParitiesPacked);
        return a;
    }
    else {
        vector<RowValue> offsetParities;
        getParity(offsetBeforeTheta, offsetParities);
        AffineSpaceOfStates a(laneSize, genValues, genParities, offsetAfterLambda, offsetParities);
        return a;
    }
}

bool KeccakFPropagation::isThetaJustAfterChi() const
{
    return parent.thetaJustAfterChi[lambdaMode];
}

void KeccakFPropagation::directLambda(const vector<SliceValue>& in, vector<SliceValue>& out) const
{
    parent.lambda(in, out, lambdaMode);
}

void KeccakFPropagation::reverseLambda(const vector<SliceValue>& in, vector<SliceValue>& out) const
{
    parent.lambda(in, out, reverseLambdaMode);
}

void KeccakFPropagation::directLambdaBeforeTheta(const vector<SliceValue>& in, vector<SliceValue>& out) const
{
    parent.lambdaBeforeTheta(in, out, lambdaMode);
}

void KeccakFPropagation::reverseLambdaBeforeTheta(const vector<SliceValue>& in, vector<SliceValue>& out) const
{
    parent.lambdaAfterTheta(in, out, reverseLambdaMode);
}

void KeccakFPropagation::directLambdaAfterTheta(const vector<SliceValue>& in, vector<SliceValue>& out) const
{
    parent.lambdaAfterTheta(in, out, lambdaMode);
}

void KeccakFPropagation::reverseLambdaAfterTheta(const vector<SliceValue>& in, vector<SliceValue>& out) const
{
    parent.lambdaBeforeTheta(in, out, reverseLambdaMode);
}

UINT64 KeccakFPropagation::displayTrailsAndCheck(const string& fileNameIn, ostream& fout, unsigned int maxWeight) const
{
    fout << parent << endl;
    if (getPropagationType() == KeccakFPropagation::DC)
        fout << "Differential cryptanalysis" << endl;
    else
        fout << "Linear cryptanalysis" << endl;
    fout << endl;
    vector<UINT64> countPerWeight, countPerLength;
    UINT64 totalCount = 0;
    unsigned int minWeight = 0;
    {
        ifstream fin(fileNameIn.c_str());
        while(!(fin.eof())) {
            try {
                Trail trail(fin);
                if (getPropagationType() == KeccakFPropagation::DC)
                    parent.checkDCTrail(trail);
                else
                    parent.checkLCTrail(trail);
                if (trail.totalWeight >= countPerWeight.size())
                    countPerWeight.resize(trail.totalWeight+1, 0);
                countPerWeight[trail.totalWeight]++;
                if (trail.states.size() >= countPerLength.size())
                    countPerLength.resize(trail.states.size()+1, 0);
                countPerLength[trail.states.size()]++;
                totalCount++;
            }
            catch(TrailException) {
            }
        }
        if (totalCount == 0) {
            fout << "No trails found in file " << fileNameIn << "!" << endl;
            return totalCount;
        }
        minWeight = 0;
        while((minWeight < countPerWeight.size()) && (countPerWeight[minWeight] == 0))
            minWeight++;
        for(unsigned int i=0; i<countPerLength.size(); i++)
            if (countPerLength[i] > 0)
                fout << dec << countPerLength[i] << " trails of length " << dec << i << " read and checked." << endl;
        fout << "Minimum weight: " << dec << minWeight << endl;
        for(unsigned int i=minWeight; i<countPerWeight.size(); i++) 
            if (countPerWeight[i] > 0) {
                fout.width(8); fout.fill(' ');
                fout << dec << countPerWeight[i] << " trails of weight ";
                fout.width(2); fout.fill(' ');
                fout << i << endl;
            }
        fout << endl;
    }
    if (maxWeight == 0) {
        const unsigned int reasonableNumber = 2000;
        maxWeight = minWeight;
        UINT64 countSoFar = countPerWeight[minWeight];
        while((maxWeight < (countPerWeight.size()-1)) && ((countSoFar+countPerWeight[maxWeight+1]) <= reasonableNumber)) {
            maxWeight++;
            countSoFar += countPerWeight[maxWeight];
        }
    }
    fout << "Showing the trails up to weight " << dec << maxWeight << " (in no particular order)." << endl;
    fout << endl;
    {
        ifstream fin(fileNameIn.c_str());
        while(!(fin.eof())) {
            try {
                Trail trail(fin);
                if (trail.totalWeight <= maxWeight) {
                    trail.display(*this, fout);
                    fout << endl;
                }
            }
            catch(TrailException) {
            }
        }
    }
    return totalCount;
}

void KeccakFPropagation::initializeChiCompatibilityTable()
{
    chiCompatibilityTable.assign(32*32, false);
    for(RowValue a=0; a<32; a++) for(RowValue b=0; b<32; b++) {
        const vector<RowValue>& values = directRowOutputListPerInput[a].values;
        chiCompatibilityTable[a+32*b] = 
            find(values.begin(), values.end(), b) != values.end();
    }
}

bool KeccakFPropagation::isChiCompatible(const vector<SliceValue>& beforeChi, const vector<SliceValue>& afterChi) const
{
    for(unsigned int z=0; z<laneSize; z++)
    for(unsigned int y=0; y<nrRowsAndColumns; y++)
        if (!isChiCompatible(getRowFromSlice(beforeChi[z], y), getRowFromSlice(afterChi[z], y)))
            return false;
    return true;
}

bool KeccakFPropagation::isRoundCompatible(const Trail& first, const Trail& second) const
{
    vector<SliceValue> stateAfterChi;
    reverseLambda(second.states[0], stateAfterChi);
    return isChiCompatible(first.states.back(), stateAfterChi);
}

string KeccakFPropagation::buildFileName(const string& suffix) const
{
    return parent.buildFileName(name, suffix);
}

string KeccakFPropagation::buildFileName(const string& prefix, const string& suffix) const
{
    return parent.buildFileName(name+prefix, suffix);
}

ReverseStateIterator KeccakFPropagation::getReverseStateIterator(const vector<SliceValue>& stateAfterChi, unsigned int maxWeight) const
{
    return ReverseStateIterator(stateAfterChi, *this, maxWeight);
}

ReverseStateIterator::ReverseStateIterator(const vector<SliceValue>& stateAfterChi, const KeccakFPropagation& DCorLC, unsigned int aMaxWeight)
{
    if (aMaxWeight == 0)
        maxWeight = (nrRowsAndColumns-1)*nrRowsAndColumns*stateAfterChi.size();
    else
        maxWeight = aMaxWeight;
    current.assign(stateAfterChi.size(), 0);
    minWeight = 0;
    index = 0;
    size = 0;
    for(unsigned int z=0; z<stateAfterChi.size(); z++)
    for(unsigned int y=0; y<nrRowsAndColumns; y++) {
        RowValue row = getRow(stateAfterChi, y, z);
        if (row != 0) {
            patterns.push_back(DCorLC.reverseRowOutputListPerInput[row]);
            Ys.push_back(y);
            Zs.push_back(z);
            indexes.push_back(0);
            setRow(current, DCorLC.reverseRowOutputListPerInput[row].values[0], y, z);
            minWeight += DCorLC.reverseRowOutputListPerInput[row].weights[0];
            size++;
        }
    }
    currentWeight = minWeight;
    end = isEmpty();
}

bool ReverseStateIterator::isEnd() const
{
    return end;
}

bool ReverseStateIterator::isEmpty() const
{
    return (minWeight > maxWeight) || (size == 0);
}

void ReverseStateIterator::operator++()
{
    next();
    index++;
}

const vector<SliceValue>& ReverseStateIterator::operator*() const
{
    return current;
}

unsigned int ReverseStateIterator::getCurrentWeight() const
{
    return currentWeight;
}

void ReverseStateIterator::next()
{
    int affordableWeight = maxWeight - currentWeight;
    unsigned int i = 0;
    while(i < size) {
        unsigned int ii = indexes[i];
        affordableWeight += patterns[i].weights[ii];
        currentWeight -= patterns[i].weights[ii];
        if (ii < (patterns[i].values.size()-1))
            if ((int)patterns[i].weights[ii+1] <= affordableWeight)
                break;
        affordableWeight -= patterns[i].minWeight;
        i++;
    }
    if (i >= size) {
        end = true;
        return;
    }
    indexes[i]++;
    currentWeight += patterns[i].weights[indexes[i]];
    setRow(current, patterns[i].values[indexes[i]], Ys[i], Zs[i]);
    for(unsigned int j=0; j<i; j++) {
        indexes[j] = 0;
        currentWeight += patterns[j].weights[0];
        setRow(current, patterns[j].values[0], Ys[j], Zs[j]);
    }
}

