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
#include <sstream>
#include "Keccak-fDCLC.h"
#include "Keccak-fPropagation.h"

using namespace std;

// -------------------------------------------------------------
//
// ListOfRowPatterns
//
// -------------------------------------------------------------

void ListOfRowPatterns::add(RowValue aValue, int aWeight)
{
    unsigned int i = 0;
    while((i < values.size()) && (weights[i] <= aWeight))
        i++;
    if (i < values.size()) {
        values.insert(values.begin()+i, aValue);
        weights.insert(weights.begin()+i, aWeight);
    }
    else {
        values.push_back(aValue);
        weights.push_back(aWeight);
    }
    if ((!minMaxInitialized) || (aWeight > maxWeight))
        maxWeight = aWeight;
    if ((!minMaxInitialized) || (aWeight < minWeight))
        minWeight = aWeight;
    minMaxInitialized = true;
}

void ListOfRowPatterns::display(ostream& fout) const
{
    for(unsigned int i=0; i<values.size(); i++) {
        fout << dec << weights[i] << "  ";
        fout.fill('0'); fout.width(2); fout << hex << (int)values[i] << endl;
    }
}


// -------------------------------------------------------------
//
// KeccakFDCLC
//
// -------------------------------------------------------------

KeccakFDCLC::KeccakFDCLC(unsigned int aWidth)
    : KeccakF(aWidth)
{
    initializeAll();
}

string KeccakFDCLC::getDescription() const
{
    return "DC/LC analysis of " + KeccakF::getDescription();
}

RowValue KeccakFDCLC::chiOnRow(RowValue a) const
{
    vector<LaneValue> A(nrRowsAndColumns*nrRowsAndColumns, 0);
    setRow(A, a);
    chi(A);
    return getRow(A);
}

RowValue KeccakFDCLC::inverseChiOnRow(RowValue a) const
{
    vector<LaneValue> A(nrRowsAndColumns*nrRowsAndColumns, 0);
    setRow(A, a);
    inverseChi(A);
    return getRow(A);
}

template<class T>
unsigned int dotProduct(const T& a, const T& b)
{
    T aandb = a&b;
    unsigned int result = 0;
    while(aandb != 0) {
        result ^= (aandb & 1);
        aandb >>= 1;
    }
    return result;
}

template<class T>
unsigned int dotProduct(const vector<T>& a, const vector<T>& b)
{
    unsigned int result = 0;
    for(unsigned int i=0; (i<a.size()) && (i<b.size()); ++i)
        result ^= dotProduct(a[i], b[i]);
    return result;
}

int computeDifferentialWeight(int count)
{
    return (int)floor(nrRowsAndColumns - log(fabs((double)count))/log(2.0)+0.5);
}

int computeLinearWeight(int correl)
{
    return 2*(int)floor(nrRowsAndColumns - log(fabs((double)correl))/log(2.0)+0.5);
}

void KeccakFDCLC::initializeAll()
{
    for(RowValue da=0; da<(1<<nrRowsAndColumns); da++) {
        vector<unsigned int> count(1<<nrRowsAndColumns, 0);
        for(RowValue a=0; a<(1<<nrRowsAndColumns); a++) {
            RowValue db = chiOnRow(a)^chiOnRow(a^da);
            count[db]++;
        }
        ListOfRowPatterns l;
        for(RowValue db=0; db<(1<<nrRowsAndColumns); db++)
            if (count[db] != 0)
                l.add(db, computeDifferentialWeight(count[db]));
        diffChi.push_back(l);
    }
    for(RowValue da=0; da<(1<<nrRowsAndColumns); da++) {
        vector<unsigned int> count(1<<nrRowsAndColumns, 0);
        for(RowValue a=0; a<(1<<nrRowsAndColumns); a++) {
            RowValue db = inverseChiOnRow(a)^inverseChiOnRow(a^da);
            count[db]++;
        }
        ListOfRowPatterns l;
        for(RowValue db=0; db<(1<<nrRowsAndColumns); db++)
            if (count[db] != 0)
                l.add(db, computeDifferentialWeight(count[db]));
        diffInvChi.push_back(l);
    }
    for(RowValue ua=0; ua<(1<<nrRowsAndColumns); ua++) {
        vector<int> correl(1<<nrRowsAndColumns, 0);
        for(RowValue ub=0; ub<(1<<nrRowsAndColumns); ub++) {
            vector<unsigned int> count(2, 0);
            for(RowValue a=0; a<(1<<nrRowsAndColumns); a++) {
                RowValue b = chiOnRow(a);
                count[dotProduct(a, ua)^dotProduct(b, ub)]++;
            }
            correl[ub] = count[0]-count[1];
        }
        ListOfRowPatterns l;
        for(RowValue ub=0; ub<(1<<nrRowsAndColumns); ub++)
            if (correl[ub] != 0)
                l.add(ub, computeLinearWeight(correl[ub]));
        corrChi.push_back(l);
    }
    for(RowValue ua=0; ua<(1<<nrRowsAndColumns); ua++) {
        vector<int> correl(1<<nrRowsAndColumns, 0);
        for(RowValue ub=0; ub<(1<<nrRowsAndColumns); ub++) {
            vector<unsigned int> count(2, 0);
            for(RowValue a=0; a<(1<<nrRowsAndColumns); a++) {
                RowValue b = inverseChiOnRow(a);
                count[dotProduct(a, ua)^dotProduct(b, ub)]++;
            }
            correl[ub] = count[0]-count[1];
        }
        ListOfRowPatterns l;
        for(RowValue ub=0; ub<(1<<nrRowsAndColumns); ub++)
            if (correl[ub] != 0)
                l.add(ub, computeLinearWeight(correl[ub]));
        corrInvChi.push_back(l);
    }
    initializeLambdaLookupTables();
}

void KeccakFDCLC::displayAll(ostream& fout, KeccakFPropagation *DC, KeccakFPropagation *LC) const
{
    for(RowValue i=0; i<(1<<nrRowsAndColumns); i++) {
        fout << "Difference of \xCF\x87 in one row: ";
        fout.fill('0'); fout.width(2); fout << hex << (int)i << endl;
        diffChi[i].display(fout);
        if (DC) {
            fout << "Affine description: ";
            DC->affinePerInput[i].display(fout);
            fout << endl;
        }
    }
    for(RowValue i=0; i<(1<<nrRowsAndColumns); i++) {
        fout << "Difference of \xCF\x87^-1 in one row: ";
        fout.fill('0'); fout.width(2); fout << hex << (int)i << endl;
        diffInvChi[i].display(fout);
        fout << endl;
    }
    for(RowValue i=0; i<(1<<nrRowsAndColumns); i++) {
        fout << "Correlation of \xCF\x87 in one row: ";
        fout.fill('0'); fout.width(2); fout << hex << (int)i << endl;
        corrChi[i].display(fout);
        fout << endl;
    }
    for(RowValue i=0; i<(1<<nrRowsAndColumns); i++) {
        fout << "Correlation of \xCF\x87^-1 in one row: ";
        fout.fill('0'); fout.width(2); fout << hex << (int)i << endl;
        corrInvChi[i].display(fout);
        if (LC) {
            fout << "Affine description: ";
            LC->affinePerInput[i].display(fout);
            fout << endl;
        }
    }
}

void KeccakFDCLC::initializeLambdaLookupTables()
{
    // lambdaRowToSlice
    {
        string fileName = buildFileName("", "-lambda.cache");
        ifstream fin(fileName.c_str(), ios::binary);
        if (!fin) {
            for(unsigned int m=0; m<KeccakFDCLC::EndOfLambdaModes ; m++) {
                vector<vector<vector<vector<SliceValue> > > > vvvv;
                for(unsigned int outputSlice=0; outputSlice<laneSize; outputSlice++) {
                    vector<vector<vector<SliceValue> > > vvv;
                    for(unsigned int inputSlice=0; inputSlice<laneSize; inputSlice++) {
                        vector<vector<SliceValue> > vv;
                        for(unsigned int y=0; y<nrRowsAndColumns; y++) {
                            vector<SliceValue> v(1<<nrRowsAndColumns);
                            for(RowValue row=0; row<(1<<nrRowsAndColumns); row++) {
                                vector<LaneValue> state(nrRowsAndColumns*nrRowsAndColumns, 0);
                                setRow(state, row, y, inputSlice);
                                lambda(state,KeccakFDCLC::LambdaMode(m));
                                v[row] = getSlice(state, outputSlice);
                            }
                            vv.push_back(v);
                        }
                        vvv.push_back(vv);
                    }
                    vvvv.push_back(vvv);
                }
                lambdaRowToSlice.push_back(vvvv);
            }
            ofstream fout(fileName.c_str(), ios::binary);
            for(unsigned int m=0; m<KeccakFDCLC::EndOfLambdaModes ; m++)
            for(unsigned int outputSlice=0; outputSlice<laneSize; outputSlice++)
            for(unsigned int inputSlice=0; inputSlice<laneSize; inputSlice++)
            for(unsigned int y=0; y<nrRowsAndColumns; y++)
            for(RowValue row=0; row<(1<<nrRowsAndColumns); row++) {
                static unsigned char tmp[4];
                SliceValue v = lambdaRowToSlice[m][outputSlice][inputSlice][y][row];
                tmp[0] =  v&0xFF;
                tmp[1] = (v>>8)&0xFF;
                tmp[2] = (v>>16)&0xFF;
                tmp[3] = (v>>24)&0xFF;
                fout.write((char *)tmp, 4);
            }
        }
        else {
            for(unsigned int m=0; m<KeccakFDCLC::EndOfLambdaModes ; m++) {
                vector<vector<vector<vector<SliceValue> > > > vvvv;
                for(unsigned int outputSlice=0; outputSlice<laneSize; outputSlice++) {
                    vector<vector<vector<SliceValue> > > vvv;
                    for(unsigned int inputSlice=0; inputSlice<laneSize; inputSlice++) {
                        vector<vector<SliceValue> > vv;
                        for(unsigned int y=0; y<nrRowsAndColumns; y++) {
                            vector<SliceValue> v(1<<nrRowsAndColumns);
                            for(RowValue row=0; row<(1<<nrRowsAndColumns); row++) {
                                static unsigned char tmp[4];
                                fin.read((char *)tmp, 4);
                                v[row]  = tmp[3];  v[row] <<= 8;
                                v[row] ^= tmp[2];  v[row] <<= 8;
                                v[row] ^= tmp[1];  v[row] <<= 8;
                                v[row] ^= tmp[0];
                            }
                            vv.push_back(v);
                        }
                        vvv.push_back(vv);
                    }
                    vvvv.push_back(vvv);
                }
                lambdaRowToSlice.push_back(vvvv);
            }
        }
    }
    // lambdaBeforeThetaRowToSlice
    for(unsigned int m=0; m<KeccakFDCLC::EndOfLambdaModes ; m++) {
        vector<vector<vector<vector<SliceValue> > > > vvvv;
        for(unsigned int outputSlice=0; outputSlice<laneSize; outputSlice++) {
            vector<vector<vector<SliceValue> > > vvv;
            for(unsigned int inputSlice=0; inputSlice<laneSize; inputSlice++) {
                vector<vector<SliceValue> > vv;
                for(unsigned int y=0; y<nrRowsAndColumns; y++) {
                    vector<SliceValue> v(1<<nrRowsAndColumns);
                    for(RowValue row=0; row<(1<<nrRowsAndColumns); row++) {
                        vector<LaneValue> state(nrRowsAndColumns*nrRowsAndColumns, 0);
                        setRow(state, row, y, inputSlice);
                        lambdaBeforeTheta(state, KeccakFDCLC::LambdaMode(m));
                        v[row] = getSlice(state, outputSlice);
                    }
                    vv.push_back(v);
                }
                vvv.push_back(vv);
            }
            vvvv.push_back(vvv);
        }
        lambdaBeforeThetaRowToSlice.push_back(vvvv);
    }
    // thetaJustAfterChi
    for(unsigned int mode=0; mode<KeccakFDCLC::EndOfLambdaModes ; mode++) {
        if ((mode == Straight) || (mode == Dual))
            thetaJustAfterChi.push_back(true);
        else if ((mode == Transpose) || (mode == Inverse))
            thetaJustAfterChi.push_back(false);
    }
    // thetaJustBeforeChi
    for(unsigned int mode=0; mode<KeccakFDCLC::EndOfLambdaModes ; mode++) {
        if ((mode == Straight) || (mode == Dual))
            thetaJustBeforeChi.push_back(false);
        else if ((mode == Transpose) || (mode == Inverse))
            thetaJustBeforeChi.push_back(true);
    }
    // lambdaAfterThetaRowToSlice
    for(unsigned int m=0; m<KeccakFDCLC::EndOfLambdaModes ; m++) {
        vector<vector<vector<vector<SliceValue> > > > vvvv;
        for(unsigned int outputSlice=0; outputSlice<laneSize; outputSlice++) {
            vector<vector<vector<SliceValue> > > vvv;
            for(unsigned int inputSlice=0; inputSlice<laneSize; inputSlice++) {
                vector<vector<SliceValue> > vv;
                for(unsigned int y=0; y<nrRowsAndColumns; y++) {
                    vector<SliceValue> v(1<<nrRowsAndColumns);
                    for(RowValue row=0; row<(1<<nrRowsAndColumns); row++) {
                        vector<LaneValue> state(nrRowsAndColumns*nrRowsAndColumns, 0);
                        setRow(state, row, y, inputSlice);
                        lambdaAfterTheta(state, KeccakFDCLC::LambdaMode(m));
                        v[row] = getSlice(state, outputSlice);
                    }
                    vv.push_back(v);
                }
                vvv.push_back(vv);
            }
            vvvv.push_back(vvv);
        }
        lambdaAfterThetaRowToSlice.push_back(vvvv);
    }
}

void KeccakFDCLC::lambda(const vector<SliceValue>& in, vector<SliceValue>& out, LambdaMode mode) const
{
    // This assumes that 'in' has size equal to 'laneSize'
    out.assign(laneSize, 0);
    for(unsigned int inputSlice=0; inputSlice<laneSize; inputSlice++)
    for(unsigned int y=0; y<nrRowsAndColumns; y++) {
        RowValue row = getRowFromSlice(in[inputSlice], y);
        for(unsigned int outputSlice=0; outputSlice<laneSize; outputSlice++)
            out[outputSlice] ^= lambdaRowToSlice[mode][outputSlice][inputSlice][y][row];
    }
}

void KeccakFDCLC::lambdaBeforeTheta(const vector<SliceValue>& in, vector<SliceValue>& out, LambdaMode mode) const
{
    if (thetaJustAfterChi[mode]) {
        out = in;
    }
    else {
        // This assumes that 'in' has size equal to 'laneSize'
        out.assign(laneSize, 0);
        for(unsigned int inputSlice=0; inputSlice<laneSize; inputSlice++)
        for(unsigned int y=0; y<nrRowsAndColumns; y++) {
            RowValue row = getRowFromSlice(in[inputSlice], y);
            for(unsigned int outputSlice=0; outputSlice<laneSize; outputSlice++)
                out[outputSlice] ^= lambdaBeforeThetaRowToSlice[mode][outputSlice][inputSlice][y][row];
        }
    }
}

void KeccakFDCLC::lambdaAfterTheta(const vector<SliceValue>& in, vector<SliceValue>& out, LambdaMode mode) const
{
    if (thetaJustBeforeChi[mode]) {
        out = in;
    }
    else {
        // This assumes that 'in' has size equal to 'laneSize'
        out.assign(laneSize, 0);
        for(unsigned int inputSlice=0; inputSlice<laneSize; inputSlice++)
        for(unsigned int y=0; y<nrRowsAndColumns; y++) {
            RowValue row = getRowFromSlice(in[inputSlice], y);
            for(unsigned int outputSlice=0; outputSlice<laneSize; outputSlice++)
                out[outputSlice] ^= lambdaAfterThetaRowToSlice[mode][outputSlice][inputSlice][y][row];
        }
    }
}

void KeccakFDCLC::checkDCTrail(const Trail& trail, KeccakFPropagation *DC) const
{
    // Check weights
    unsigned int totalWeight = 0;
    unsigned int offsetIndex = (trail.firstStateSpecified ? 0 : 1);
    if ((!trail.firstStateSpecified) && (trail.weights.size() >= 1))
        totalWeight += trail.weights[0];
    for(unsigned int i=offsetIndex; i<trail.weights.size(); i++) {
        unsigned int weight = 0;
        for(unsigned int z=0; z<laneSize; z++)
            for(unsigned int y=0; y<nrRowsAndColumns; y++)
                weight += diffChi[getRowFromSlice(trail.states[i][z], y)].minWeight;
        if (weight != trail.weights[i]) {
            if (DC) trail.display(*DC, cerr);
            cerr << "The weight of state at round " << dec << i << " is incorrect; it should be " << weight << "." << endl;
            throw KeccakException("The weights in the trail are incorrect!");
        }
        totalWeight += weight;
    }
    if (totalWeight != trail.totalWeight) {
        if (DC) trail.display(*DC, cerr);
        cerr << "The total weight of the trail is incorrect; it should be " << totalWeight << "." << endl;
        throw KeccakException("The total weight in the trail is incorrect!");
    }

    // Check compatibility between consecutive states
    for(unsigned int i=1+offsetIndex; i<trail.states.size(); i++) {
        vector<SliceValue> stateAfterChi;
        lambda(trail.states[i], stateAfterChi, Inverse);
        for(unsigned int z=0; z<laneSize; z++)
            for(unsigned int y=0; y<nrRowsAndColumns; y++) {
                const vector<RowValue>& values = diffChi[getRowFromSlice(trail.states[i-1][z], y)].values;
                if (find(values.begin(), values.end(), getRowFromSlice(stateAfterChi[z], y)) == values.end()) {
                    if (DC) trail.display(*DC, cerr);
                    cerr << "The state at round " << dec << i-1 << " is incompatible with that at round " << dec << i << "." << endl;
                    throw KeccakException("Incompatible states found in the trail.");
                }
            }
    }
    if (trail.stateAfterLastChiSpecified) {
        for(unsigned int z=0; z<laneSize; z++)
            for(unsigned int y=0; y<nrRowsAndColumns; y++) {
                const vector<RowValue>& values = diffChi[getRowFromSlice(trail.states.back()[z], y)].values;
                if (find(values.begin(), values.end(), getRowFromSlice(trail.stateAfterLastChi[z], y)) == values.end()) {
                    if (DC) trail.display(*DC, cerr);
                    cerr << "The state after the last \xCF\x87 is incompatible with that of the last round." << endl;
                    throw KeccakException("Incompatible states found in the trail.");
                }
            }
    }
}

void KeccakFDCLC::checkLCTrail(const Trail& trail, KeccakFPropagation *LC) const
{
    // Check weights
    unsigned int totalWeight = 0;
    unsigned int offsetIndex = (trail.firstStateSpecified ? 0 : 1);
    if ((!trail.firstStateSpecified) && (trail.weights.size() >= 1))
        totalWeight += trail.weights[0];
    for(unsigned int i=offsetIndex; i<trail.weights.size(); i++) {
        unsigned int weight = 0;
        for(unsigned int z=0; z<laneSize; z++)
            for(unsigned int y=0; y<nrRowsAndColumns; y++)
                weight += corrInvChi[getRowFromSlice(trail.states[i][z], y)].minWeight;
        if (weight != trail.weights[i]) {
            if (LC) trail.display(*LC, cerr);
            cerr << "The weight of state at round " << dec << i << " is incorrect; it should be " << weight << "." << endl;
            throw KeccakException("The weights in the trail are incorrect!");
        }
        totalWeight += weight;
    }
    if (totalWeight != trail.totalWeight) {
        if (LC) trail.display(*LC, cerr);
        cerr << "The total weight of the trail is incorrect; it should be " << totalWeight << "." << endl;
        throw KeccakException("The total weight in the trail is incorrect!");
    }

    // Check compatibility between consecutive states
    for(unsigned int i=1+offsetIndex; i<trail.states.size(); i++) {
        vector<SliceValue> stateAfterChi;
        lambda(trail.states[i], stateAfterChi, Dual);
        for(unsigned int z=0; z<laneSize; z++)
            for(unsigned int y=0; y<nrRowsAndColumns; y++) {
                const vector<RowValue>& values = corrInvChi[getRowFromSlice(trail.states[i-1][z], y)].values;
                if (find(values.begin(), values.end(), getRowFromSlice(stateAfterChi[z], y)) == values.end()) {
                    if (LC) trail.display(*LC, cerr);
                    cerr << "The state at round " << dec << i-1 << " is incompatible with that at round " << dec << i << "." << endl;
                    throw KeccakException("Incompatible states found in the trail.");
                }
            }
    }
    if (trail.stateAfterLastChiSpecified) {
        for(unsigned int z=0; z<laneSize; z++)
            for(unsigned int y=0; y<nrRowsAndColumns; y++) {
                const vector<RowValue>& values = corrInvChi[getRowFromSlice(trail.states.back()[z], y)].values;
                if (find(values.begin(), values.end(), getRowFromSlice(trail.stateAfterLastChi[z], y)) == values.end()) {
                    if (LC) trail.display(*LC, cerr);
                    cerr << "The state after the last \xCF\x87 is incompatible with that of the last round." << endl;
                    throw KeccakException("Incompatible states found in the trail.");
                }
            }
    }
}

void KeccakFDCLC::thetaTransEnvelope(vector<LaneValue>& state) const
{
    // We invert the order of bits in lanes and of lanes in planes
    LaneValue tmp, tmp2;
    for(unsigned int y=0; y<nrRowsAndColumns; y++) {
        tmp = state[index(1,y)];
        state[index(1,y)] = state[index(4,y)];
        state[index(4,y)] = tmp;
        tmp = state[index(2,y)];
        state[index(2,y)] = state[index(3,y)];
        state[index(3,y)] = tmp;
    }
    for(unsigned int x=0; x<nrRowsAndColumns; x++) {
        for(unsigned int y=0; y<nrRowsAndColumns; y++) {
            tmp = state[index(x,y)];
            tmp2 = 0;
            for(unsigned int z=0; z<laneSize; z++) {
                tmp2 ^= tmp&1;
                ROL(tmp,1);
                ROL(tmp2,-1);
            }
            state[index(x,y)] = tmp2;
        }
    }
}

unsigned int KeccakFDCLC::getThetaGap(const vector<LaneValue>& state) const
{
    vector<LaneValue> C(5);
    for(unsigned int x=0; x<5; x++) {
        C[x] = state[index(x,0)];
        for(unsigned int y=1; y<5; y++)
            C[x] ^= state[index(x,y)];
    }
    return getThetaGapFromParity(C);
}

unsigned int KeccakFDCLC::getThetaGapFromParity(const vector<LaneValue>& parities) const
{
    vector<LaneValue> D(5);
    getThetaEffectFromParity(parities, D);
    return
         (getHammingWeightLane(D[0])
        + getHammingWeightLane(D[1])
        + getHammingWeightLane(D[2])
        + getHammingWeightLane(D[3])
        + getHammingWeightLane(D[4]))/2;
}

string KeccakFDCLC::getName() const
{
    stringstream a;
    a << "KeccakF-" << dec << width;
    return a.str();
}
