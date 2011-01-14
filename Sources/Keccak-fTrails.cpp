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

#include <fstream>
#include "Keccak-fPropagation.h"
#include "Keccak-fTrails.h"

using namespace std;

Trail::Trail()
    : totalWeight(0)
{
}

Trail::Trail(istream& fin)
{
    load(fin);
}

void Trail::append(const vector<SliceValue>& state, unsigned int weight)
{
    states.push_back(state);
    weights.push_back(weight);
    totalWeight += weight;
}

void Trail::prepend(const vector<SliceValue>& state, unsigned int weight)
{
    states.insert(states.begin(), state);
    weights.insert(weights.begin(), weight);
    totalWeight += weight;
}

void Trail::display(const KeccakFPropagation& DCorLC, ostream& fout) const
{
    if (states.size() == 0) {
        fout << "This trail is empty." << endl;
        return;
    }

    fout << dec << states.size() << "-round ";
    if (DCorLC.getPropagationType() == KeccakFPropagation::DC)
        fout << "differential ";
    else
        fout << "linear ";
    fout << "trail of total weight " << dec << totalWeight << endl;

    bool thetaJustAfterChi = DCorLC.isThetaJustAfterChi();
    vector<vector<SliceValue> > allAfterPreviousChi;
    vector<vector<SliceValue> > allBeforeTheta;

    string trailKernelType;
    vector<unsigned int> thetaGaps, activeRows;
    for(unsigned int i=0; i<states.size(); i++) {
        vector<SliceValue> stateAfterChi;
        DCorLC.reverseLambda(states[i], stateAfterChi);
        allAfterPreviousChi.push_back(stateAfterChi);
        vector<SliceValue> stateBeforeTheta;
        DCorLC.directLambdaBeforeTheta(stateAfterChi, stateBeforeTheta);
        allBeforeTheta.push_back(stateBeforeTheta);

        bool kernel = true;
        for(unsigned int z=0; z<stateBeforeTheta.size(); z++) {
            if (getParity(stateBeforeTheta[z]) != 0) {
                kernel = false;
                break;
            }
        }
        trailKernelType += (kernel ? 'k' : 'N');

        vector<LaneValue> stateBeforeThetaLanes;
        fromSlicesToLanes(stateBeforeTheta, stateBeforeThetaLanes);
        thetaGaps.push_back(DCorLC.parent.getThetaGap(stateBeforeThetaLanes));

        activeRows.push_back(getNrActiveRows(states[i]));
    }

    fout << "* Profile related to \xCF\x87:" << endl;
    fout << "Propagation weights: ";
    for(unsigned int i=0; i<weights.size(); i++) {
        fout.width(4); fout.fill(' ');
        fout << dec << weights[i];
    }
    fout << endl;
    fout << "Active rows:         ";
    for(unsigned int i=0; i<activeRows.size(); i++) {
        fout.width(4); fout.fill(' ');
        fout << dec << activeRows[i];
    }
    fout << endl;

    fout << "* Profile related to \xCE\xB8:" << endl;
    fout << "Gaps:   ";
    for(unsigned int i=0; i<thetaGaps.size(); i++) {
        fout.width(4); fout.fill(' ');
        fout << dec << thetaGaps[i];
    }
    fout << endl;
    fout << "Kernel: ";
    for(unsigned int i=0; i<trailKernelType.size(); i++) {
        fout << "   " << trailKernelType[i];
    }
    fout << endl;

    unsigned int minReverseWeight = 0;
    for(unsigned int z=0; z<allAfterPreviousChi[0].size(); z++)
    for(unsigned int y=0; y<nrRowsAndColumns; y++) {
        RowValue row = getRowFromSlice(allAfterPreviousChi[0][z], y);
        minReverseWeight += DCorLC.reverseRowOutputListPerInput[row].minWeight;
    }
    fout << "Previous round would have weight at least " << dec << minReverseWeight << endl;
    for(unsigned int i=0; i<states.size(); i++) {
        fout << "Round " << dec << i << " (weight " << weights[i];
        if (thetaJustAfterChi)
            fout << ", \xCE\xB8-gap " << thetaGaps[i];
        fout << ") after previous \xCF\x87";
        if (!thetaJustAfterChi)
            fout << ", then before \xCE\xB8 of gap " << thetaGaps[i];
        fout << ", then before \xCF\x87";
        fout << " (" << activeRows[i] << " active rows):" << endl;
        if (thetaJustAfterChi)
            displayStates(fout, allAfterPreviousChi[i], true, states[i], false);
        else
            displayStates(fout, allAfterPreviousChi[i], false, allBeforeTheta[i], true, states[i], false);
    }
}

void Trail::save(ostream& fout) const
{
    fout << hex;
    if (states.size() > 0)
        fout << states[0].size() << " ";
    else
        fout << 0 << " ";
    fout << totalWeight << " ";
    fout << weights.size() << " ";
    for(unsigned int i=0; i<weights.size(); i++)
        fout << weights[i] << " ";
    for(unsigned int i=0; i<states.size(); i++)
    for(unsigned int j=0; j<states[i].size(); j++)
        fout << states[i][j] << " ";
    fout << endl;
}

void Trail::load(istream& fin)
{
    unsigned int laneSize = 0;
    fin >> hex;
    fin >> laneSize;
    if ((laneSize == 0) || (fin.eof()))
        throw TrailException();
    fin >> totalWeight;
    unsigned int size;
    fin >> size;
    weights.resize(size);
    for(unsigned int i=0; i<size; i++)
        fin >> weights[i];
    states.clear();
    for(unsigned int i=0; i<size; i++) {
        vector<SliceValue> state(laneSize);
        for(unsigned int j=0; j<laneSize; j++)
            fin >> state[j];
        states.push_back(state);    
    }
}

void Trail::append(const Trail& otherTrail)
{
    for(unsigned int i=0; i<otherTrail.weights.size(); i++)
        append(otherTrail.states[i], otherTrail.weights[i]);
}

UINT64 Trail::produceHumanReadableFile(const KeccakFPropagation& DCorLC, const string& fileName, bool verbose, unsigned int maxWeight)
{
    string fileName2 = fileName+".txt";
    ofstream fout(fileName2.c_str());
    if (verbose)
        cout << "Writing " << fileName2 << flush;
    UINT64 count = DCorLC.displayTrailsAndCheck(fileName, fout, maxWeight);
    if (verbose)
        cout << endl;
    return count;
}
