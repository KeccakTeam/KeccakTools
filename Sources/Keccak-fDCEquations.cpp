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

#include "Keccak-fDCEquations.h"

KeccakFDCEquations::KeccakFDCEquations(unsigned int aWidth)
    : KeccakFDCLC(aWidth)
{
}

void KeccakFDCEquations::buildDCTrailFromPair(vector<SliceValue>& a1, vector<SliceValue>& a2, Trail& trail) const
{
    KeccakFPropagation DC(*this, KeccakFPropagation::DC);

    vector<LaneValue> state1;
    fromSlicesToLanes(a1, state1);
    vector<LaneValue> state2;
    fromSlicesToLanes(a2, state2);

    for(unsigned int i=0; i<nrRounds+1; i++) {
        lambda(state1, KeccakFDCLC::Straight);
        lambda(state2, KeccakFDCLC::Straight);

        vector<LaneValue> stateDiff(state1);
        for(unsigned int xy=0; xy<nrRowsAndColumns*nrRowsAndColumns; xy++)
            stateDiff[xy] ^= state2[xy];
        vector<SliceValue> aDiff;
        fromLanesToSlices(stateDiff, aDiff);
        trail.append(aDiff, DC.getWeight(aDiff));

        chi(state1);       chi(state2);
        iota(state1, i);   iota(state2, i);

        if (i == (nrRounds-1)) {
            fromLanesToSlices(state1, a1);
            fromLanesToSlices(state2, a2);
        }
    }
}

void KeccakFDCEquations::getDCEquations(RowValue diffIn, RowValue diffOut,
    const vector<SymbolicBit>& inputVariables, vector<SymbolicBit>& inputRelations) const
{
    RowValue diffOutXorChiDiffIn = diffOut ^ chiOnRow(diffIn);
    // Note: when diffIn=11111, 5 equations are generated, the fifth being redundant with the 4 first
    for(unsigned int i=0; i<nrRowsAndColumns; i++) {
        RowValue t = translateRowSafely(diffIn, -(int)i);
        bool addOne = ((diffOutXorChiDiffIn >> i) & 1) != 0;
        // in:  .X-
        // out: a_i+2 = ...
        if ((t & 0x6) == 0x2) {
            SymbolicBit relation = inputVariables[(i+2)%nrRowsAndColumns];
            if (addOne)
                relation.complement();
            inputRelations.push_back(relation);
        }
        // in:  .XX
        // out: a_i+1 + a_i+2 = ...
        if ((t & 0x6) == 0x6) {
            SymbolicBit relation = inputVariables[(i+1)%nrRowsAndColumns];
            relation.add(inputVariables[(i+2)%nrRowsAndColumns]);
            if (addOne)
                relation.complement();
            inputRelations.push_back(relation);
        }
        // in:  --X
        // out: a_i+1 = ...
        if ((t & 0x7) == 0x4) {
            SymbolicBit relation = inputVariables[(i+1)%nrRowsAndColumns];
            if (addOne)
                relation.complement();
            inputRelations.push_back(relation);
        }
    }
}

void KeccakFDCEquations::getDCEquations(const vector<SliceValue>& diffIn,
    const vector<SliceValue>& diffOut, const vector<SymbolicLane>& input,
    vector<SymbolicBit>& inputRelations) const
{
    for(unsigned int z=0; z<laneSize; z++)
    for(unsigned int y=0; y<nrRowsAndColumns; y++) {
        RowValue diffInRow = getRowFromSlice(diffIn[z], y);
        if (diffInRow != 0) {
            RowValue diffOutRow = getRowFromSlice(diffOut[z], y);
            vector<SymbolicBit> inputVariables(nrRowsAndColumns);
            for(unsigned int x=0; x<nrRowsAndColumns; x++)
                inputVariables[x] = input[index(x, y)].values[z];
            getDCEquations(diffInRow, diffOutRow, inputVariables, inputRelations);
        }
    }
}

void KeccakFDCEquations::genDCEquations(ostream& fout, const Trail& trail, bool forSage) const
{
    char input = 'A';
    char output = 'B';
    for(unsigned int r=0; r<trail.states.size()-1; r++) {
        string inputName(1, input);
        string outputName(1, output);

        fout << "// Round " << dec << r << endl;
        fout << "// Conditions at input of \xCF\x87" << endl;
        vector<SliceValue> stateAfterChi;
        lambda(trail.states[r+1], stateAfterChi, KeccakFDCLC::Inverse);
        vector<SymbolicLane> variables;
        KeccakFEquations::initializeState(variables, inputName, laneSize);
        vector<SymbolicBit> relations;
        getDCEquations(trail.states[r], stateAfterChi, variables, relations);
        for(unsigned int i=0; i<relations.size(); i++) {
            fout << relations[i].value;
            if (forSage)
                fout << ",";
            else
                fout << " = 0";
            fout << endl;
        }

        fout << "// Linking to next round: " << outputName << " = \xCF\x80(\xCF\x81(\xCE\xB8(\xCE\xB9(\xCF\x87(" << inputName << ")))))" << endl;
        vector<SymbolicLane> variablesAfterRound;
        KeccakFEquations::initializeState(variablesAfterRound, inputName, laneSize);
        chi(variablesAfterRound);
        iota(variablesAfterRound, r);
        theta(variablesAfterRound);
        rho(variablesAfterRound);
        pi(variablesAfterRound);
        displayEquations(fout, variablesAfterRound, outputName, forSage);

        input = output;
        output++;
    }
}

void KeccakFDCEquations::displayEquations(ostream& fout, const vector<SymbolicLane>& state, const string& prefixOutput, bool forSage) const
{
    for(unsigned int y=0; y<nrRowsAndColumns; y++)
    for(unsigned int x=0; x<nrRowsAndColumns; x++)
    for(unsigned int z=0; z<laneSize; z++) {
        string outputBit = bitName(prefixOutput, x, y, z);
        fout << outputBit;
        if (forSage)
            fout << " + ";
        else
            fout << " = ";
        fout << state[index(x,y)].values[z].value;
        if (forSage)
            fout << ", ";
        fout << endl;
    }
}
