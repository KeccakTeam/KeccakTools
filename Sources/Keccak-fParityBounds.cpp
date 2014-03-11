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
#include "Keccak-fParityBounds.h"
#include "Keccak-fPositions.h"
#include "progress.h"
#include "translationsymmetry.h"

using namespace std;

unsigned int getBoundOfTotalWeightGivenTotalHammingWeight(const KeccakFPropagation& DCorLC, unsigned int totalHW)
{
    return DCorLC.getLowerBoundOnReverseWeightGivenHammingWeight(totalHW);
    // It has been explicitly checked that putting all the Hamming weight 
    // before λ always gives the minimum lower bound.
    // The check was done both for DC and LC.
}

unsigned int getLowerBoundTotalActiveRowsFromACandUOC(const KeccakFPropagation& DCorLC, 
    const vector<ColumnPosition>& xzAC, const vector<ColumnPosition>& xzUOC)
{
    unsigned int activeRows = 0;
    vector<vector<bool> > rowTakenLeft(5, vector<bool>(DCorLC.laneSize, false));
    vector<vector<bool> > rowTakenRight(5, vector<bool>(DCorLC.laneSize, false));
    
    for(unsigned int i=0; i<xzAC.size(); i++) {
        unsigned int x = xzAC[i].x;
        unsigned int z = xzAC[i].z;
        for(unsigned int y=0; y<5; y++) {
            BitPosition left(x, y, z);
            DCorLC.reverseRhoPiBeforeTheta(left);
            BitPosition right(x, y, z);
            DCorLC.directRhoPiAfterTheta(right);
            if ((!rowTakenLeft[left.y][left.z]) && (!rowTakenRight[right.y][right.z])) {
                activeRows++;
                rowTakenLeft[left.y][left.z] = true;
                rowTakenRight[right.y][right.z] = true;
            }
        }
    }
    
    for(unsigned int i=0; i<xzUOC.size(); i++) {
        unsigned int x = xzUOC[i].x;
        unsigned int z = xzUOC[i].z;
        bool takenLeft = false;
        bool takenRight = false;
        for(unsigned int y=0; y<5; y++) {
            BitPosition left(x, y, z);
            DCorLC.reverseRhoPiBeforeTheta(left);
            BitPosition right(x, y, z);
            DCorLC.directRhoPiAfterTheta(right);
            takenLeft |= rowTakenLeft[left.y][left.z];
            takenRight |= rowTakenRight[right.y][right.z];
            rowTakenLeft[left.y][left.z] = true;
            rowTakenRight[right.y][right.z] = true;
        }
        if (!takenLeft)
            activeRows++;
        if (!takenRight)
            activeRows++;
        for(unsigned int y=0; y<5; y++) {
            BitPosition left(x, y, z);
            DCorLC.reverseRhoPiBeforeTheta(left);
            BitPosition  right(x, y, z);
            DCorLC.directRhoPiAfterTheta(right);
            if (takenLeft)
                rowTakenLeft[left.y][left.z] = true;
            if (takenRight)
                rowTakenRight[right.y][right.z] = true;
        }
    }

    return activeRows;
}

unsigned int getLowerBoundTotalActiveRows(const KeccakFPropagation& DCorLC, 
    const vector<RowValue>& C, const vector<RowValue>& D)
{
    vector<ColumnPosition> xzAC, xzUOC;
    for(unsigned int x=0; x<5; x++)
    for(unsigned int z=0; z<DCorLC.laneSize; z++) {
        bool odd = (getBit(C, x, z) != 0);
        bool affected = (getBit(D, x, z) != 0);
        if (affected)
            xzAC.push_back(ColumnPosition(x, z));
        else  {
            if (odd)
                xzUOC.push_back(ColumnPosition(x, z));
        }
    }
    return getLowerBoundTotalActiveRowsFromACandUOC(DCorLC, xzAC, xzUOC);
}

string Run::display() const
{
    stringstream str;
    str << "[" << dec << tStart;
    if (length == 1)
        str << "]";
    else
        str << "-" << (tStart+length-1) << "]";
    return str.str();
}

string ParityAsRuns::display() const
{
    string result;
    for(unsigned int i=0; i<runs.size(); i++) {
        result += runs[i].display();
        if (i < (runs.size()-1))
            result += " ";
    }
    return result;
}

void ParityAsRuns::toParityAndParityEffect(const KeccakFPropagation& DCorLC, vector<RowValue>& C, vector<RowValue>& D) const
{
    C.assign(DCorLC.laneSize, 0);
    D.assign(DCorLC.laneSize, 0);
    for(unsigned int i=0; i<runs.size(); i++) {
        unsigned int x, z;
        DCorLC.getXandZfromT(DCorLC.translateAlongXinT(runs[i].tStart), x, z);
        setBitToOne(D, x, z);
        DCorLC.getXandZfromT(DCorLC.translateAlongXinT(runs[i].tStart + runs[i].length), x, z);
        setBitToOne(D, x, z);
        for(unsigned int t=runs[i].tStart; t<runs[i].tStart+runs[i].length; t++) {
            unsigned int x, z;
            DCorLC.getXandZfromT(t, x, z);
            setBitToOne(C, x, z);
        }
    }
}

unsigned int ParityAsRuns::getLowerBoundTotalHammingWeight(const KeccakFPropagation& DCorLC) const
{
    vector<bool> affected(DCorLC.laneSize*5, false);
    for(unsigned int i=0; i<runs.size(); i++) {
        affected[DCorLC.translateAlongXinT(runs[i].tStart)] = true;
        affected[DCorLC.translateAlongXinT(runs[i].tStart + runs[i].length)] = true;
    }
    unsigned int total = 5*2*runs.size();
    for(unsigned int i=0; i<runs.size(); i++) {
        for(unsigned int t=runs[i].tStart; t<runs[i].tStart+runs[i].length; t++) {
            if (!affected[t % (DCorLC.laneSize*5)])
                total += 2;
        }
    }
    return total;
}

unsigned int ParityAsRuns::getLowerBoundTotalActiveRowsUsingOnlyAC(const KeccakFPropagation& DCorLC) const
{
    vector<ColumnPosition> xzAC, xzUOC;
    for(unsigned int i=0; i<runs.size(); i++) {
        unsigned int x, z;
        DCorLC.getXandZfromT(DCorLC.translateAlongXinT(runs[i].tStart), x, z);
        xzAC.push_back(ColumnPosition(x, z));
        DCorLC.getXandZfromT(DCorLC.translateAlongXinT(runs[i].tStart + runs[i].length), x, z);
        xzAC.push_back(ColumnPosition(x, z));
    }
    return getLowerBoundTotalActiveRowsFromACandUOC(DCorLC, xzAC, xzUOC);
}

unsigned int ParityAsRuns::getLowerBoundTotalActiveRows(const KeccakFPropagation& DCorLC) const
{
    vector<ColumnPosition> xzAC, xzUOC;
    vector<bool> affected(DCorLC.laneSize*5, false);
    for(unsigned int i=0; i<runs.size(); i++) {
        unsigned int x, z;
        DCorLC.getXandZfromT(DCorLC.translateAlongXinT(runs[i].tStart), x, z);
        xzAC.push_back(ColumnPosition(x, z));
        DCorLC.getXandZfromT(DCorLC.translateAlongXinT(runs[i].tStart + runs[i].length), x, z);
        xzAC.push_back(ColumnPosition(x, z));
        affected[DCorLC.translateAlongXinT(runs[i].tStart)] = true;
        affected[DCorLC.translateAlongXinT(runs[i].tStart + runs[i].length)] = true;
    }
    for(unsigned int i=0; i<runs.size(); i++) {
        for(unsigned int t=runs[i].tStart; t<runs[i].tStart+runs[i].length; t++) {
            if (!affected[t % (DCorLC.laneSize*5)]) {
                unsigned int x, z;
                DCorLC.getXandZfromT(t, x, z);
                xzUOC.push_back(ColumnPosition(x, z));
            }
        }
    }
    return getLowerBoundTotalActiveRowsFromACandUOC(DCorLC, xzAC, xzUOC);
}

void lookForRunsBelowTargetWeight(const KeccakFPropagation& DCorLC, ostream& out,
    unsigned int targetWeight, ParityAsRuns& parity, ProgressMeter& progress, bool verbose)
{
    unsigned int lowerBound;
    unsigned int weightBoundBasedOnTotalHammingWeight = getBoundOfTotalWeightGivenTotalHammingWeight(DCorLC,
        parity.getLowerBoundTotalHammingWeight(DCorLC));
    if (weightBoundBasedOnTotalHammingWeight <= targetWeight) {
        unsigned int minActiveRows = parity.getLowerBoundTotalActiveRowsUsingOnlyAC(DCorLC);
        lowerBound = max(minActiveRows*2, weightBoundBasedOnTotalHammingWeight);
    }
    else
        lowerBound = weightBoundBasedOnTotalHammingWeight;
    if (lowerBound <= targetWeight) {
        unsigned int thisOneLowerBound = parity.getLowerBoundTotalActiveRows(DCorLC)*2;
        if (thisOneLowerBound <= targetWeight) {
            vector<RowValue> C, D;
            parity.toParityAndParityEffect(DCorLC, C, D);
            unsigned int thisOneLowerBoundAgain = getLowerBoundTotalActiveRows(DCorLC, C, D)*2;
            if (thisOneLowerBoundAgain <= targetWeight) {
                if (verbose) {
                    displayParity(cout, C, D);
                    cout << "Lower bound = " << dec << max(thisOneLowerBound, thisOneLowerBoundAgain) << endl;
                    cout << endl;
                }
                vector<RowValue> Cmin;
                getSymmetricMinimum(C, Cmin);
                writeParity(out, Cmin);
            }
        }
        progress.stack("Adding runs to "+parity.display());
        for(unsigned int tStart=parity.runs.back().tStart+parity.runs.back().length+1; 
                tStart<DCorLC.laneSize*5; tStart++) {
            unsigned int maxLength = DCorLC.laneSize*5-1-tStart+parity.runs[0].tStart;
            for(unsigned int length=1; length<=maxLength; length++) {
                Run run;
                run.tStart = tStart;
                run.length = length;
                parity.runs.push_back(run);
                lookForRunsBelowTargetWeight(DCorLC, out, targetWeight, parity, progress, verbose);
                parity.runs.pop_back();
                ++progress;
            }
        }
        progress.unstack();
    }
}

void lookForRunsBelowTargetWeight(const KeccakFPropagation& DCorLC, ostream& out,  unsigned int targetWeight, bool verbose)
{
    ProgressMeter progress;
    progress.stack("Initial run starting point", 5);
    for(unsigned int tStart=0; tStart<5; tStart++) {
        progress.stack("Initial run length", DCorLC.laneSize*5-1);
        for(unsigned int length=1; length<=DCorLC.laneSize*5-1; length++) {
            ParityAsRuns parity;
            Run run;
            run.tStart = tStart;
            run.length = length;
            parity.runs.push_back(run);
            lookForRunsBelowTargetWeight(DCorLC, out, targetWeight, parity, progress, verbose);
            ++progress;
        }
        progress.unstack();
        ++progress;
    }
    progress.unstack();
}    
