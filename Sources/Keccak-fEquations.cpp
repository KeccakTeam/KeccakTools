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
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include "Keccak-fEquations.h"

using namespace std;

SymbolicBit::SymbolicBit()
    : value("0"), containsAddition(false)
{
}

SymbolicBit::SymbolicBit(bool aValue)
    : value(aValue ? "1" : "0"), containsAddition(false)
{
}

SymbolicBit::SymbolicBit(const string& aValue, bool aContainsAddition)
    : value(aValue), containsAddition(aContainsAddition)
{
}

void SymbolicBit::complement()
{
    if (value == "0")
        value = "1";
    else {
        value += " + 1";
        containsAddition = true;
    }
}

void SymbolicBit::add(const SymbolicBit& a)
{
    if (value == "0")
        value = a.value;
    else if (a.value != "0") {
        value += " + ";
        value += a.value;
        containsAddition = true;
    }
}

void SymbolicBit::multiply(const SymbolicBit& a)
{
    string temp;
    if (containsAddition)
        temp = "("+value+")";
    else
        temp = value;
    if (a.containsAddition)
        value = temp+"*("+a.value+")";
    else
        value = temp+"*"+a.value;
    containsAddition = false;
}


SymbolicLane::SymbolicLane()
{
}

SymbolicLane::SymbolicLane(LaneValue aValues)
{
    for(unsigned int i=0; i<64; i++)
        values.push_back(SymbolicBit((aValues & ((LaneValue)1 << i)) != 0));
}

SymbolicLane::SymbolicLane(unsigned int laneSize, const string& prefixSymbol)
{
    for(unsigned int z=0; z<laneSize; z++)
        values.push_back(SymbolicBit(KeccakF::buildBitName(prefixSymbol, laneSize, z)));
}

void SymbolicLane::ROL(int offset, unsigned int laneSize)
{
    if (laneSize <= values.size()) {
        offset %= laneSize;
        if (offset < 0) offset += laneSize;
        vector<SymbolicBit> newValues;
        for(unsigned int i=0; i<laneSize; i++)
            newValues.push_back(values[(laneSize+i-offset)%laneSize]);
        values = newValues;
    }
    else
        throw KeccakException("Incorrect usage of SymbolicLane::ROL.");
}

SymbolicLane operator~(const SymbolicLane& lane)
{
    SymbolicLane result = lane;
    for(unsigned int i=0; i<result.values.size(); i++)
        result.values[i].complement();
    return result;
}

SymbolicLane operator^(const SymbolicLane& a, LaneValue b)
{
    SymbolicLane result = a;
    for(unsigned int i=0; i<result.values.size(); i++)
        if (b & ((LaneValue)1 << i))
            result.values[i].complement();
    return result;
}

SymbolicLane operator^(const SymbolicLane& a, const SymbolicLane& b)
{
    SymbolicLane result = a;
    for(unsigned int i=0; i<result.values.size(); i++)
        if (i < b.values.size())
            result.values[i].add(b.values[i]);
    return result;
}

SymbolicLane operator&(const SymbolicLane& a, const SymbolicLane& b)
{
    SymbolicLane result = a;
    for(unsigned int i=0; i<result.values.size(); i++)
        if (i < b.values.size())
            result.values[i].multiply(b.values[i]);
    return result;
}

SymbolicLane& SymbolicLane::operator^=(const SymbolicLane& b)
{
    for(unsigned int i=0; i<values.size(); i++)
        if (i < b.values.size())
            values[i].add(b.values[i]);
    return *this;
}


KeccakFEquations::KeccakFEquations(unsigned int aWidth)
    : KeccakF(aWidth)
{
}

void KeccakFEquations::genEquations(ostream& fout, const vector<SymbolicLane>& state, const string& prefixOutput, bool forSage) const
{
    for(unsigned int y=0; y<5; y++)
    for(unsigned int x=0; x<5; x++) {
        for(unsigned int z=0; z<laneSize; z++) {
            string outputBit = bitName(prefixOutput, x, y, z);
            if (forSage)
                fout << "    '";
            fout << outputBit;
            if (forSage)
                fout << " + ";
            else
                fout << " = ";
            fout << state[index(x,y)].values[z].value;
            if (forSage)
                fout << "',";
            fout << endl;
        }
    }
}

void KeccakFEquations::initializeState(vector<SymbolicLane>& state, const string& prefix) const
{
    initializeState(state, prefix, laneSize);
}

void KeccakFEquations::initializeState(vector<SymbolicLane>& state, const string& prefix,
    unsigned int laneSize)
{
    state.resize(25);
    for(unsigned int x=0; x<5; x++)
    for(unsigned int y=0; y<5; y++)
        state[index(x,y)] = SymbolicLane(laneSize, laneName(prefix, x, y));
}

void KeccakFEquations::genRoundEquations(ostream& fout, int startRoundIndex, unsigned nrRounds, bool forSage) const
{
    char input = 'A';
    char output = 'B';
    for(unsigned int i=startRoundIndex; i<startRoundIndex+nrRounds; i++) {
        if (!forSage)
            fout << "// --- Round " << dec << i << endl;
        string inputName(1, input);
        string outputName(1, output);
        vector<SymbolicLane> state;
        initializeState(state, inputName);
        round(state, i);
        genEquations(fout, state, outputName, forSage);
        input = output;
        output++;
    }
}

void KeccakFEquations::genComponentEquations(ostream& fout, const string& prefixInput, const string& prefixOutput) const
{
    { // θ
        fout << "// --- \xCE\xB8" << endl;
        vector<SymbolicLane> state;
        initializeState(state, prefixInput);
        theta(state);
        genEquations(fout, state, prefixOutput);
    }
    { // θ^-1
        fout << "// --- \xCE\xB8^-1" << endl;
        vector<SymbolicLane> state;
        initializeState(state, prefixInput);
        inverseTheta(state);
        genEquations(fout, state, prefixOutput);
    }
    { // ρ then π
        fout << "// --- \xCF\x81 then \xCF\x80" << endl;
        vector<SymbolicLane> state;
        initializeState(state, prefixInput);
        rho(state);
        pi(state);
        genEquations(fout, state, prefixOutput);
    }
    { // π^-1 then ρ^-1
        fout << "// --- \xCF\x80^-1 then \xCF\x81^-1" << endl;
        vector<SymbolicLane> state;
        initializeState(state, prefixInput);
        inversePi(state);
        inverseRho(state);
        genEquations(fout, state, prefixOutput);
    }
    { // χ
        fout << "// --- \xCF\x87" << endl;
        vector<SymbolicLane> state;
        initializeState(state, prefixInput);
        chi(state);
        genEquations(fout, state, prefixOutput);
    }
    { // χ^-1
        fout << "// --- \xCF\x87^-1" << endl;
        vector<SymbolicLane> state;
        initializeState(state, prefixInput);
        inverseChi(state);
        genEquations(fout, state, prefixOutput);
    }
    for(unsigned int i=0; i<roundConstants.size(); i++) { // ι
        fout << "// --- \xCE\xB9 for round " << dec << i << endl;
        vector<SymbolicLane> state;
        initializeState(state, prefixInput);
        iota(state, i);
        genEquations(fout, state, prefixOutput);
    }
}

void KeccakFEquations::genAbsoluteValuesBeforeChi(ostream& fout, const vector<LaneValue>& input, const string& prefix, int startRoundIndex, unsigned nrRounds) const
{
    vector<LaneValue> state(input);
    for(unsigned int i=startRoundIndex; i<startRoundIndex+nrRounds; i++) {
        theta(state);
        rho(state);
        pi(state);

        fout << "// Round " << dec << i << endl;
        for(unsigned int y=0; y<5; y++)
        for(unsigned int x=0; x<5; x++)
        for(unsigned int z=0; z<laneSize; z++)
            fout << bitName(prefix, x, y, z) << " = " << ((state[index(x, y)] >> z) & 1) << endl;

        chi(state);
        iota(state, i);
    }
}
