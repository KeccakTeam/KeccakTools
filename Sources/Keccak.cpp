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

#include <iostream>
#include <sstream>
#include "Keccak.h"

using namespace std;

Keccak::Keccak(unsigned int aRate, unsigned int aCapacity)
    : Sponge(new KeccakF(aRate+aCapacity), new MultiRatePadding(), aRate)
{
}

Keccak::~Keccak()
{
    delete f;
    delete pad;
}

string Keccak::getDescription() const
{
    stringstream a;
    a << "Keccak[r=" << dec << rate << ", c=" << dec << capacity << "]";
    return a.str();
}

ReducedRoundKeccak::ReducedRoundKeccak(unsigned int aRate, unsigned int aCapacity, int aStartRoundIndex, unsigned int aNrRounds)
    : Sponge(new KeccakFanyRounds(aRate+aCapacity, aStartRoundIndex, aNrRounds), new MultiRatePadding(), aRate),
    nrRounds(aNrRounds),
    startRoundIndex(aStartRoundIndex)
{
}

ReducedRoundKeccak::~ReducedRoundKeccak()
{
    delete f;
    delete pad;
}

string ReducedRoundKeccak::getDescription() const
{
    stringstream a;
    a << "Keccak[r=" << dec << rate << ", c=" << dec << capacity << ", " << dec << nrRounds << " ";
    a << ((nrRounds > 1) ? "rounds" : "round");
    a << " from " << dec << startRoundIndex << " to " << startRoundIndex+nrRounds-1 << "]";
    return a.str();
}
