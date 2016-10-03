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

#include "spongetree.h"

ShortLeafInterleavedSpongeTree::ShortLeafInterleavedSpongeTree(const SpongeFactory& aFactory, const int aD, const int aB)
    : factory(aFactory), absorbQueue(aB), D(aD), B(aB)
{
    for(int i=0; i<D; i++)
        leaves.push_back(factory.newSponge());
    final = factory.newSponge();
    C = 8*((final->getCapacity()+7)/8);
    leafIndex = 0;
    squeezing = false;
}

ShortLeafInterleavedSpongeTree::~ShortLeafInterleavedSpongeTree()
{
    for(unsigned int i=0; i<leaves.size(); i++) {
        delete leaves[i];
        leaves[i] = (Sponge*)0;
    }
    delete final;
}

void ShortLeafInterleavedSpongeTree::absorb(const UINT8 *input, unsigned int lengthInBits)
{
    unsigned int lengthInBytes = (lengthInBits+7)/8;
    vector<UINT8> inputAsVector(input, input+lengthInBytes);
    absorb(inputAsVector, lengthInBits);
}

void ShortLeafInterleavedSpongeTree::absorb(const vector<UINT8>& input, unsigned int lengthInBits)
{
    if (squeezing)
        throw SpongeException("The absorbing phase is over.");
    absorbQueue.append(input.begin(), lengthInBits);
    while(absorbQueue.firstBlockIsWhole()) {
        leaves[leafIndex]->absorb(absorbQueue.firstBlock(), B);
        absorbQueue.removeFirstBlock();
        leafIndex = (leafIndex+1)%D;
    }
}

void ShortLeafInterleavedSpongeTree::flushAndSwitchToSqueezingPhase()
{
    leaves[leafIndex]->absorb(absorbQueue.firstBlock(), absorbQueue.lastBlockSize());
    absorbQueue.clear();
    for(int i=0; i<D; i++) {
        leaves[i]->absorb((const UINT8*)"\x00", 1); // non-final bit
        vector<UINT8> chaining;
        leaves[i]->squeeze(chaining, C);
        final->absorb(chaining, C);
    }
    UINT8 parametersAndFinalBit[5];
    parametersAndFinalBit[0] = (UINT8)B;
    parametersAndFinalBit[1] = (UINT8)(B/256);
    parametersAndFinalBit[2] = (UINT8)(B/65536);
    parametersAndFinalBit[3] = (UINT8)(B/16777216);
    parametersAndFinalBit[4] = 0x00 | 0x20; // 5 bits set to zero (reserved for future use, for other layouts), plus final bit
    final->absorb(parametersAndFinalBit, 32+5+1);
    squeezing = true;
}

void ShortLeafInterleavedSpongeTree::squeeze(UINT8 *output, unsigned int desiredLengthInBits)
{
    if (!squeezing)
        flushAndSwitchToSqueezingPhase();
    final->squeeze(output, desiredLengthInBits);
}

void ShortLeafInterleavedSpongeTree::squeeze(vector<UINT8>& output, unsigned int desiredLengthInBits)
{
    if (!squeezing)
        flushAndSwitchToSqueezingPhase();
    final->squeeze(output, desiredLengthInBits);
}
