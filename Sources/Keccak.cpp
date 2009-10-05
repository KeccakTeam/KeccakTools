/*
Tools for the Keccak sponge function family.
Authors: Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche

This code is hereby put in the public domain. It is given as is,
without any guarantee.

For more information, feedback or questions, please refer to our website:
http://keccak.noekeon.org/
*/

#include <iostream>
#include <sstream>
#include "Keccak.h"

using namespace std;

Keccak::Keccak(unsigned int aRate, unsigned int aCapacity, unsigned char aDiversifier)
    : Sponge(new KeccakF(aRate+aCapacity), aCapacity), diversifier(aDiversifier)
{
    if ((rate % 8) != 0)
        throw KeccakException("The rate must be a multiple of 8 bits.");
}

Keccak::~Keccak()
{
    delete f;
}

string Keccak::getDescription() const
{
    stringstream a;
    a << "Keccak[r=" << dec << rate << ", c=" << dec << capacity 
        << ", d=" << dec << (int)diversifier << "]";
    return a.str();
}

void Keccak::pad()
{
    UINT8 outerPadding[2];

    padCurrentBlock(8);
    outerPadding[0] = diversifier;
    outerPadding[1] = rate/8;
    absorb(outerPadding, 16);
    Sponge::pad();
}
