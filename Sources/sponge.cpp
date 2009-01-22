/*
Tools for the Keccak sponge function family.
Authors: Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche

This code is hereby put in the public domain. It is given as is, 
without any guarantee.

For more information, feedback or questions, please refer to our website:
http://keccak.noekeon.org/
*/

#include <sstream>
#include <vector>
#include "sponge.h"

using namespace std;

SpongeException::SpongeException(const string& aReason)
    : reason(aReason)
{
}


Sponge::Sponge(Transformation *aF, unsigned int aCapacity)
    : f(aF), capacity(aCapacity), squeezing(false), validInput(false), bitsInCurrentBlock(0)
{
    unsigned int width = f->getWidth();
    if (capacity >= width)
        throw SpongeException("The requested capacity is too large when using this function.");
    rate = width-capacity;
    if ((rate % 8) != 0)
        throw SpongeException("The rate must be a multiple of 8.");
    state.reset(new UINT8[(width+7)/8]);
    for(unsigned int i=0; i<(width+7)/8; i++)
        state.get()[i] = 0;
}

unsigned int Sponge::getCapacity()
{
    return capacity;
}

unsigned int Sponge::getRate()
{
    return rate;
}

void Sponge::absorb(const UINT8 *input, unsigned int lengthInBits)
{
    if ((bitsInCurrentBlock % 8) != 0)
        throw SpongeException("The length of the input must be a multiple of 8, except for the last call to absorb().");
    if (squeezing)
        throw SpongeException("The absorbing phase is over.");
    unsigned int i = 0;
    while(i < lengthInBits) {
        unsigned int partialBlock = lengthInBits - i;
        if (partialBlock+bitsInCurrentBlock > rate)
            partialBlock = rate-bitsInCurrentBlock;
        for(unsigned int j=0; j<partialBlock; j+=8)
            currentBlock.push_back(input[(i+j)/8]);
        bitsInCurrentBlock += partialBlock;
        if (bitsInCurrentBlock == rate)
            absorbCurrentBlock();
        i += partialBlock;
    }
}

void Sponge::padCurrentBlock(unsigned int desiredLengthFactor)
{
    if (desiredLengthFactor == 0)
        throw SpongeException("The desired length factor must be greater than zero.");
    if ((desiredLengthFactor % 8) != 0)
        throw SpongeException("The desired length factor must be a multiple of 8.");
    if ((rate % desiredLengthFactor) != 0)
        throw SpongeException("The desired length factor must divide the rate.");
    if ((bitsInCurrentBlock % 8) != 0) {
        // Here the bits are numbered from 0=LSB to 7=MSB
        UINT8 padByte = 1 << (bitsInCurrentBlock % 8);
        currentBlock.back() |= padByte;
        bitsInCurrentBlock++;
        if ((bitsInCurrentBlock % 8) != 0) {
            UINT8 maskByte = (UINT8)0xFF << (bitsInCurrentBlock % 8);
            currentBlock.back() &= ~maskByte;
            bitsInCurrentBlock += 8-(bitsInCurrentBlock % 8);
        }
    }
    else {
        currentBlock.push_back(0x01);
        bitsInCurrentBlock += 8;
    }
    while((bitsInCurrentBlock % desiredLengthFactor) != 0) {
        currentBlock.push_back(0x00);
        bitsInCurrentBlock += 8;
    }
    if (bitsInCurrentBlock == rate)
        absorbCurrentBlock();
}

void Sponge::pad()
{
    padCurrentBlock(rate);
}

void Sponge::squeeze(UINT8 *output, unsigned int desiredLengthInBits)
{
    if ((desiredLengthInBits % 8) != 0)
        throw SpongeException("The desired output length must be a multiple of 8.");
    if (!squeezing)
        flushAndSwitchToSqueezingPhase();
    unsigned int i = 0;
    while(i < desiredLengthInBits) {
        if (bitsInCurrentBlock == 0)
            squeezeIntoCurrentBlock();
        unsigned int partialBlock = desiredLengthInBits-i;
        if (partialBlock > bitsInCurrentBlock)
            partialBlock = bitsInCurrentBlock;
        for(unsigned int j=0; j<partialBlock; j+=8) {
            output[(i+j)/8] = currentBlock.front();
            currentBlock.pop_front();
        }
        bitsInCurrentBlock -= partialBlock;
        i += partialBlock;
    }
}

void Sponge::absorbCurrentBlock()
{
    // bitsInCurrentBlock is assumed to be equal to rate
    for(vector<UINT8>::size_type i=0; i<currentBlock.size(); ++i)
        state.get()[i] ^= currentBlock[i];
    validInput = false;
    for(vector<UINT8>::size_type i=0; i<currentBlock.size(); ++i)
        if (currentBlock[i] != 0) {
            validInput = true;
            break;
        }
    currentBlock.clear();
    bitsInCurrentBlock = 0;
    (*f)(state.get());
}

void Sponge::squeezeIntoCurrentBlock()
{
    // bitsInCurrentBlock is assumed to be zero
    for(unsigned int i=0; i<rate/8; ++i)
        currentBlock.push_back(state.get()[i]);
    bitsInCurrentBlock = rate;
    (*f)(state.get());
}

void Sponge::flushAndSwitchToSqueezingPhase()
{
    if (bitsInCurrentBlock == rate)
        absorbCurrentBlock();
    if (bitsInCurrentBlock > 0)
        throw SpongeException("The input must contain a whole number of blocks.");
    if (!validInput)
        throw SpongeException("The input is not a valid sponge input.");
    squeezing = true;
}

string Sponge::getDescription() const
{
    stringstream a;
    a << "Sponge(r=" << dec << rate << ", c=" << capacity << ") using " << (*f);
    return a.str();
}

ostream& operator<<(ostream& a, const Sponge& sponge)
{
    return a << sponge.getDescription();
}
