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
#include <string.h>
#include <vector>
#include "duplex.h"

using namespace std;

Duplex::Duplex(const Transformation *aF, const PaddingRule *aPad, unsigned int aRate)
    : f(aF), pad(aPad), rate(aRate)
{
    unsigned int width = f->getWidth();
    if (rate <= 0)
        throw DuplexException("The requested rate must be strictly positive.");
    if (rate > width)
        throw DuplexException("The requested rate is too large when using this function.");
    if (!pad->isRateValid(rate))
        throw DuplexException("The requested rate is incompatible with the padding rule.");
    capacity = width-rate;
    state.reset(new UINT8[(width+7)/8]);
    for(unsigned int i=0; i<(width+7)/8; i++)
        state.get()[i] = 0;
    computeRhoMax();
}

Duplex::Duplex(const Duplex& other)
    : f(other.f), pad(other.pad), capacity(other.capacity), rate(other.rate), rho_max(other.rho_max)
{
    unsigned int width = f->getWidth();
    state.reset(new UINT8[(width+7)/8]);
    memcpy(state.get(), other.state.get(), (width+7)/8);
}

void Duplex::computeRhoMax()
{
    rho_max = 0;
    unsigned int inputSize = 0;
    while(pad->getPaddedSize(rate, inputSize) <= rate) {
        rho_max = inputSize;
        inputSize++;
    }
}

unsigned int Duplex::getCapacity() const
{
    return capacity;
}

unsigned int Duplex::getMaximumInputLength() const
{
    return rho_max;
}

unsigned int Duplex::getMaximumOutputLength() const
{
    return rate;
}

void Duplex::duplexing(const UINT8 *input, unsigned int inputLengthInBits, UINT8 *output, unsigned int desiredOutputLengthInBits)
{
    vector<UINT8> outputAsVector;
    duplexing(input, inputLengthInBits, outputAsVector, desiredOutputLengthInBits);
    for(unsigned int i=0; i<outputAsVector.size(); i++)
        output[i] = outputAsVector[i];
}

const UINT8* Duplex::processDuplexing(MessageQueue& queue, UINT8 delimitedSigmaEnd)
{
    if (delimitedSigmaEnd == 0x00)
        throw DuplexException("delimitedSigmaEnd has an invalid coding.");
    while (delimitedSigmaEnd != 0x01) {
        queue.appendBit(delimitedSigmaEnd & 1);
        delimitedSigmaEnd >>= 1;
    }
    queue.pad(*pad);
    if ((queue.blockCount() != 1) || (!queue.firstBlockIsWhole()))
        throw DuplexException("The given input length must be such that it spans exactly one block after padding.");
    const vector<UINT8>& block = queue.firstBlock();
    for(vector<UINT8>::size_type i=0; i<block.size(); ++i)
        state.get()[i] ^= block[i];
    (*f)(state.get());
    return state.get();
}

void Duplex::duplexingBytes(const UINT8 *sigmaBegin, unsigned int sigmaBeginByteLen, UINT8 delimitedSigmaEnd, UINT8 *Z, unsigned int ZByteLen)
{
    vector<UINT8> ZAsVector;
    duplexingBytes(sigmaBegin, sigmaBegin+sigmaBeginByteLen, delimitedSigmaEnd, ZAsVector, ZByteLen);
    for(unsigned int i=0; i<ZAsVector.size(); i++)
        Z[i] = ZAsVector[i];
}

string Duplex::getDescription() const
{
    stringstream a;
    a << "Duplex[f=" << (*f) << ", pad=" << (*pad)
        << ", r=" << dec << rate
        << ", c=" << capacity
        << ", \xCF\x81max=" << rho_max
        << "]";
    return a.str();
}

ostream& operator<<(ostream& a, const Duplex& duplex)
{
    return a << duplex.getDescription();
}
