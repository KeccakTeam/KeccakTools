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
#include <vector>
#include "duplex.h"

using namespace std;

DuplexException::DuplexException(const string& aReason)
    : reason(aReason)
{
}


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

void Duplex::computeRhoMax()
{
    rho_max = 0;
    unsigned int inputSize = 0;
    while(pad->getPaddedSize(rate, inputSize) <= rate) {
        rho_max = inputSize; 
        inputSize++;
    }
}

unsigned int Duplex::getCapacity()
{
    return capacity;
}

unsigned int Duplex::getMaximumInputLength()
{
    return rho_max;
}

unsigned int Duplex::getMaximumOutputLength()
{
    return rate;
}

void Duplex::duplexing(const UINT8 *input, unsigned int inputLengthInBits, UINT8 *output, unsigned int desiredOutputLengthInBits)
{
    MessageQueue queue(rate);
    queue.append(input, inputLengthInBits);
    pad->pad(rate, queue);
    if ((queue.blockCount() != 1) || (!queue.firstBlockIsWhole()))
        throw DuplexException("The given input length must be such that it spans exactly one block after padding.");
    const vector<UINT8>& block = queue.firstBlock();
    for(vector<UINT8>::size_type i=0; i<block.size(); ++i)
        state.get()[i] ^= block[i];
    (*f)(state.get());

    if (desiredOutputLengthInBits > rate)
        throw DuplexException("The given output length must be less than or equal to the rate.");
    for(unsigned int i=0; i<desiredOutputLengthInBits/8; ++i) {
        *output = state.get()[i];
        output++;
    }
    if ((desiredOutputLengthInBits % 8) != 0) {
        UINT8 lastByte = (state.get())[desiredOutputLengthInBits/8];
        UINT8 mask = (1 << (desiredOutputLengthInBits % 8)) - 1;
        *output = lastByte & mask;
    }
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
