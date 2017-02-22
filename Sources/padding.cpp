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
#include "padding.h"

using namespace std;

MessageBlock::MessageBlock()
    : bitsInBlock(0)
{
}

void MessageBlock::appendBit(int bitValue)
{
    unsigned int bitsInByte = bitsInBlock % 8;
    if (bitsInByte == 0)
        block.push_back(bitValue & 1);
    else if (bitValue & 1)
        block.back() |= (UINT8)(1 << bitsInByte);
    bitsInBlock++;
}

void MessageBlock::appendByte(UINT8 byteValue)
{
    if ((bitsInBlock % 8) == 0) {
        block.push_back(byteValue);
        bitsInBlock += 8;
    }
    else {
        for(unsigned int i=0; i<8; i++) {
            appendBit(byteValue % 2);
            byteValue = byteValue / 2;
        }
    }
}

void MessageBlock::appendZeroes(unsigned int count)
{
    while((count > 0) && ((bitsInBlock % 8) != 0)) {
        appendBit(0);
        count--;
    }
    while(count >= 8) {
        block.push_back(0);
        bitsInBlock += 8;
        count -= 8;
    }
    while(count > 0) {
        appendBit(0);
        count--;
    }
}

unsigned int MessageBlock::size() const
{
    return bitsInBlock;
}

const vector<UINT8>& MessageBlock::get() const
{
    return block;
}

MessageQueue::MessageQueue(unsigned int aBlockSize)
    : blockSize(aBlockSize)
{
    queue.push_back(MessageBlock());
}

unsigned int MessageQueue::lastBlockSize() const
{
    return queue.back().size();
}

unsigned int MessageQueue::blockCount() const
{
    if (queue.size() == 0)
        return 0;
    else if (queue.back().size() == 0)
        return queue.size()-1;
    else
        return queue.size();
}

void MessageQueue::adjustLastBlock()
{
    if (lastBlockSize() >= blockSize)
        queue.push_back(MessageBlock());
}

void MessageQueue::appendBit(int bitValue)
{
    adjustLastBlock();
    queue.back().appendBit(bitValue);
}

void MessageQueue::appendByte(UINT8 byteValue)
{
    adjustLastBlock();
    if (lastBlockSize()+8 <= blockSize)
        queue.back().appendByte(byteValue);
    else {
        for(unsigned int i=0; i<8; i++) {
            appendBit(byteValue % 2);
            byteValue = byteValue / 2;
        }
    }
}

void MessageQueue::appendZeroes(unsigned int count)
{
    while(count > 0) {
        adjustLastBlock();
        unsigned int countInBlock = count;
        if (countInBlock > (blockSize - lastBlockSize()))
            countInBlock = blockSize - lastBlockSize();
        queue.back().appendZeroes(countInBlock);
        count -= countInBlock;
    }
}

void MessageQueue::pad(const PaddingRule &pad)
{
    pad.pad(blockSize, *this);
}

bool MessageQueue::firstBlockIsWhole() const
{
    return queue.front().size() == blockSize;
}

const vector<UINT8>& MessageQueue::firstBlock() const
{
    return queue.front().get();
}

void MessageQueue::removeFirstBlock()
{
    queue.pop_front();
    if (queue.size() == 0)
        queue.push_back(MessageBlock());
}

void MessageQueue::clear()
{
    queue.clear();
    queue.push_back(MessageBlock());
}


bool PaddingRule::isRateValid(unsigned int rate) const
{
    (void)rate;
    return true;
}

void PaddingRule::append10star(unsigned int blockSize, MessageQueue& queue) const
{
    queue.appendBit(1);
    if ((queue.lastBlockSize() % blockSize) != 0)
        queue.appendZeroes(blockSize - (queue.lastBlockSize() % blockSize));
}

unsigned int PaddingRule::getDuplexRate(unsigned int rho_max) const
{
    unsigned int rate = 0;
    bool rateIsSufficient = false;
    while ( rateIsSufficient == false ) {
        rate++;
        rateIsSufficient = true;
        unsigned int inputSize = 0;
        while ( ( rateIsSufficient == true ) && (inputSize <= rho_max) ) {
            rateIsSufficient = rateIsSufficient && (rate == getPaddedSize(rate,inputSize));
            inputSize++;
        }
    }
    return rate;
}


ostream& operator<<(ostream& a, const PaddingRule& pad)
{
    return a << pad.getDescription();
}

void SimplePadding::pad(unsigned int rate, MessageQueue& queue) const
{
    append10star(rate, queue);
}

unsigned int SimplePadding::getPaddedSize(unsigned int rate, unsigned int inputSize) const
{
    inputSize++;
    if ((inputSize % rate) != 0)
        inputSize += rate - (inputSize % rate);
    return inputSize;
}

string SimplePadding::getDescription() const
{
    return "pad10*";
}

void MultiRatePadding::pad(unsigned int rate, MessageQueue& queue) const
{
    queue.appendBit(1);
    queue.appendZeroes(rate - 1 - (queue.lastBlockSize() % rate));
    queue.appendBit(1);
}

unsigned int MultiRatePadding::getPaddedSize(unsigned int rate, unsigned int inputSize) const
{
    inputSize++;
    inputSize += rate - 1 - (inputSize % rate);
    inputSize++;
    return inputSize;
}

string MultiRatePadding::getDescription() const
{
    return "pad10*1";
}

void OldDiversifiedKeccakPadding::pad(unsigned int rate, MessageQueue& queue) const
{
    append10star(8, queue);
    queue.appendByte(diversifier);
    queue.appendByte(rate/8);
    append10star(rate, queue);
}

unsigned int OldDiversifiedKeccakPadding::getPaddedSize(unsigned int rate, unsigned int inputSize) const
{
    inputSize++;
    if ((inputSize % 8) != 0)
        inputSize += 8 - (inputSize % 8);
    inputSize += 8; // diversifier
    inputSize += 8; // rate
    inputSize++;
    if ((inputSize % rate) != 0)
        inputSize += rate - (inputSize % rate);
    return inputSize;
}

string OldDiversifiedKeccakPadding::getDescription() const
{
    stringstream str;
    str << "oldKeccakPadding[d=" << dec << (int)diversifier << "]";
    return str.str();
}

bool OldDiversifiedKeccakPadding::isRateValid(unsigned int rate) const
{
    return ((rate % 8) == 0);
}
