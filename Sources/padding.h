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

#ifndef _PADDING_H_
#define _PADDING_H_

#include <deque>
#include <iostream>
#include <vector>
#include "types.h"

using namespace std;

/** Class representing a message block whose size is not necessarily
  * a multiple of 8 bits.
  * The block can contain any number of bits.
  */
class MessageBlock {
private:
    /** The current block stored as a sequence of bytes.
      * If the number of bits is not a multiple of 8, the last byte contains the
      * last few bits in its least significant bits.
      */
    vector<UINT8> block;
    /** The number of bits in the block. */
    unsigned int bitsInBlock;
public:
    /** The constructor. */
    MessageBlock();
    /** Method to append one bit to the block.
      * @param  bitValue     The value (0 or 1) of the bit to append.
      */
    void appendBit(int bitValue);
    /** Method to append one byte to the block.
      * @param  byteValue     The value (0x00…0xFF) of the byte to append.
      */
    void appendByte(UINT8 byteValue);
    /** Method to append a series of bits with value '0'.
      * @param  count     The number of zeroes to append.
      */
    void appendZeroes(unsigned int count);
    /** Method that returns the number of bits in the block.
      * @return  The number of bits in the block.
      */
    unsigned int size() const;
    /** Method that returns a reference to the block.
      * @return  A constant reference to the block represented as a vector of bytes.
      * If the number of bits is not a multiple of 8, the last byte contains the
      * last few bits in its least significant bits.
      */
    const vector<UINT8>& get() const;
};

class PaddingRule;

/** Class representing a sequence of fixed-size blocks, except the last one,
  * whose size can be smaller.
  */
class MessageQueue {
private:
    /** The sequence of blocks of type MessageBlock . */
    deque<MessageBlock> queue;
    /** The size of the blocks in bits. */
    unsigned int blockSize;
public:
    /** The constructor.
      * @param  aBlockSize    The desired block size in bits.
      */
    MessageQueue(unsigned int aBlockSize);
    /** Method to append one bit to the sequence.
      * @param  bitValue     The value (0 or 1) of the bit to append.
      */
    void appendBit(int bitValue);
    /** Method to append one byte to the sequence.
      * @param  byteValue     The value (0x00…0xFF) of the byte to append.
      */
    void appendByte(UINT8 byteValue);
    /** Method to append a series of bits with value '0'.
      * @param  count     The number of zeroes to append.
      */
    void appendZeroes(unsigned int count);
    /** Method to append a number of bits to the sequence.
      * @param  inputStart  Constant iterator for the sequence of bytes to append.
      * If the number of bits is not a multiple of 8, the last byte contains the
      * last few bits in its least significant bits.
      * @param  lengthInBits The number of bits to append.
      */
    template<class InputIterator>
    void append(InputIterator inputStart, unsigned int lengthInBits)
    {
        if (lengthInBits == 0)
            return;
        InputIterator i = inputStart;
        while(lengthInBits >= 8) {
            appendByte(*i);
            ++i;
            lengthInBits -= 8;
        }
        if (lengthInBits > 0) {
            UINT8 lastByte = *i;
            for(unsigned int i=0; i<lengthInBits; i++) {
                appendBit(lastByte & 1);
                lastByte = lastByte >> 1;
            }
        }
    }
    /** Method to append a number of bytes to the sequence.
      * @param  inputStart  Constant iterator for the begin of the sequence of bytes to append.
      * @param  inputStop   Constant iterator for the end of the sequence of bytes to append.
      */
    template<class InputIterator>
    void append(InputIterator inputStart, InputIterator inputStop)
    {
        for(InputIterator i = inputStart; i != inputStop; ++i)
            appendByte(*i);
    }
    /** Method that performs the padding up to the block length.
      * @param  pad The padding rule.
      */
    void pad(const PaddingRule &pad);
    /** Method that returns the size of the last block.
      * @return  The size of the last block in bits.
      */
    unsigned int lastBlockSize() const;
    /** Method that returns the number of blocks in the sequence.
      * @return  The number of blocks in the sequence.
      */
    unsigned int blockCount() const;
    /** Method that tells whether the first block has exactly blockSize bits.
      * @return  True iff the first block has exactly blockSize bits.
      */
    bool firstBlockIsWhole() const;
    /** Method that returns a reference to the first block of the sequence.
      * @return  A constant reference to the block represented as a vector of bytes.
      * If the number of bits is not a multiple of 8, the last byte contains the
      * last few bits in its least significant bits.
      */
    const vector<UINT8>& firstBlock() const;
    /** Method that removes the first block of the sequence.
      */
    void removeFirstBlock();
    /** Method to empty the sequence. */
    void clear();
private:
    /** Private method to create a new block if the last one has reached blockSize bits. */
    void adjustLastBlock();
};

/**
  * Abstract class implementing a padding rule.
  */
class PaddingRule {
public:
    /**
      * Virtual destructor - necessary because this is an abstract class.
      */
    virtual ~PaddingRule() {};
    /** Abstract method to apply the padding on the given message queue.
      * @param  rate The block size in bits to which the padding must align.
      * @param  queue    The message bits to which the padding must be appended.
      */
    virtual void pad(unsigned int rate, MessageQueue& queue) const = 0;
    /** Abstract method to compute the size after padding.
      * @param  rate The block size in bits to which the padding must align.
      * @param  inputSize   The size in bits of the input message before padding.
      */
    virtual unsigned int getPaddedSize(unsigned int rate, unsigned int inputSize) const = 0;
    /** Abstract method to compute the minimum rate of a duplex object with the
      * given padding, given the maximum duplex rate rho_max.
      * @param  rate The block size in bits to which the padding must align.
      * @param  inputSize   The size in bits of the input message before padding.
      */
    unsigned int getDuplexRate(unsigned int rho_max) const;
    /**
      * Abstract method that returns a string with a description of itself.
      */
    virtual string getDescription() const = 0;
    /** Method to determine whether a given block size (or rate) is valid for the padding.
      * @param  rate The block size in bits.
      * @return  True iff the rate is valid.
      */
    virtual bool isRateValid(unsigned int rate) const;
    /**
      * Method that prints a brief description of the padding rule.
      */
    friend ostream& operator<<(ostream& a, const PaddingRule& pad);
protected:
    /** Method that appends a bit '1' then the minimum number of bits '0'
      * to get a whole number of blocks. After calling this method,
      * @a queue should contain a multiple of @a blockSize bits.
      * @param  blockSize  The block size in bits to which to align.
      */
    void append10star(unsigned int blockSize, MessageQueue& queue) const;
};

/** Class that implements the simple padding rule.
  * A bit '1' is appended, then the minimum number of bits '0'
  * to get a whole number of blocks.
  */
class SimplePadding : public PaddingRule {
public:
    /** The constructor. */
    SimplePadding() {}
    /** Actual method for SimplePadding, see PaddingRule::pad(). */
    void pad(unsigned int rate, MessageQueue& queue) const;
    /** Actual method for SimplePadding, see PaddingRule::getPaddedSize(). */
    unsigned int getPaddedSize(unsigned int rate, unsigned int inputSize) const;
    /** Actual method for SimplePadding, see PaddingRule::getDescription(). */
    string getDescription() const;
};

/** Class that implements the simple padding rule.
  * A bit '1' is appended, then the minimum number of bits '0'
  * to get a whole number of blocks after appending a final bit '1'.
  */
class MultiRatePadding : public PaddingRule {
public:
    /** The constructor. */
    MultiRatePadding() {}
    /** Actual method for MultiRatePadding, see PaddingRule::pad(). */
    void pad(unsigned int rate, MessageQueue& queue) const;
    /** Actual method for MultiRatePadding, see PaddingRule::getPaddedSize(). */
    unsigned int getPaddedSize(unsigned int rate, unsigned int inputSize) const;
    /** Actual method for MultiRatePadding, see PaddingRule::getDescription(). */
    string getDescription() const;
};

/** Class that implements the padding rule used by Keccak versions 1 and 2.
  */
class OldDiversifiedKeccakPadding : public PaddingRule {
protected:
    /** The 8-bit diversifier. */
    unsigned char diversifier;
public:
    /** The constructor. */
    OldDiversifiedKeccakPadding(unsigned char aDiversifier)
        : diversifier(aDiversifier) {}
    /** Actual method for OldDiversifiedKeccakPadding, see PaddingRule::pad(). */
    void pad(unsigned int rate, MessageQueue& queue) const;
    /** Actual method for OldDiversifiedKeccakPadding, see PaddingRule::getPaddedSize(). */
    unsigned int getPaddedSize(unsigned int rate, unsigned int inputSize) const;
    /** Actual method for OldDiversifiedKeccakPadding, see PaddingRule::getDescription(). */
    string getDescription() const;
    /** Actual method for OldDiversifiedKeccakPadding, see PaddingRule::isRateValid().
      * The rate must be a multiple of 8.
      */
    bool isRateValid(unsigned int rate) const;
};

template<class InputIterator>
vector<UINT8> getKeyPack(InputIterator keyStart, unsigned int keyLengthInBits, unsigned int packLengthInBits)
{
    if ((packLengthInBits%8) != 0)
        throw Exception("The pack length must be a multiple of 8 bits");
    if (packLengthInBits > 255*8)
        throw Exception("The pack cannot be longer than 255 bytes");

    MessageQueue keyPackQueue(packLengthInBits);
    keyPackQueue.appendByte((UINT8)(packLengthInBits/8));
    keyPackQueue.append(keyStart, keyLengthInBits);
    keyPackQueue.pad(SimplePadding());

    if (keyPackQueue.blockCount() != 1)
        throw Exception("The pack length is not large enough to make the key fit");

    return keyPackQueue.firstBlock();
}

#endif
