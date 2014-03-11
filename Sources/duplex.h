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

#ifndef _DUPLEX_H_
#define _DUPLEX_H_

#include <iostream>
#include <memory>
#include "padding.h"
#include "transformations.h"
#include "types.h"

using namespace std;

typedef Exception DuplexException;

/**
  * Class implementing the duplex construction.
  */
class Duplex {
protected:
    /** The transformation (or permutation) used by the sponge construction.
      * The memory allocated by f is assumed to belong to the caller;
      * this class does not free the allocated memory.
      */
    const Transformation *f;
    /** The padding rule used by the sponge construction.
      * The memory allocated by pad is assumed to belong to the caller;
      * this class does not free the allocated memory.
      */
    const PaddingRule *pad;
    /** The capacity of the duplex construction. */
    unsigned int capacity;
    /** The rate of the duplex construction. */
    unsigned int rate;
    /** The maximum input length of the duplex construction. */
    unsigned int rho_max;
    /** The state of the duplex construction. */
    auto_ptr<UINT8> state;
public:
    /**
      * The constructor. The transformation, padding rule and rate are given to the
      * constructor, while the capacity is computed from the function width
      * and the requested rate.
      * The duplex construction is initialized.
      *
      * @param  aF          A pointer to the transformation used in the
      *                     sponge construction.
      * @param  aPad        A pointer to the padding rule used in the
      *                     sponge construction.
      * @param  aRate      The desired value of the rate (in bits),
      *                    not necessarily a multiple of 8.
      */
    Duplex(const Transformation *aF, const PaddingRule *aPad, unsigned int aRate);

    /**
      * The copy constructor.
      * @param  other   A reference to the object to copy from.
      */
    Duplex(const Duplex& other);

    /**
      * Method that performs a duplexing call.
      * The input data is given as a sequence of bytes.
      * Within each byte, the bits are understood to be ordered from the
      * least significant bit to the most significant bit.
      * The output data are given using the same structure.
      *
      * @param  input       The input data. When inputLengthInBits is not
      *                     a multiple of 8, the last bits of data must be
      *                     in the least significant bits of the last byte.
      * @param  inputLengthInBits   The length in bits of the data provided in input.
      *                             This value does not need to be a multiple of 8.
      * @param  output      The buffer where to store the output data.
      * @param  desiredOutputLengthInBits    The length in bits of the output.
      *
      * @pre    inputLengthInBits ≤ getMaximumInputLength()
      * @pre    desiredOutputLengthInBits ≤ getMaximumOutputLength()
      */
    void duplexing(const UINT8 *input, unsigned int inputLengthInBits, UINT8 *output, unsigned int desiredOutputLengthInBits);

    /**
      * Method that performs a duplexing call.
      * The input data is given as a sequence of bytes.
      * Within each byte, the bits are understood to be ordered from the
      * least significant bit to the most significant bit.
      * The output data are given using the same structure.
      *
      * @param  inputStart  Constant iterator for the input data.
      *                     When inputLengthInBits is not
      *                     a multiple of 8, the last bits of data must be
      *                     in the least significant bits of the last byte.
      * @param  inputLengthInBits   The length in bits of the data provided in input.
      *                             This value does not need to be a multiple of 8.
      * @param  output      The buffer where to store the output data.
      *                     The output is appended to the existing content of the vector.
      * @param  desiredOutputLengthInBits    The length in bits of the output.
      *
      * @pre    inputLengthInBits ≤ getMaximumInputLength()
      * @pre    desiredOutputLengthInBits ≤ getMaximumOutputLength()
      */
    template<class InputIterator, class OutputContainer>
    void duplexing(InputIterator inputStart, unsigned int inputLengthInBits, OutputContainer& output, unsigned int desiredOutputLengthInBits)
    {
        MessageQueue queue(rate);
        queue.append(inputStart, inputLengthInBits);
        const UINT8* state = processDuplexing(queue);
        outputDuplexing(state, output, desiredOutputLengthInBits);
    }

    /**
      * Method that performs a duplexing call.
      * @param  sigmaBegin  Pointer to the first part of the input σ given as bytes.
      *                     Trailing bits are given in @a delimitedSigmaEnd.
      * @param  sigmaBeginByteLen   The number of input bytes provided in @a sigmaBegin.
      * @param  delimitedSigmaEnd   Byte containing from 0 to 7 trailing bits that must be
      *                     appended to the input data in @a sigmaBegin.
      *                     These <i>n</i>=|σ| mod 8 bits must be in the least significant bit positions.
      *                     These bits must be delimited with a bit 1 at position <i>n</i>
      *                     (counting from 0=LSB to 7=MSB) and followed by bits 0
      *                     from position <i>n</i>+1 to position 7.
      *                     Some examples:
      *                         - If |σ| is a multiple of 8, then @a delimitedSigmaEnd must be 0x01.
      *                         - If |σ| mod 8 is 1 and the last bit is 1 then @a delimitedSigmaEnd must be 0x03.
      *                         - If |σ| mod 8 is 4 and the last 4 bits are 0,0,0,1 then @a delimitedSigmaEnd must be 0x18.
      *                         - If |σ| mod 8 is 6 and the last 6 bits are 1,1,1,0,0,1 then @a delimitedSigmaEnd must be 0x67.
      *                     .
      * @param  Z           Pointer to the buffer where to store the output data Z.
      * @param  ZByteLen    The number of output bytes desired for Z.
      * @note   The input bits σ are the result of the concatenation of the bytes in @a sigmaBegin
      *                     and the bits in @a delimitedSigmaEnd before the delimiter.
      * @pre    @a delimitedSigmaEnd ≠ 0x00
      * @pre    @a sigmaBeginByteLen*8+<i>n</i> ≤ getMaximumInputLength()
      * @pre    @a ZByteLen*8 ≤ getMaximumOutputLength()
      */
    void duplexingBytes(const UINT8 *sigmaBegin, unsigned int sigmaBeginByteLen, UINT8 delimitedSigmaEnd, UINT8 *Z, unsigned int ZByteLen);

    /**
      * Method that performs a duplexing call.
      * @param  sigmaBegin_start    Constant iterator marking the begin of the first part of the input σ given as bytes.
      * @param  sigmaBegin_stop     Constant iterator marking the end of the first part of the input σ given as bytes.
      *                     Trailing bits are given in @a delimitedSigmaEnd.
      * @param  delimitedSigmaEnd   Byte containing from 0 to 7 trailing bits that must be
      *                     appended to the input data in @a sigmaBegin.
      *                     These <i>n</i>=|σ| mod 8 bits must be in the least significant bit positions.
      *                     These bits must be delimited with a bit 1 at position <i>n</i>
      *                     (counting from 0=LSB to 7=MSB) and followed by bits 0
      *                     from position <i>n</i>+1 to position 7.
      * @param  Z           Buffer to store the output data Z.
      *                     The output is appended to the existing content.
      * @param  ZByteLen    The number of output bytes desired for Z.
      * @pre    @a delimitedSigmaEnd ≠ 0x00
      * @pre    @a sigmaBegin.size()*8+<i>n</i> ≤ getMaximumInputLength()
      * @pre    @a ZByteLen*8 ≤ getMaximumOutputLength()
      */
    template<class InputIterator, class OutputContainer>
    void duplexingBytes(InputIterator sigmaBegin_start, InputIterator sigmaBegin_stop, UINT8 delimitedSigmaEnd, OutputContainer& Z, unsigned int ZByteLen)
    {
        MessageQueue queue(rate);
        queue.append(sigmaBegin_start, sigmaBegin_stop);
        const UINT8* state = processDuplexing(queue, delimitedSigmaEnd);
        outputDuplexing(state, Z, ZByteLen*8);
    }

    /**
      * Method that performs a duplexing call without any requested output.
      * @param  sigmaBegin_start    Constant iterator marking the begin of the first part of the input σ given as bytes.
      * @param  sigmaBegin_stop     Constant iterator marking the end of the first part of the input σ given as bytes.
      *                     Trailing bits are given in @a delimitedSigmaEnd.
      * @param  delimitedSigmaEnd   Byte containing from 0 to 7 trailing bits that must be
      *                     appended to the input data in @a sigmaBegin.
      *                     These <i>n</i>=|σ| mod 8 bits must be in the least significant bit positions.
      *                     These bits must be delimited with a bit 1 at position <i>n</i>
      *                     (counting from 0=LSB to 7=MSB) and followed by bits 0
      *                     from position <i>n</i>+1 to position 7.
      * @pre    @a delimitedSigmaEnd ≠ 0x00
      * @pre    @a sigmaBegin.size()*8+<i>n</i> ≤ getMaximumInputLength()
      */
    template<class InputIterator>
    void duplexingBytes(InputIterator sigmaBegin_start, InputIterator sigmaBegin_stop, UINT8 delimitedSigmaEnd)
    {
        MessageQueue queue(rate);
        queue.append(sigmaBegin_start, sigmaBegin_stop);
        processDuplexing(queue, delimitedSigmaEnd);
    }

    /**
      * Method that returns the capacity of the sponge function.
      */
    unsigned int getCapacity() const;
    /**
      * Method that returns the maximum input length in bits of the duplexing() method.
      */
    unsigned int getMaximumInputLength() const;
    /**
      * Method that returns the maximum output length in bits of the duplexing() method.
      */
    unsigned int getMaximumOutputLength() const;
    /**
      * Method that returns a string with a description of itself.
      */
    virtual string getDescription() const;
    /**
      * Method that prints a brief description of the sponge function.
      */
    friend ostream& operator<<(ostream& a, const Duplex& duplex);
protected:
    /** Internal method to compute ρ_max, the maximum input length.
      * This value is such that an input message fits in one block after padding.
      */
    void computeRhoMax();
    const UINT8* processDuplexing(MessageQueue& queue, UINT8 delimitedSigmaEnd = 0x01);
    template<class OutputContainer>
    void outputDuplexing(const UINT8* state, OutputContainer& output, unsigned int desiredOutputLengthInBits)
    {
        if (desiredOutputLengthInBits > rate)
            throw DuplexException("The given output length must be less than or equal to the rate.");
        for(unsigned int i=0; i<desiredOutputLengthInBits/8; ++i)
            output.push_back(state[i]);
        if ((desiredOutputLengthInBits % 8) != 0) {
            UINT8 lastByte = state[desiredOutputLengthInBits/8];
            UINT8 mask = (1 << (desiredOutputLengthInBits % 8)) - 1;
            output.push_back(lastByte & mask);
        }
    }
};

#endif
