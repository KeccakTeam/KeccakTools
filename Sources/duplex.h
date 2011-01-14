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

/**
 * Exception that can be thrown by the class Duplex.
 */
class DuplexException {
public:
    /** A string expressing the reason for the exception. */
    string reason;
    /**
     * The constructor.
     * @param   aReason     A string giving the reason of the exception.
     */
    DuplexException(const string& aReason);
};

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
      * Method that returns the capacity of the sponge function.
      */
    unsigned int getCapacity();
    /**
      * Method that returns the maximum input length in bits of the duplexing() method.
      */
    unsigned int getMaximumInputLength();
    /**
      * Method that returns the maximum output length in bits of the duplexing() method.
      */
    unsigned int getMaximumOutputLength();
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
};

#endif
