/*
Tools for the Keccak sponge function family.
Authors: Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche

This code is hereby put in the public domain. It is given as is, 
without any guarantee.

For more information, feedback or questions, please refer to our website:
http://keccak.noekeon.org/
*/

#ifndef _SPONGE_H_
#define _SPONGE_H_

#include <deque>
#include <iostream>
#include <memory>
#include "transformations.h"
#include "types.h"

using namespace std;

/**
 * Exception that can be thrown by the class Sponge.
 */
class SpongeException {
public:
    /** A string expressing the reason for the exception. */
    string reason;
    /**
     * The constructor.
     * @param   aReason     A string giving the reason of the exception.
     */
    SpongeException(const string& aReason);
};

/**
  * Class implementing the sponge construction.
  * This class uses a given transformation (inherited from class 
  * Transformation or Permutation). 
  * The following restrictions are made:
  * - the rate must be a multiple of 8 bits;
  * - the input message can have any length (not necessarily a multiple
  *     of 8), but all calls to absorb() <em>except the last one</em>
  *     must give a multiple of 8 bits;
  * - the output stream can be given per byte only.
  * 
  * The capacity and the function width can have any value, provided that the
  * rate is a multiple of 8 bits.
  */
class Sponge {
protected:
    /** The transformation (or permutation) used by the sponge construction. 
      * The memory allocated by f is assumed to belong to the caller;
      * this class does not free the allocated memory.
      */
    Transformation *f;
    /** The capacity of the sponge construction. */
    unsigned int capacity;
    /** The rate of the sponge construction. */
    unsigned int rate;
    /** Boolean indicating whether the sponge is in the squeezing phase 
      * (true) or in the absorbing phase (false). */
    bool squeezing;
    /** Boolean indicating whether a valid sponge input has been absorbed. */
    bool validInput;
    /** The state of the sponge construction. */
    auto_ptr<UINT8> state;
    /** Buffer containing the partial block before it is absorbed or the
      * partial block that is being squeezed. */
    deque<UINT8> currentBlock;
    /** Number of bits in currentBlock. */
    unsigned int bitsInCurrentBlock;
public:
    /**
      * The constructor. The transformation and capacity are given to the 
      * constructor, while the rate is computed from the function width
      * and the requested capacity. The rate must be a multiple of 8.
      * The sponge construction is set to the absorbing phase.
      *
      * @param  aF          A pointer to the transformation used in the 
      *                     sponge construction.
      * @param  aCapacity   The desired value of the capacity.
      */
    Sponge(Transformation *aF, unsigned int aCapacity);
    /**
      * Method that absorbs data. The data is given as a sequence of bytes.
      * Within each byte, the bits are understood to be ordered from the
      * least significant bit to the most significant bit.
      * The last incomplete block is kept in attribute currentBlock.
      *
      * @pre This function must be used in the absorbing phase only.
      *
      * @param  input       The data to absorb. When lengthInBits is not
      *                     a multiple of 8, the last bits of data must be
      *                     in the least significant bits of the last byte.
      * @param  lengthInBits    The length in bits of the data provided in
      *                     input. When lengthInBits is not a multiple
      *                     of 8, this function can no longer be used.
      */
    void absorb(const UINT8 *input, unsigned int lengthInBits);
    /**
      * Method to pad the data absorbed so far. By default, it calls
      * padCurrentBlock(rate), but it can be overridden to provide
      * other valid padding methods.
      *
      * @pre This function must be used in the absorbing phase only.
      */
    virtual void pad();
    /**
      * Method to extract data from the squeezing phase. If in the 
      * absorbing phase, this function also switches to the squeezing phase.
      *
      * @param  output      The buffer where to store the squeezed data.
      * @param  desiredLengthInBits    The length in bits of the output.
      *                     It must be a multiple of 8.
      */
    void squeeze(UINT8 *output, unsigned int desiredLengthInBits);
    /**
      * Method that returns the capacity of the sponge function.
      */
    unsigned int getCapacity();
    /**
      * Method that returns the rate of the sponge function.
      */
    unsigned int getRate();
    /**
      * Method that returns a string with a description of itself.
      */
    virtual string getDescription() const;
    /**
      * Method that prints a brief description of the sponge function.
      */
    friend ostream& operator<<(ostream& a, const Sponge& sponge);
protected:
    /**
      * Internal method to pad the data absorbed so far, by applying the simple
      * padding. It appends a single bit 1, then the smallest number of
      * bits 0 such that the data absorbed so far have a length that is a
      * multiple of @a desiredLengthFactor. 
      *
      * @param  desiredLengthFactor     The number of bits that divides
      *                     the length of @a currentBlock after padding.
      *                     It must be a strictly positive multiple of 8,
      *                     and it must divide @a rate.
      *
      * @pre This function must be used in the absorbing phase only.
      * @pre @a desiredLengthFactor must be a strictly positive multiple of 8,
      * and it must divide @a rate.
      */
    void padCurrentBlock(unsigned int desiredLengthFactor);
    /**
      * Internal method that does the actual absorbing of the whole block 
      * in @a currentBlock.
      * This method is called only when bitsInCurrentBlock is equal to the rate.
      */
    void absorbCurrentBlock();
    /**
      * Internal method that absorbs the data still in currentBlock, 
      * and then switches the sponge function to the squeezing phase.
      */
    void flushAndSwitchToSqueezingPhase();
    /**
      * Internal method that does the actual squeezing and stores the whole 
      * squeezed block into currentBlock.
      */
    void squeezeIntoCurrentBlock();
};

#endif
