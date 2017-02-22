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

#ifndef _SPONGE_H_
#define _SPONGE_H_

#include <deque>
#include <iostream>
#include <memory>
#include "padding.h"
#include "transformations.h"
#include "types.h"

using namespace std;

typedef Exception SpongeException;

/**
  * Class implementing the sponge construction.
  * This class uses a given transformation (inherited from class
  * Transformation or Permutation).
  * The following restrictions are made:
  * - the input message can have any length (not necessarily a multiple
  *     of 8), but all calls to absorb() <em>except the last one</em>
  *     must give a multiple of 8 bits;
  * - if the rate is a multiple of 8, the output stream can be given per byte only;
  * - if the rate is not a multiple of 8, the output stream can be given per whole blocks only.
  */
class Sponge {
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
    /** The capacity of the sponge construction. */
    unsigned int capacity;
    /** The rate of the sponge construction. */
    unsigned int rate;
    /** Boolean indicating whether the sponge is in the squeezing phase
      * (true) or in the absorbing phase (false). */
    bool squeezing;
    /** The state of the sponge construction. */
    auto_ptr<UINT8> state;
    /** The message blocks not yet absorbed. */
    MessageQueue absorbQueue;
    /** Buffer containing the partial block that is being squeezed. */
    deque<UINT8> squeezeBuffer;
public:
    /**
      * The constructor. The transformation, padding rule and rate are given to the
      * constructor, while the capacity is computed from the function width
      * and the requested rate.
      * The sponge construction is set to the absorbing phase.
      *
      * @param  aF          A pointer to the transformation used in the
      *                     sponge construction.
      * @param  aPad        A pointer to the padding rule used in the
      *                     sponge construction.
      * @param  aRate      The desired value of the rate (in bits),
      *                    not necessarily a multiple of 8.
      */
    Sponge(const Transformation *aF, const PaddingRule *aPad, unsigned int aRate);
    /**
      * The copy constructor.
      * @param  other   A reference to the object to copy from.
      */
    Sponge(const Sponge& other);
    /**
      * The (virtual) destructor.
      */
    virtual ~Sponge(void) {}
    /** Method to reset the sponge to the initial state. */
    void reset();
    /**
      * Method that absorbs data. The data is given as a sequence of bytes.
      * Within each byte, the bits are understood to be ordered from the
      * least significant bit to the most significant bit.
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
      * Method that absorbs data. The data is given as a sequence of bytes.
      * Within each byte, the bits are understood to be ordered from the
      * least significant bit to the most significant bit.
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
    void absorb(const vector<UINT8>& input, unsigned int lengthInBits);
    /**
      * Method to extract data from the squeezing phase. If in the
      * absorbing phase, this function also switches to the squeezing phase.
      *
      * @param  output      The buffer where to store the squeezed data.
      * @param  desiredLengthInBits     The length in bits of the output.
      *                     If the rate of the sponge is a multiple of 8,
      *                     @a desiredOutputLength must be a multiple of 8.
      *                     Otherwise, @a desiredOutputLength must be equal to the rate.
      */
    void squeeze(UINT8 *output, unsigned int desiredLengthInBits);
    /**
      * Method to extract data from the squeezing phase. If in the
      * absorbing phase, this function also switches to the squeezing phase.
      *
      * @param  output      The buffer where to store the squeezed data.
      * @param  desiredLengthInBits     The length in bits of the output.
      *                     If the rate of the sponge is a multiple of 8,
      *                     @a desiredOutputLength must be a multiple of 8.
      *                     Otherwise, @a desiredOutputLength must be equal to the rate.
      */
    void squeeze(vector<UINT8>& output, unsigned int desiredLengthInBits);
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
      * Internal method that does the actual absorbing of the whole block
      * in @a block.
      * @param  block  A block to absorb.
      */
    void absorbBlock(const vector<UINT8>& block);
    /**
      * Internal method that absorbs the data still in absorbQueue,
      * and then switches the sponge function to the squeezing phase.
      */
    void flushAndSwitchToSqueezingPhase();
    /**
      * Internal method that does the actual squeezing and stores the whole
      * squeezed block into squeezeBuffer.
      */
    void squeezeIntoBuffer();
    /**
      * Internal method that copies a block from state to squeezeBuffer.
      */
    void fromStateToSqueezeBuffer();
};

#endif
