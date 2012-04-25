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

#ifndef _SPONGETREE_H_
#define _SPONGETREE_H_

#include "sponge.h"

using namespace std;

/** This abstract class implements the creation of sponge functions.
  * This is for use with ShortLeafInterleavedSpongeTree.
  */
class SpongeFactory {
public:
    /** This abstract method allocates a new Sponge object.
      */
    virtual Sponge* newSponge() const = 0;
};

/** This class implements a tree hashing mode.
  * In particular, it implements the second tree hashing mode described 
  * in Section 7.2 of [Bertoni et al., Sufficient conditions for sound tree 
  * and sequential hashing modes, IACR ePrint 2009/210], 
  * where the tree has fixed height <i>H</i>=2, 
  * the final node has (parameterized) degree <i>D</i>,
  * and the message bits are interleaved by blocks of (parameterized) 
  * <i>B</i> bits onto the <i>D</i> leaves.
  * When the leaves are done absorbing the message blocks, a bit 0 is appended
  * for domain separation with the final node.
  * The chaining value size <i>C</i> is equal to the capacity of the
  * underlying sponge function, rounded up to the next multiple of 8 bits.
  * The final node contains:
  * - the concatenation of the <i>D</i> chaining values each of size <i>C</i>;
  * - the value <i>B</i> coded on 32-bit in a little-endian fashion;
  * - 5 bits set to zero (reserved for future use, for other layouts);
  * - the final bit set to 1 for domain separation with leaves.
  */
class ShortLeafInterleavedSpongeTree {
protected:
    /** This attributes allows the class to create new sponge functions
      * as needed.
      */
    const SpongeFactory& factory;
    /** This attribute contains the message blocks not yet processed.
      */
    MessageQueue absorbQueue;
    /** Boolean indicating whether the sponge is in the squeezing phase 
      * (true) or in the absorbing phase (false). */
    bool squeezing;
    /** This attribute contains the index of the leaf being filled.
      */
    int leafIndex;
    /** This attribute contains the sponge functions for the <i>D</i> leaves.
      */
    vector<Sponge*> leaves;
    /** This attribute contains the sponge function for the final node.
      */
    Sponge* final;
    /** The degree of the final node.
      */
    int D;
    /** The block size.
      */
    int B;
    /** The chaining value size.
      */
    int C;
public:
    /** The constructor.
      * @param  aFactory    The object that will create the right sponge functions 
      *                     for the leaves and for the final node.
      * @param  aD  The value for <i>D</i>.
      * @param  aB  The value for <i>B</i>.
      */
    ShortLeafInterleavedSpongeTree(const SpongeFactory& aFactory, const int aD, const int aB);
    /** The desctructor. It deallocates any allocated object.
      */
    ~ShortLeafInterleavedSpongeTree();
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
protected:
    /**
      * Internal method that absorbs the data still in absorbQueue, 
      * and then switches the sponge function of the final node 
      * to the squeezing phase.
      */
    void flushAndSwitchToSqueezingPhase();
};

#endif
