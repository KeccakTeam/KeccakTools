/*
Tools for the Keccak sponge function family.
Authors: Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche

This code is hereby put in the public domain. It is given as is, 
without any guarantee.

For more information, feedback or questions, please refer to our website:
http://keccak.noekeon.org/
*/

#ifndef _KECCAK_H_
#define _KECCAK_H_

#include "Keccak-f.h"
#include "sponge.h"

/**
  * Class that implements the Keccak sponge function family.
  */
class Keccak : public Sponge {
protected:
    /** The 8-bit diversifier. */
    unsigned char diversifier;
public:
    /**
      * The constructor. It dynamically allocates a KeccakF permutation.
      *
      * @param  aRate       The desired rate (in bits) of the Keccak sponge 
      *                     function. This must be a multiple of 8.
      * @param  aCapacity   The desired capacity (in bits) of the Keccak sponge 
      *                     function. The sum of the rate and capacity must
      *                     be equal to the width of one of the Keccak-f 
      *                     permutations.
      * @param  aDiversifier    The desired diversifier. This can be any 
      *                     value between 0 and 255.
      */
    Keccak(unsigned int aRate = 1024, unsigned int aCapacity = 576, unsigned char aDiversifier = 0);
    /**
      * The destructor. It frees the allocated KeccakF permutation.
      */
    virtual ~Keccak();
    /**
      * Method that returns a string describing the instance of the Keccak sponge 
      * function.
      */
    string getDescription() const;
    /**
      * Method that completes the last incomplete (possibly empty) input block 
      * according by padding according to the Keccak specifications.
      *
      * @pre This function must be used in the absorbing phase only.
      */
    void pad();
};

#endif
