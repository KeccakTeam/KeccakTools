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

#ifndef _KECCAK_H_
#define _KECCAK_H_

#include "Keccak-f.h"
#include "sponge.h"

/**
  * Class that implements the Keccak sponge function family.
  */
class Keccak : public Sponge {
public:
    /**
      * The constructor. It dynamically allocates a KeccakF permutation.
      *
      * @param  aRate       The desired rate (in bits) of the Keccak sponge 
      *                     function.
      * @param  aCapacity   The desired capacity (in bits) of the Keccak sponge 
      *                     function. The sum of the rate and capacity must
      *                     be equal to the width of one of the Keccak-<i>f</i> 
      *                     permutations.
      */
    Keccak(unsigned int aRate = 1024, unsigned int aCapacity = 576);
    /**
      * The destructor. It frees the allocated KeccakF permutation.
      */
    virtual ~Keccak();
    /**
      * Method that returns a string describing the instance of the Keccak sponge 
      * function.
      */
    string getDescription() const;
};

/**
  * Class that implements the (potentially) reduced-round Keccak sponge function family.
  */
class ReducedRoundKeccak : public Sponge {
protected:
    unsigned int nrRounds;
public:
    /**
      * The constructor. It dynamically allocates a KeccakF permutation.
      *
      * @param  aRate       The desired rate (in bits) of the Keccak sponge 
      *                     function.
      * @param  aCapacity   The desired capacity (in bits) of the Keccak sponge 
      *                     function. The sum of the rate and capacity must
      *                     be equal to the width of one of the Keccak-<i>f</i> 
      *                     permutations.
      * @param  aNrRounds     The desired number of rounds used in Keccak-<i>f</i>.
      */
    ReducedRoundKeccak(unsigned int aRate, unsigned int aCapacity, unsigned int aNrRounds);
    /**
      * The destructor. It frees the allocated KeccakF permutation.
      */
    virtual ~ReducedRoundKeccak();
    /**
      * Method that returns a string describing the instance of the Keccak sponge 
      * function.
      */
    string getDescription() const;
};

#endif
