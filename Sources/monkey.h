/*
 * KeccakTools
 *
 * The Keccak sponge function, designed by Guido Bertoni, Joan Daemen,
 * Michaël Peeters and Gilles Van Assche. For more information, feedback or
 * questions, please refer to our website: http://keccak.noekeon.org/
 *
 * Implementation by the designers,
 * hereby denoted as "the implementer".
 *
 * To the extent possible under law, the implementer has waived all copyright
 * and related or neighboring rights to the source code in this file.
 * http://creativecommons.org/publicdomain/zero/1.0/
 */

#ifndef _MONKEY_H_
#define _MONKEY_H_

#include <iostream>
#include <memory>
#include "padding.h"
#include "transformations.h"
#include "types.h"
#include "bitstring.h"

using namespace std;

/**
 * Class implementing an iterable permutation
 */
class BaseIterableTransformation {
public:
    const unsigned int            width;

    BaseIterableTransformation(unsigned int width) : width(width) {}

    virtual const Transformation &operator[](unsigned int n) = 0;
};

template<class T>
class IterableTransformation: public BaseIterableTransformation {
protected:
    T                     f;
public:
    IterableTransformation(unsigned int width) : BaseIterableTransformation(width), f(width, 1) {}

    const Transformation &operator[](unsigned int n)
    {
        f = T(width, n);
        return f;
    }
};

typedef Exception  MonkeyDuplexException;

/**
 * Class implementing the monkeyDuplex construction.
 */
class MonkeyDuplex {
protected:
    BaseIterableTransformation &f;
    const unsigned int          r;
    const unsigned int          nStart;
    const unsigned int          nStep;
    const unsigned int          nStride;
    BitString                   s;
public:
    MonkeyDuplex(BaseIterableTransformation &f, unsigned int r, unsigned int nStart, unsigned int nStep, unsigned int nStride);
    void       start(const BitString &I);
    BitString  step(const BitString &sigma, unsigned int ell);
    BitString  stride(const BitString &sigma, unsigned int ell);
};

typedef Exception  MonkeyWrapException;

/**
 * Class implementing the MonkeyWrap mode of use for authenticated encryption
 */
class MonkeyWrap {
protected:
    BaseIterableTransformation &f;
    const int                   rho;
    MonkeyDuplex                D;
public:
    MonkeyWrap(BaseIterableTransformation &f, unsigned int rho, unsigned int nStart, unsigned int nStep, unsigned int nStride);
    void       initialize(const BitString &K, const BitString &N);
    BitString  wrap(const BitString &A, const BitString &B, unsigned int ell, BitString &T);
    BitString  unwrap(const BitString &A, const BitString &C, const BitString &T);
};

#endif
