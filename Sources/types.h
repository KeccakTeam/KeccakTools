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

#ifndef _TYPES_H_
#define _TYPES_H_

#include <string>

typedef unsigned char UINT8;
typedef unsigned short UINT16;
typedef unsigned int UINT32;
typedef unsigned long long UINT64;

/**
 * Exception with a string expressing the reason.
 */
class Exception {
public:
    /** A string expressing the reason for the exception. */
    std::string reason;

public:
    /**
     * The default constructor.
     */
    Exception()
        : reason() {}

    /**
     * The constructor.
     * @param   aReason     A string giving the reason of the exception.
     */
    Exception(const std::string& aReason)
        : reason(aReason) {}
};

#endif
