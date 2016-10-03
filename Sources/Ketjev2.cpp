/*
 * Implementation by the Keccak, Keyak and Ketje Teams, namely, Guido Bertoni,
 * Joan Daemen, Michaël Peeters, Gilles Van Assche and Ronny Van Keer, hereby
 * denoted as "the implementer".
 *
 * For more information, feedback or questions, please refer to our websites:
 * http://keccak.noekeon.org/
 * http://keyak.noekeon.org/
 * http://ketje.noekeon.org/
 *
 * To the extent possible under law, the implementer has waived all copyright
 * and related or neighboring rights to the source code in this file.
 * http://creativecommons.org/publicdomain/zero/1.0/
 */

#include "Ketjev2.h"

Ketje::Ketje(unsigned int width, unsigned int rho)
    : f(width), monkeyWrap(f, rho, 12, 1, 6)
{}

void Ketje::initialize(const string &K, const string &N)
{
    monkeyWrap.initialize(BitString(K), BitString(N));
}

string Ketje::wrap(const string &A, const string &B, unsigned int ell, string &T)
{
    if ( ell % 8 != 0 ) {
        throw Exception("This implementation restricts ell to multiple of 8."); // Actually a limitation of the interface (string class)
    }
    BitString  Tbits(T);                                                        // Tbits takes T by reference, and will automatically update it
    return monkeyWrap.wrap(BitString(A), BitString(B), ell, Tbits).str();
}

string Ketje::unwrap(const string &A, const string &C, const string &T)
{
    return monkeyWrap.unwrap(BitString(A), BitString(C), BitString(T)).str();
}

unsigned int Ketje::getWidth() const
{
    return f.width;
}
