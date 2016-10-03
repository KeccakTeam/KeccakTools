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
 *
 * and related or neighboring rights to the source code in this file.
 * http://creativecommons.org/publicdomain/zero/1.0/
 */

#include <string>
#include <vector>
#include "monkey.h"

using namespace std;

UINT8  enc8(unsigned int x);

static void _assert(bool condition, const string &synopsis, const char *fct)
{
    if ( !condition ) {
        throw Exception(string(fct) + "(): " + synopsis);
    }
}

#if defined(__GNUC__)
#define assert(cond, msg)  _assert(cond, msg, __PRETTY_FUNCTION__)
#else
#define assert(cond, msg)  _assert(cond, msg, __FUNCTION__)
#endif

MonkeyDuplex::MonkeyDuplex(BaseIterableTransformation &f,
                           unsigned int                r,
                           unsigned int                nStart,
                           unsigned int                nStep,
                           unsigned int                nStride)
    : f(f), r(r), nStart(nStart), nStep(nStep), nStride(nStride)
{
    assert((f.width % 8) == 0, "This implementation only supports permutation width that are multiple of 8."); // Limitation of Transformation class
    assert(2 < r,              "r must be greater than 2.");
    assert(r < f.width,        "r must be less than the permutation width.");
    assert(nStep < nStride,    "nStep must be less than nStride.");
}

void MonkeyDuplex::start(const BitString &I)
{
    assert(I.size() + 2 <= f.width, "I length must be less than or equal to the permutation width minus 2.");
    s = I || BitString::pad101(f.width, I.size());
    f[nStart](s.array());
}

BitString MonkeyDuplex::step(const BitString &sigma, unsigned int ell)
{
    assert(ell <= r,              "ell must be less than or equal to r.");
    assert(sigma.size() + 2 <= r, "sigma length must be less than or equal to r minus 2.");

    BitString  P;

    P = sigma || BitString::pad101(r, sigma.size());
    s = s ^ (P || BitString::zeroes(f.width - r));
    f[nStep](s.array());

    return BitString(s).truncate(ell);
}

BitString MonkeyDuplex::stride(const BitString &sigma, unsigned int ell)
{
    assert(ell <= r,              "ell must be less than or equal to r.");
    assert(sigma.size() + 2 <= r, "sigma length must be less than or equal to r minus 2.");

    BitString  P;

    P = sigma || BitString::pad101(r, sigma.size());
    s = s ^ (P || BitString::zeroes(f.width - r));
    f[nStride](s.array());

    return BitString(s).truncate(ell);
}

MonkeyWrap::MonkeyWrap(BaseIterableTransformation &f,
                       unsigned int                rho,
                       unsigned int                nStart,
                       unsigned int                nStep,
                       unsigned int                nStride)
    : f(f), rho(rho), D(f, rho + 4, nStart, nStep, nStride)
{
    assert(rho + 4 <= f.width, "rho must be lower than or equal to the permutation width minus 4.");
}

void MonkeyWrap::initialize(const BitString &K, const BitString &N)
{
    assert(K.size() + 18 <= f.width, "K length must be lower than or equal to the permutation width minus 18.");
    assert((K.size() % 8) == 0,      "K length must be a multiple of 8.");
    assert(
        N.size() + K.size() + 18 <= f.width,
        "N length must be lower than or equal to the permutation width minus K length and 18.");

    D.start(BitString::keypack(K, K.size() + 16) || N);
}

BitString MonkeyWrap::wrap(const BitString &Abits, const BitString &Bbits, unsigned int ell, BitString &T)
{
    BitString  Cbits, Z;
    Blocks     A(Abits, rho), B(Bbits, rho), C(Cbits, rho);

    for ( unsigned int i = 0; i + 2 <= A.size(); ++i ) {
        D.step(A[i] || 0 || 0, 0);
    }
    Z    = D.step(A[A.size() - 1] || 0 || 1, B[0].size());
    C[0] = B[0] ^ Z;
    for ( unsigned int i = 0; i + 2 <= B.size(); ++i ) {
        Z        = D.step(B[i] || 1 || 1, B[i + 1].size());
        C[i + 1] = B[i + 1] ^ Z;
    }
    T    = D.stride(B[B.size() - 1] || 1 || 0, rho);
    while ( T.size() < ell ) {
        T = T || D.step(BitString(0), rho);
    }
    T.truncate(ell);

    return Cbits;
}

BitString MonkeyWrap::unwrap(const BitString &Abits, const BitString &Cbits, const BitString &T)
{
    BitString  Bbits, Z, Tprime;
    Blocks     A(Abits, rho), B(Bbits, rho), C(Cbits, rho);

    for ( unsigned int i = 0; i + 2 <= A.size(); ++i ) {
        D.step(A[i] || 0 || 0, 0);
    }
    Z      = D.step(A[A.size() - 1] || 0 || 1, C[0].size());
    B[0]   = C[0] ^ Z;
    for ( unsigned int i = 0; i + 2 <= C.size(); ++i ) {
        Z        = D.step(B[i] || 1 || 1, C[i + 1].size());
        B[i + 1] = C[i + 1] ^ Z;
    }
    Tprime = D.stride(B[C.size() - 1] || 1 || 0, rho);
    while ( Tprime.size() < T.size()) {
        Tprime = Tprime || D.step(BitString(0), rho);
    }
    Tprime.truncate(T.size());
    if ( T == Tprime ) {
        return Bbits;
    }
    else {
        throw Exception("Tags do not match after unwrap.");
    }
}
