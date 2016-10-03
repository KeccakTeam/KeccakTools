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

#ifndef _KETJEV2_H_
#define _KETJEV2_H_

#include <string>
#include "monkey.h"
#include "Keccak-f.h"
#include "types.h"

using namespace std;

class Ketje {
protected:
    IterableTransformation<KeccakPStar>  f;
    MonkeyWrap                           monkeyWrap;
public:
    Ketje(unsigned int width, unsigned int rho);
    void          initialize(const string &K, const string &N);
    string        wrap(const string &A, const string &B, unsigned int ell, string &T);
    string        unwrap(const string &A, const string &C, const string &T);
    unsigned int  getWidth() const;
};

class KetjeJr: public Ketje {
public:
    KetjeJr() : Ketje(200, 16) {}
};

class KetjeSr: public Ketje {
public:
    KetjeSr() : Ketje(400, 32) {}
};

class KetjeMinor: public Ketje {
public:
    KetjeMinor() : Ketje(800, 128) {}
};

class KetjeMajor: public Ketje {
public:
    KetjeMajor() : Ketje(1600, 256) {}
};

#endif
