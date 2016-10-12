/*
Implementation by the Keccak, Keyak and Ketje Teams, namely, Guido Bertoni,
Joan Daemen, Michaël Peeters, Gilles Van Assche and Ronny Van Keer, hereby
denoted as "the implementer".

For more information, feedback or questions, please refer to our websites:
http://keccak.noekeon.org/
http://keyak.noekeon.org/
http://ketje.noekeon.org/

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/
*/

#ifndef _KEYAKV2_H_
#define _KEYAKV2_H_

#include <sstream>
#include "Keccak-f.h"
#include "Motorist.h"
#include "types.h"

using namespace std;

class Keyak {
protected:
    KeccakP f;
    unsigned int W;
    unsigned int Pi;
    unsigned int c;
    unsigned int tau;
    Motorist motorist;
public:
    Keyak(unsigned int b, unsigned int nr, unsigned int Pi, unsigned int c, unsigned int tau);
    Keyak(const Keyak& keyak);
    bool StartEngine(const string& K, const string& N, bool tagFlag, stringstream& T, bool unwrapFlag, bool forgetFlag);
    bool Wrap(istream& I, stringstream& O, istream& A, stringstream& T, bool unwrapFlag, bool forgetFlag);
    friend ostream& operator<<(ostream& a, const Keyak& piston);
    unsigned int getWidth() const;
    unsigned int getPi() const;
};

class RiverKeyak : public Keyak {
public:
    RiverKeyak() : Keyak(800, 12, 1, 256, 128) {}
};

class LakeKeyak : public Keyak {
public:
    LakeKeyak() : Keyak(1600, 12, 1, 256, 128) {}
};

class SeaKeyak : public Keyak {
public:
    SeaKeyak() : Keyak(1600, 12, 2, 256, 128) {}
};

class OceanKeyak : public Keyak {
public:
    OceanKeyak() : Keyak(1600, 12, 4, 256, 128) {}
};

class LunarKeyak : public Keyak {
public:
    LunarKeyak() : Keyak(1600, 12, 8, 256, 128) {}
};

#endif
