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

#include "Keyakv2.h"

UINT8 enc8(unsigned int x);

Keyak::Keyak(unsigned int b, unsigned int nr, unsigned int Pi, unsigned int c, unsigned int tau)
    : f(b, nr), W(max((int)b/25, 8)), Pi(Pi), c(c), tau(tau), motorist(&f, Pi, W, c, tau)
{
}

Keyak::Keyak(const Keyak& keyak)
    : f(keyak.f), W(keyak.W), Pi(keyak.Pi), c(keyak.c), tau(keyak.tau), motorist(&f, Pi, W, c, tau)
{
}

string keypack(const string& K, unsigned int l)
{
    if (K.size() + 2 > l)
        throw Exception("The key is too big and does not fit in the key pack.");
    string result = char(enc8(l)) + K + char(1);
    while(result.size() < l)
        result += char(0);
    return result;
}

bool Keyak::StartEngine(const string& K, const string& N, bool tagFlag, stringstream& T, bool unwrapFlag, bool forgetFlag)
{
    unsigned int lk = W/8*((c+9+W-1)/W);
    stringstream SUV(keypack(K, lk) + N);
    return motorist.StartEngine(SUV, tagFlag, T, unwrapFlag, forgetFlag);
}

bool Keyak::Wrap(istream& I, stringstream& O, istream& A, stringstream& T, bool unwrapFlag, bool forgetFlag)
{
    return motorist.Wrap(I, O, A, T, unwrapFlag, forgetFlag);
}

ostream& operator<<(ostream& a, const Keyak& keyak)
{
    return a << "Keyak[b=" << dec << keyak.f.getWidth()
        << ", nr=" << keyak.f.getNumberOfRounds()
        << ", \316\240=" << keyak.Pi
        << ", c=" << keyak.c
        << ", \317\204=" << keyak.tau << "]";
}

unsigned int Keyak::getWidth() const
{
    return f.getWidth();
}

unsigned int Keyak::getPi() const
{
    return Pi;
}

