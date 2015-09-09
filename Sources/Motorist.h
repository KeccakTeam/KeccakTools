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

#ifndef _MOTORIST_H_
#define _MOTORIST_H_

#include <iostream>
#include <memory>
#include <sstream>
#include <vector>
#include "Keccak-f.h"
#include "types.h"

using namespace std;

class Piston {
protected:
    const Permutation *f;
    auto_ptr<UINT8> state;
    unsigned int Rs, Ra;
    unsigned int EOM, CryptEnd, InjectStart, InjectEnd;
public:
    Piston(const Permutation *aF, unsigned int aRs, unsigned int aRa);
    Piston(const Piston& other);
    void Crypt(istream& I, ostream& O, unsigned int omega, bool unwrapFlag);
    void Inject(istream& X, bool cryptingFlag);
    void Spark(bool eomFlag, unsigned int l);
    void GetTag(ostream& T, unsigned int l) const;
    friend ostream& operator<<(ostream& a, const Piston& piston);
};

class Engine {
protected:
    unsigned int Pi;
    vector<Piston>& Pistons;
    enum { fresh, crypted, endOfCrypt, endOfMessage } phase;
    vector<unsigned int> Et;
public:
    Engine(vector<Piston>& aPistons);
protected:
    void Spark(bool eomFlag, const vector<unsigned int>& l);
public:
    void Crypt(istream& I, ostream& O, bool unwrapFlag);
    void Inject(istream& A);
    void GetTags(ostream& T, const vector<unsigned int>& l);
    void InjectCollective(istream& X, bool diversifyFlag);
    friend ostream& operator<<(ostream& a, const Engine& engine);
};

class Motorist {
protected:
    unsigned int Pi;
    vector<Piston> Pistons;
    Engine engine;
    unsigned int W;
    unsigned int c;
    unsigned int cprime;
    unsigned int tau;
    enum { ready, riding, failed } phase;
public:
    Motorist(const Permutation *aF, unsigned int aPi, unsigned int aW, unsigned int ac, unsigned int atau);
    bool StartEngine(istream& SUV, bool tagFlag, stringstream& T, bool unwrapFlag, bool forgetFlag);
    bool Wrap(istream& I, stringstream& O, istream& A, stringstream& T, bool unwrapFlag, bool forgetFlag);
protected:
    void MakeKnot();
    bool HandleTag(bool tagFlag, stringstream& T, bool unwrapFlag);
public:
    friend ostream& operator<<(ostream& a, const Motorist& motorist);
};

#endif
