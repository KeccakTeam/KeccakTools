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
	unsigned int OmegaC, OmegaI;
public:
    Piston(const Permutation *f, unsigned int Rs, unsigned int Ra);
    Piston(const Piston& other);
    void Crypt(istream& I, ostream& O, bool decryptFlag);
    void Inject(istream& X);
    void Spark(void);
    void GetTag(ostream& T, unsigned int l);
    friend ostream& operator<<(ostream& a, const Piston& piston);
};

class Engine {
protected:
    vector<Piston>& Pistons;
public:
    Engine(vector<Piston>& Pistons);
public:
    void Wrap(istream& I, ostream& O, istream& A, bool decryptFlag);
    void GetTags(ostream& T, const vector<unsigned int>& l);
    void InjectCollective(istream& X, bool diversifyFlag);
    friend ostream& operator<<(ostream& a, const Engine& engine);
};

class Motorist {
protected:
    unsigned int Pi;
    unsigned int W;
    unsigned int c;
    unsigned int cprime;
    unsigned int tau;
    vector<Piston> Pistons;
    Engine engine;
    enum { ready, riding, failed } phase;
public:
    Motorist(const Permutation *f, unsigned int Pi, unsigned int W, unsigned int c, unsigned int tau);
    bool StartEngine(istream& SUV, bool tagFlag, stringstream& T, bool decryptFlag, bool forgetFlag);
    bool Wrap(istream& I, stringstream& O, istream& A, stringstream& T, bool decryptFlag, bool forgetFlag);
protected:
    void MakeKnot(void);
    bool HandleTag(bool tagFlag, stringstream& T, bool decryptFlag);
public:
    friend ostream& operator<<(ostream& a, const Motorist& motorist);
};

#endif
