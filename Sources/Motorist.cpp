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

#include <string.h>
#include "Motorist.h"

bool hasMore(istream& in)
{
    return in.peek() >= 0;
}

UINT8 enc8(unsigned int x)
{
    if (x > 255) {
        stringstream s;
        s << "The integer " << x << " cannot be encoded on 8 bits.";
        throw Exception(s.str());
    }
    else
        return UINT8(x);
}

Piston::Piston(const Permutation *f, unsigned int Rs, unsigned int Ra)
    : f(f), Rs(Rs), Ra(Ra), OmegaC(0), OmegaI(0)
{
    unsigned int b = f->getWidth();
    if ((b%8) != 0)
        throw Exception("b is not a multiple of 8.");
    if (((b-32)/8) >= 248)
        throw Exception("(b-32)/8 is larger or equal to 248.");
    if (Rs > Ra)
        throw Exception("Rs is larger than Ra.");
    if (Ra > ((b-32)/8))
        throw Exception("Ra is larger than (b-32)/8.");
    state.reset(new UINT8[b/8]);
    memset(state.get(), 0, b/8);
    EOM = Ra;
    CryptEnd = Ra+1;
    InjectStart = Ra+2;
    InjectEnd = Ra+3;
}

Piston::Piston(const Piston& other)
    : f(other.f), Rs(other.Rs), Ra(other.Ra),
        EOM(other.EOM), CryptEnd(other.CryptEnd),
        InjectStart(other.InjectStart), InjectEnd(other.InjectEnd),
        OmegaC(0), OmegaI(0)
{
    unsigned int b = f->getWidth();
    state.reset(new UINT8[b/8]);
    memcpy(state.get(), other.state.get(), b/8);
}

void Piston::Crypt(istream& I, ostream& O, bool decryptFlag)
{
    while(hasMore(I) && (OmegaC < Rs)) {
        UINT8 x = I.get();
        O.put(state.get()[OmegaC] ^ x);
        if (decryptFlag)
            state.get()[OmegaC] = x;
        else
            state.get()[OmegaC] ^= x;
        OmegaC++;
    }
    state.get()[CryptEnd] ^= enc8(OmegaC);
    OmegaC = 0;
    OmegaI = Rs;
}

void Piston::Inject(istream& X)
{
    state.get()[InjectStart] ^= enc8(OmegaI);
    while(hasMore(X) && (OmegaI < Ra)) {
        state.get()[OmegaI] ^= X.get();
        OmegaI++;
    }
    state.get()[InjectEnd] ^= enc8(OmegaI);
    OmegaC = 0;
    OmegaI = 0;
}

void Piston::Spark(void)
{
    (*f)(state.get());
}

void Piston::GetTag(ostream& T, unsigned int l)
{
    if (l > Rs)
        throw Exception("The requested tag is too long.");
    if (l == 0)
        state.get()[EOM] ^= enc8(255);
    else
        state.get()[EOM] ^= enc8(l);
    Spark();
    for(unsigned int i=0; i<l; i++)
        T.put(state.get()[i]);
    OmegaC = l;
}

ostream& operator<<(ostream& a, const Piston& piston)
{
    return a << "Piston[f=" << (*piston.f) << ", Rs=" << dec << piston.Rs << ", Ra=" << piston.Ra << "]";
}

Engine::Engine(vector<Piston>& aPistons)
    : Pistons(aPistons)
{
}

void Engine::Wrap(istream& I, ostream& O, istream& A, bool decryptFlag)
{
    unsigned int Pi = Pistons.size();
    if(hasMore(I))
        for(unsigned int i=0; i<Pi; i++)
            Pistons[i].Crypt(I, O, decryptFlag);
    for(unsigned int i=0; i<Pi; i++)
        Pistons[i].Inject(A);
    if (hasMore(I) || hasMore(A))
        for(unsigned int i=0; i<Pi; i++)
            Pistons[i].Spark();
}

void Engine::GetTags(ostream& T, const vector<unsigned int>& l)
{
    unsigned int Pi = Pistons.size();
    for(unsigned int i=0; i<Pi; i++)
        Pistons[i].GetTag(T, l[i]);
}

void Engine::InjectCollective(istream& X, bool diversifyFlag)
{
    unsigned int Pi = Pistons.size();
    stringstream *Y = new stringstream[Pi];
    while(hasMore(X)) {
        UINT8 x = X.get();
        for(unsigned int i=0; i<Pi; i++)
            Y[i].put(x);
    }
    if (diversifyFlag) {
        for(unsigned int i=0; i<Pi; i++) {
            Y[i].put(enc8(Pi));
            Y[i].put(enc8(i));
        }
    }
    for(unsigned int i=0; i<Pi; i++)
        Y[i].seekg(0, ios_base::beg);
    while(hasMore(Y[0])) {
        for(unsigned int i=0; i<Pi; i++)
            Pistons[i].Inject(Y[i]);
        if (hasMore(Y[0]))
            for(unsigned int i=0; i<Pi; i++)
                Pistons[i].Spark();
    }
    delete[] Y;
}

ostream& operator<<(ostream& a, const Engine& engine)
{
    return a << "Engine[" << dec << engine.Pistons.size() << "\303\227" << engine.Pistons[0] << "]";
}

Motorist::Motorist(const Permutation *f, unsigned int Pi, unsigned int W, unsigned int c, unsigned int tau):
    Pi(Pi), W(W), c(c), cprime(W*((c + W - 1)/W)), tau(tau),
    Pistons(Pi, Piston(f, W/8*((f->getWidth() - max(c, (unsigned int)32))/W), W/8*((f->getWidth() - 32)/W))),
    engine(Pistons), phase(ready)
{
}

bool Motorist::StartEngine(istream& SUV, bool tagFlag, stringstream& T, bool decryptFlag, bool forgetFlag)
{
    if (phase != ready)
        throw Exception("The phase must be ready to call Motorist::StartEngine().");
    engine.InjectCollective(SUV, true);
    if (forgetFlag)
        MakeKnot();
    phase = riding;
    return HandleTag(tagFlag, T, decryptFlag);
}

bool Motorist::Wrap(istream& I, stringstream& O, istream& A, stringstream& T, bool decryptFlag, bool forgetFlag)
{
    if (phase != riding)
        throw Exception("The phase must be riding to call Motorist::Wrap().");
    do {
        engine.Wrap(I, O, A, decryptFlag);
    } while(hasMore(I) || hasMore(A));
    if ((Pi > 1) || forgetFlag)
        MakeKnot();
    bool res = HandleTag(true, T, decryptFlag);
    if (!res)
        O.str(string(""));
    return res;
}

void Motorist::MakeKnot(void)
{
    stringstream Tprime;
    engine.GetTags(Tprime, vector<unsigned int>(Pi, cprime/8));
    Tprime.seekg(0, ios_base::beg);
    engine.InjectCollective(Tprime, false);
}

bool Motorist::HandleTag(bool tagFlag, stringstream& T, bool decryptFlag)
{
    stringstream Tprime;
    if (!tagFlag)
        engine.GetTags(Tprime, vector<unsigned int>(Pi, 0));
    else {
        vector<unsigned int> l(Pi, 0);
        l[0] = tau/8;
        engine.GetTags(Tprime, l);
        if (!decryptFlag)
            T.str(Tprime.str());
        else if (Tprime.str() != T.str()) {
            phase = failed;
            return false;
        }
    }
    return true;
}

ostream& operator<<(ostream& a, const Motorist& motorist)
{
    return a << "Motorist[" << motorist.engine
        << ", W=" << dec << motorist.W << ", c=" << motorist.c
        << ", \317\204=" << motorist.tau << "]";
}
