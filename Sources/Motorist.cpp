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

Piston::Piston(const Permutation *aF, unsigned int aRs, unsigned int aRa)
    : f(aF), Rs(aRs), Ra(aRa)
{
    unsigned int b = f->getWidth();
    if (Rs > Ra)
        throw Exception("Rs is larger than Ra.");
    if (Ra > ((b-32)/8))
        throw Exception("Ra is larger than (b-32)/8.");
    state.reset(new UINT8[(b+7)/8]);
    memset(state.get(), 0, (b+7)/8);
    EOM = Ra;
    CryptEnd = Ra+1;
    InjectStart = Ra+2;
    InjectEnd = Ra+3;
}

Piston::Piston(const Piston& other)
    : f(other.f), Rs(other.Rs), Ra(other.Ra),
        EOM(other.EOM), CryptEnd(other.CryptEnd),
        InjectStart(other.InjectStart), InjectEnd(other.InjectEnd)
{
    unsigned int b = f->getWidth();
    state.reset(new UINT8[(b+7)/8]);
    memcpy(state.get(), other.state.get(), (b+7)/8);
}

void Piston::Crypt(istream& I, ostream& O, unsigned int omega, bool unwrapFlag)
{
    while(hasMore(I) && (omega < Rs)) {
        UINT8 x = I.get();
        O.put(state.get()[omega] ^ x);
        if (unwrapFlag)
            state.get()[omega] = x;
        else
            state.get()[omega] ^= x;
        omega++;
    }
    state.get()[CryptEnd] ^= enc8(omega);
}

void Piston::Inject(istream& X, bool cryptingFlag)
{
    unsigned int omega;
    if (cryptingFlag)
        omega = Rs;
    else
        omega = 0;
    state.get()[InjectStart] ^= enc8(omega);
    while(hasMore(X) && (omega < Ra)) {
        state.get()[omega] ^= X.get();
        omega++;
    }
    state.get()[InjectEnd] ^= enc8(omega);
}

void Piston::Spark(bool eomFlag, unsigned int l)
{
    if (eomFlag) {
        if (l == 0)
            state.get()[EOM] ^= enc8(255);
        else
            state.get()[EOM] ^= enc8(l);
    }
    else
        state.get()[EOM] ^= enc8(0);
    (*f)(state.get());
}

void Piston::GetTag(ostream& T, unsigned int l) const
{
    if (l > Rs)
        throw Exception("The requested tag is too long.");
    for(unsigned int i=0; i<l; i++)
        T.put(state.get()[i]);
}

ostream& operator<<(ostream& a, const Piston& piston)
{
    return a << "Piston[f=" << (*piston.f) << ", Rs=" << dec << piston.Rs << ", Ra=" << piston.Ra << "]";
}

Engine::Engine(vector<Piston>& aPistons)
    : Pi(aPistons.size()), Pistons(aPistons), phase(fresh), Et(Pi, 0)
{
}

void Engine::Spark(bool eomFlag, const vector<unsigned int>& l)
{
    for(unsigned int i=0; i<Pi; i++)
        Pistons[i].Spark(eomFlag, l[i]);
    Et = l;
}

void Engine::Crypt(istream& I, ostream& O, bool unwrapFlag)
{
    if (phase != fresh)
        throw Exception("The phase must be fresh to call Engine::Crypt().");
    for(unsigned int i=0; i<Pi; i++)
        Pistons[i].Crypt(I, O, Et[i], unwrapFlag);
    if (hasMore(I))
        phase = crypted;
    else
        phase = endOfCrypt;
}

void Engine::Inject(istream& A)
{
    if ((phase != fresh) && (phase != crypted) && (phase != endOfCrypt))
        throw Exception("The phase must be fresh, crypted or endOfCrypt to call Engine::Inject().");
    bool cryptingFlag = (phase == crypted) || (phase == endOfCrypt);
    for(unsigned int i=0; i<Pi; i++)
        Pistons[i].Inject(A, cryptingFlag);
    if ((phase == crypted) || hasMore(A)) {
        Spark(false, vector<unsigned int>(Pi, 0));
        phase = fresh;
    }
    else
        phase = endOfMessage;
}

void Engine::GetTags(ostream& T, const vector<unsigned int>& l)
{
    if (phase != endOfMessage)
        throw Exception("The phase must be endOfMessage to call Engine::GetTags().");
    Spark(true, l);
    for(unsigned int i=0; i<Pi; i++)
        Pistons[i].GetTag(T, l[i]);
    phase = fresh;
}

void Engine::InjectCollective(istream& X, bool diversifyFlag)
{
    if (phase != fresh)
        throw Exception("The phase must be fresh to call Engine::InjectCollective().");
    stringstream *Xt = new stringstream[Pi];
    while(hasMore(X)) {
        UINT8 x = X.get();
        for(unsigned int i=0; i<Pi; i++)
            Xt[i].put(x);
    }
    if (diversifyFlag) {
        for(unsigned int i=0; i<Pi; i++) {
            Xt[i].put(enc8(Pi));
            Xt[i].put(enc8(i));
        }
    }
    for(unsigned int i=0; i<Pi; i++)
        Xt[i].seekg(0, ios_base::beg);
    while(hasMore(Xt[0])) {
        for(unsigned int i=0; i<Pi; i++)
            Pistons[i].Inject(Xt[i], 0);
        if (hasMore(Xt[0]))
            Spark(false, vector<unsigned int>(Pi, 0));
    }
    phase = endOfMessage;
    delete[] Xt;
}

ostream& operator<<(ostream& a, const Engine& engine)
{
    return a << "Engine[" << engine.Pi << "\303\227" << engine.Pistons[0] << "]";
}

Motorist::Motorist(const Permutation *f, unsigned int aPi, unsigned int aW, unsigned int ac, unsigned int atau):
    Pi(aPi), W(aW), c(ac),
    Pistons(aPi, Piston(f, aW/8*((f->getWidth() - max(ac, (unsigned int)32))/aW), aW/8*((f->getWidth() - 32)/aW))),
    engine(Pistons), cprime(aW*((ac + aW - 1)/aW)), tau(atau), phase(ready)
{
}

bool Motorist::StartEngine(istream& SUV, bool tagFlag, stringstream& T, bool unwrapFlag, bool forgetFlag)
{
    if (phase != ready)
        throw Exception("The phase must be ready to call Motorist::StartEngine().");
    engine.InjectCollective(SUV, true);
    if (forgetFlag)
        MakeKnot();
    bool res = HandleTag(tagFlag, T, unwrapFlag);
    if (res)
        phase = riding;
    return res;
}

bool Motorist::Wrap(istream& I, stringstream& O, istream& A, stringstream& T, bool unwrapFlag, bool forgetFlag)
{
    if (phase != riding)
        throw Exception("The phase must be riding to call Motorist::Wrap().");
    if ((!hasMore(I)) && (!hasMore(A)))
        engine.Inject(A);
    while(hasMore(I)) {
        engine.Crypt(I, O, unwrapFlag);
        engine.Inject(A);
    }
    while(hasMore(A))
        engine.Inject(A);
    if ((Pi > 1) || forgetFlag)
        MakeKnot();
    bool res = HandleTag(true, T, unwrapFlag);
    if (!res)
        O.str(string(""));
    return res;
}

void Motorist::MakeKnot()
{
    stringstream Tprime;
    engine.GetTags(Tprime, vector<unsigned int>(Pi, cprime/8));
    Tprime.seekg(0, ios_base::beg);
    engine.InjectCollective(Tprime, false);
}

bool Motorist::HandleTag(bool tagFlag, stringstream& T, bool unwrapFlag)
{
    stringstream Tprime;
    if (!tagFlag)
        engine.GetTags(Tprime, vector<unsigned int>(Pi, 0));
    else {
        vector<unsigned int> l(Pi, 0);
        l[0] = tau/8;
        engine.GetTags(Tprime, l);
        if (!unwrapFlag)
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
