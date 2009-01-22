/*
Tools for the Keccak sponge function family.
Authors: Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche

This code is hereby put in the public domain. It is given as is, 
without any guarantee.

For more information, feedback or questions, please refer to our website:
http://keccak.noekeon.org/
*/

#ifndef _KECCAKF_H_
#define _KECCAKF_H_

#include <iostream>
#include <string>
#include <vector>
#include "transformations.h"

using namespace std;

/**
 * Exception that can be thrown by the classes Keccak and KeccakF.
 */
class KeccakException {
public:
    /** A string expressing the reason for the exception. */
    string reason;
    /**
     * The constructor.
     * @param   aReason     A string giving the reason of the exception.
     */
    KeccakException(const string& aReason);
};

/**
  * Class implementing the 7 Keccak-f permutations, as well as their inverses.
  */
class KeccakF : public Permutation {
protected:
    /** The width of the permutation. */
    unsigned int width;
    /** The size of the lanes. */
    unsigned int laneSize;
    /** The nominal number of rounds, as function of the width. */
    unsigned int nominalNrRounds;
    /** The actual number of rounds. */
    unsigned int nrRounds;
    /** The translation offsets for ρ. */
    vector<int> rhoOffsets; 
    /** The round constants for ι. */
    vector<UINT64> roundConstants;
    /** A 64-bit word whose first laneSize bits are 1 and all others 0. */
    UINT64 mask;
public:
    /**
      * The constructor. The width and number of rounds are given to the 
      * constructor. If omitted or set to zero, nrRounds is set to 
      * nominalNrRounds.
      *
      * @param  aWidth      The width of the Keccak-f permutation. 
      *                     It must be one of the valid Keccak-f widths, namely
      *                     25, 50, 100, 200, 400, 800 or 1600.
      * @param  aNrRounds   The desired number of rounds. By omitting or 
      *                     setting this parameter to 0, the nominal number
      *                     of rounds is taken.
      */
    KeccakF(unsigned int aWidth, unsigned int aNrRounds = 0);
    /**
      * Method that returns the number of bits of its domain and range.
      */
    unsigned int getWidth() const;
    /**
      * Method that applies the Keccak-f permutation onto the parameter 
      * @a state.
      */
    void operator()(UINT8 * state) const;
    /**
      * Method that applies the inverse of the Keccak-f permutation onto 
      * the parameter @a state.
      */
    void inverse(UINT8 * state) const;
    /**
      * Method that returns a string describing the instance of the Keccak-f 
      * permutation.
      */
    string getDescription() const;
    /**
      * Method that returns a short string that uniquely identifies the
      * Keccak-f instance.
      */
    string getName() const;
    /** 
      * Method that builds a file name by prepending a prefix and appending 
      * a suffix to getName().
      */
    string buildFileName(const string& prefix, const string& suffix) const;
    /**
      * Method that maps the coordinates (x, y) onto the lanes 
      * numbered from 0 to 24. The formula is (x mod 5)+5*(y mod 5), 
      * so that the lanes are ordered in line with the bit ordering defined 
      * in the specifications.
      *
      * @param  x           The x coordinate. It can be any signed integer,
      *                     as it will be reduced modulo 5.
      * @param  y           The y coordinate. It can be any signed integer,
      *                     as it will be reduced modulo 5.
      */
    static unsigned int index(int x, int y);
    /**
      * Method that reduces modulo 5 the coordinate x expressed as a signed
      * integer.
      *
      * @param  x           The x coordinate. It can be any signed integer,
      *                     as it will be reduced modulo 5.
      */
    static unsigned int index(int x);
    /**
      * Template method that applies the permutation.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  state   The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a T can be a UINT64 for an
      *                 actual evaluation.
      */
    template<class T> void forward(vector<T>& state) const;
    /**
      * Template method that applies the inverse permutation.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  state   The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a T can be a UINT64 for an
      *                 actual evaluation.
      */
    template<class T> void inverse(vector<T>& state) const;
    /**
      * Template method that applies the round function.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a T can be a UINT64 for an
      *                 actual evaluation.
      * @param  roundNumber     The round number, from 0 to nrRounds - 1.
      */
    template<class T> void round(vector<T>& A, unsigned int roundNumber) const;
    /**
      * Template method that applies the inverse of the round function. 
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a T can be a UINT64 for an
      *                 actual evaluation.
      * @param  roundNumber     The round number, from 0 to nrRounds - 1.
      */
    template<class T> void inverseRound(vector<T>& A, unsigned int roundNumber) const;
    /**
      * Template method that applies χ.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a T can be a UINT64 for an
      *                 actual evaluation.
      */
    template<class T> void chi(vector<T>& A) const;
    /**
      * Template method that applies the inverse of χ.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a T can be a UINT64 for an
      *                 actual evaluation.
      */
    template<class T> void inverseChi(vector<T>& A) const;
    /**
      * Template method that applies θ.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a T can be a UINT64 for an
      *                 actual evaluation.
      */
    template<class T> void theta(vector<T>& A) const;
    /**
      * Template method that applies the inverse of θ.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a T can be a UINT64 for an
      *                 actual evaluation.
      */
    template<class T> void inverseTheta(vector<T>& A) const;
    /**
      * Template method that applies π.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a T can be a UINT64 for an
      *                 actual evaluation.
      */
    template<class T> void pi(vector<T>& A) const;
    /**
      * Method that applies the π coordinate transformation to a lane position (x,y).
      *
      * @param[in]  x   The input coordinate x.
      * @param[in]  y   The input coordinate y.
      * @param[out] X   The output coordinate x.
      * @param[out] Y   The output coordinate y.
      */
    void pi(unsigned int x, unsigned int y, unsigned int& X, unsigned int& Y) const;
    /**
      * Template method that applies the inverse of π.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a T can be a UINT64 for an
      *                 actual evaluation.
      */
    template<class T> void inversePi(vector<T>& A) const;
    /**
      * Method that applies the inverse π coordinate transformation to a lane position (X,Y).
      *
      * @param[in]  X   The input coordinate x.
      * @param[in]  Y   The input coordinate y.
      * @param[out] x   The output coordinate x.
      * @param[out] y   The output coordinate y.
      */
    void inversePi(unsigned int X, unsigned int Y, unsigned int& x, unsigned int& y) const;
    /**
      * Template method that applies ρ.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a T can be a UINT64 for an
      *                 actual evaluation.
      */
    template<class T> void rho(vector<T>& A) const;
    /**
      * Template method that applies the inverse of ρ.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a T can be a UINT64 for an
      *                 actual evaluation.
      */
    template<class T> void inverseRho(vector<T>& A) const;
    /**
      * Template method that applies ι, which is its own inverse.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a T can be a UINT64 for an
      *                 actual evaluation.
      * @param  roundNumber     The round number, from 0 to nrRounds - 1.
      */
    template<class T> void iota(vector<T>& A, unsigned int roundNumber) const;
    /**
      * Template method that translates a lane along the z-axis.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  L       The given lane.
      * @param offset   The translation offset. It can be any signed
      *                 integer, as it will be reduced modulo laneSize.
      */
    template<class T> void ROL(T& L, int offset) const;
    /**
      * Method that implementats ROL when the lane is in a 64-bit word UINT64.
      *
      * @param  L       The given lane.
      * @param offset   The translation offset. It can be any signed
      *                 integer, as it will be reduced modulo laneSize.
      */
    void ROL(UINT64& L, int offset) const;
    /**
      * Method that converts a state given as an array of bytes into a vector
      * of lanes in 64-bit words.
      *
      * @param  in      The state as an array of bytes.
      *                 The array @a in must have a size of at least 
      *                 ceil(getWidth()/8.0) bytes.
      * @param  out     The state as a vector of lanes.
      *                 It will be resized to 25 if necessary.
      */
    void fromBytesToLanes(UINT8 *in, vector<UINT64>& out) const;
    /**
      * Method that converts a vector of lanes in 64-bit words into a state 
      * given as an array of bytes.
      *
      * @param  in      The state as a vector of lanes.
      *                 in.size() must be equal to 25.
      * @param  out     The state as an array of bytes.
      *                 The array @a out must have a size of at least 
      *                 ceil(getWidth()/8.0) bytes.
      */
    void fromLanesToBytes(const vector<UINT64>& in, UINT8 *out) const;
    /**
      * Method that initializes the nrRounds round constants according to the 
      * specifications.
      */
    void initializeRoundConstants();
    /**
      * Method that initializes the 25 lane translation offsets for ρ according to 
      * the specifications.
      */
    void initializeRhoOffsets();
    /**
      * Function that appends the z coordinate to the given prefix. 
      * If the lane size is 1, the z coordinate is not appended.
      * If the lane size is at most 10, the z coordinate is appended using one digit.
      * Otherwise, the z coordinate is appended using two digits.
          *
          * @param  prefixSymbol  The name prefix.
          * @param laneSize The lane size.
          * @param z        The z coordinate.
      */
    static string buildBitName(const string& prefixSymbol, unsigned int laneSize, unsigned int z);
    /**
      * Method that constructs a variable name for a particular bit in the state.
      * The resulting variable name consists of the prefix, a consonant coding the y coordinate, 
      * a vowel coding the x coordinate and finally a digit (or pair of digits) coding the z coordinate 
      * (using buildBitName(), i.e., if the lane size is greater than 1).
      * This naming convention is such that the alphabetic order is also
      * the bit ordering for the message input at sponge level.
      *
      * @param  prefix  The name prefix.
      * @param x        The x coordinate.
      * @param y        The y coordinate.
      * @param z        The z coordinate.
      */
    string bitName(const string& prefix, unsigned int x, unsigned int y, unsigned int z) const;
    /**
      * Method that constructs a variable name for a particular lane in the state.
      * The resulting variable name consists of the prefix, a consonant coding the y coordinate 
      * and a vowel coding the x coordinate.
      *
      * @param  prefix  The name prefix.
      * @param x        The x coordinate.
      * @param y        The y coordinate.
      */
    static string laneName(const string& prefix, unsigned int x, unsigned int y);
    /**
      * Method that constructs a variable name for a particular sheet in the state.
      * The resulting variable name consists of the prefix and a vowel coding the x coordinate.
      *
      * @param  prefix  The name prefix.
      * @param x        The x coordinate.
      */
    static string sheetName(const string& prefix, unsigned int x);
};

template<class T>
void KeccakF::forward(vector<T>& state) const
{
    for(unsigned int i=0; i<nrRounds; i++)
        round(state, i);
}

template<class T>
void KeccakF::inverse(vector<T>& state) const
{
    for(int i=nrRounds-1; i>=0; i--)
        inverseRound(state, i);
}

template<class T>
void KeccakF::round(vector<T>& state, unsigned int roundNumber) const
{
    theta(state);
    rho(state);
    pi(state);
    chi(state);
    iota(state, roundNumber);
}

template<class T>
void KeccakF::inverseRound(vector<T>& state, unsigned int roundNumber) const
{
    iota(state, roundNumber);
    inverseChi(state);
    inversePi(state);
    inverseRho(state);
    inverseTheta(state);
}

template<class T>
void KeccakF::chi(vector<T>& A) const
{
    vector<T> C(5);
    for(unsigned int y=0; y<5; y++) { 
        for(unsigned int x=0; x<5; x++)
            C[x] = A[index(x,y)] ^ ((~A[index(x+1,y)]) & A[index(x+2,y)]);
        for(unsigned int x=0; x<5; x++)
            A[index(x,y)] = C[index(x)];
    }
}

template<class T>
void KeccakF::inverseChi(vector<T>& A) const
{
    vector<T> C(5);
    for(unsigned int y=0; y<5; y++) { 
        for(unsigned int x=0; x<5; x++) 
            C[x] = A[index(x,y)];
        for(unsigned int x=0; x<(5+1); x++) {
            unsigned int X = x*3;
            A[index(X,y)] = C[index(X)] ^ (A[index(X+2,y)] & (~C[index(X+1)]));
        }
    }
}

template<class T>
void KeccakF::theta(vector<T>& A) const
{
    vector<T> C(5);
    for(unsigned int x=0; x<5; x++) {
        C[x] = A[index(x,0)]; 
        for(unsigned int y=1; y<5; y++) 
            C[x] ^= A[index(x,y)];
    }
    vector<T> D(C);
    for(unsigned int x=0; x<5; x++)
        ROL(D[x], 1);
    for(int x=0; x<5; x++)
        for(unsigned int y=0; y<5; y++)
            A[index(x,y)] ^= D[index(x+1)] ^ C[index(x-1)];
}

template<class T>
void KeccakF::inverseTheta(vector<T>& A) const
{
    vector<T> C(5);
    for(unsigned int x=0; x<5; x++) {
        C[x] = A[index(x,0)];
        for(unsigned int y=1; y<5; y++){ 
            C[x] ^= A[index(x,y)];
        }
    }
    const UINT64 inversePositions64[5] = {
        0xDE26BC4D789AF134ULL,
        0x09AF135E26BC4D78ULL,
        0xEBC4D789AF135E26ULL,
        0x7135E26BC4D789AFULL,
        0xCD789AF135E26BC4ULL };
    vector<UINT64> inversePositions(5, 0);
    for(unsigned int z=0; z<64; z+=laneSize)
        for(unsigned int x=0; x<5; x++)
            inversePositions[x] ^= inversePositions64[x] >> z;
    for(unsigned int z=0; z<laneSize; z++) {
        for(unsigned int xOff=0; xOff<5; xOff++)
           for(int x=0; x<5; x++)
               for(unsigned int y=0; y<5; y++)
                   if ((inversePositions[xOff] & 1) != 0) 
                      A[index(x, y)] ^= C[index(x-xOff)];
        for(unsigned int xOff=0; xOff<5; xOff++) {
            ROL(C[xOff], 1);
            inversePositions[xOff] >>= 1;
        }
    }
}

template<class T>
void KeccakF::pi(vector<T>& A) const
{
    vector<T> a(A);
    for(unsigned int x=0; x<5; x++) 
    for(unsigned int y=0; y<5; y++) {
        unsigned int X, Y;
        pi(x, y, X, Y);
        A[index(X,Y)] = a[index(x,y)];
    }
}

template<class T>
void KeccakF::inversePi(vector<T>& A) const
{
    vector<T> a(A);
    for(unsigned int X=0; X<5; X++) 
    for(unsigned int Y=0; Y<5; Y++) {
        unsigned int x, y;
        inversePi(X, Y, x, y);
        A[index(x,y)] = a[index(X,Y)];
    }
}

template<class T>
void KeccakF::rho(vector<T>& A) const
{
    for(unsigned int x=0; x<5; x++) 
    for(unsigned int y=0; y<5; y++)
        ROL(A[index(x,y)], rhoOffsets[index(x,y)]);
}

template<class T>
void KeccakF::inverseRho(vector<T>& A) const
{
    for(unsigned int x=0; x<5; x++) 
    for(unsigned int y=0; y<5; y++)
        ROL(A[index(x,y)], -rhoOffsets[index(x,y)]);
}

template<class T>
void KeccakF::iota(vector<T>& A, unsigned int roundNumber) const
{
    if (roundNumber < roundConstants.size())
        A[index(0,0)] ^= (roundConstants[roundNumber] & mask);
}

template<class T> 
void KeccakF::ROL(T& L, int offset) const
{
    L.ROL(offset, laneSize);
}

#endif
