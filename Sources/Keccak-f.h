/*
KeccakTools

The Keccak sponge function, designed by Guido Bertoni, Joan Daemen,
Michaël Peeters and Gilles Van Assche. For more information, feedback or
questions, please refer to our website: http://keccak.noekeon.org/

Implementation by the designers,
hereby denoted as "the implementer".

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/
*/

#ifndef _KECCAKF_H_
#define _KECCAKF_H_

#include <iostream>
#include <string>
#include <vector>
#include "transformations.h"

using namespace std;

typedef Exception KeccakException;

typedef UINT64 LaneValue;

/**
  * Class implementing the 7 Keccak-<i>f</i> permutations, as well as their inverses.
  */
class KeccakF : public Permutation {
protected:
    /** The width of the permutation. */
    unsigned int width;
    /** The size of the lanes. */
    unsigned int laneSize;
    /** The nominal number of rounds, as function of the width. */
    unsigned int nominalNrRounds;
    /** The round index to start with (for experiments on reduced-round version). */
    int startRoundIndex;
    /** The actual number of rounds (for experiments on reduced-round version). */
    unsigned int nrRounds;
    /** The translation offsets for ρ. */
    vector<int> rhoOffsets; 
    /** The round constants for ι. */
    vector<LaneValue> roundConstants;
    /** A 64-bit word whose first laneSize bits are 1 and all others 0. */
    LaneValue mask;
protected:
    /**
      * The constructor. The width and the range of rounds are
      * given as parameters.
      *
      * @param  aWidth      The width of the Keccak-<i>f</i> permutation.
      *                     It must be one of the valid Keccak-<i>f</i> widths, namely
      *                     25, 50, 100, 200, 400, 800 or 1600.
      * @param  aStartRoundIndex    The index of the first round to perform.
      * @param  aNrRounds   The desired number of rounds.
      */
    KeccakF(unsigned int aWidth, int aStartRoundIndex, unsigned int aNrRounds);
public:
    /**
      * The constructor. The width and number of rounds are given to the 
      * constructor. The number of rounds is the nominal number of rounds
      * as defined by the specifications.
      *
      * @param  aWidth      The width of the Keccak-<i>f</i> permutation. 
      *                     It must be one of the valid Keccak-<i>f</i> widths, namely
      *                     25, 50, 100, 200, 400, 800 or 1600.
      */
    KeccakF(unsigned int aWidth);
    /**
      * Method that returns the number of bits of its domain and range.
      */
    unsigned int getWidth() const;
    /** 
      * Method that retuns the lane size of the Keccak-<i>f</i> instance.
      */
    unsigned int getLaneSize() const;
    /**
      * Method that returns the number of rounds of this instance.
      * This should be the nominal number of rounds for Keccak-<i>f</i> itself.
      */
    unsigned int getNumberOfRounds() const;
    /**
      * Method that returns the nominal number of rounds of Keccak-<i>f</i>
      * according to the specifications for the width of this instance.
      */
    unsigned int getNominalNumberOfRounds() const;
    /**
      * Method that returns the index of the first round of this instance.
      * This should be zero for Keccak-<i>f</i> itself.
      */
    int getIndexOfFirstRound() const;
    /**
      * Method that applies the Keccak-<i>f</i> permutation onto the parameter 
      * @a state.
      */
    void operator()(UINT8 * state) const;
    /**
      * Method that applies the inverse of the Keccak-<i>f</i> permutation onto 
      * the parameter @a state.
      */
    void inverse(UINT8 * state) const;
    /**
      * Method that returns a string describing the instance of the Keccak-<i>f</i> 
      * permutation.
      */
    string getDescription() const;
    /**
      * Method that returns a short string that uniquely identifies the
      * Keccak-<i>f</i> instance.
      */
    virtual string getName() const;
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
      * Method that extracts the x coordinate from a lane numbered according
      * to the index() method.
      *
      * @param  index       The index of a lane, between 0 and 24.
      */
    inline static unsigned int getX(unsigned int index)
    {
        return index % 5;
    }
    /**
      * Method that extracts the y coordinate from a lane numbered according
      * to the index() method.
      *
      * @param  index       The index of a lane, between 0 and 24.
      */
    inline static unsigned int getY(unsigned int index)
    {
        return index / 5;
    }
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
      *                 The parameter type @a Lane can be a LaneValue for an
      *                 actual evaluation.
      */
    template<class Lane> void forward(vector<Lane>& state) const;
    /**
      * Template method that applies the inverse permutation.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  state   The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a Lane can be a LaneValue for an
      *                 actual evaluation.
      */
    template<class Lane> void inverse(vector<Lane>& state) const;
    /**
      * Template method that applies the round function.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a Lane can be a LaneValue for an
      *                 actual evaluation.
      * @param  roundIndex  The round index.
      */
    template<class Lane> void round(vector<Lane>& A, int roundIndex) const;
    /**
      * Template method that applies the inverse of the round function. 
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a Lane can be a LaneValue for an
      *                 actual evaluation.
      * @param  roundIndex  The round index.
      */
    template<class Lane> void inverseRound(vector<Lane>& A, int roundIndex) const;
    /**
      * Template method that applies χ.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a Lane can be a LaneValue for an
      *                 actual evaluation.
      */
    template<class Lane> void chi(vector<Lane>& A) const;
    /**
      * Template method that applies the inverse of χ.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a Lane can be a LaneValue for an
      *                 actual evaluation.
      */
    template<class Lane> void inverseChi(vector<Lane>& A) const;
    /**
      * Template method that applies θ.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a Lane can be a LaneValue for an
      *                 actual evaluation.
      */
    template<class Lane> void theta(vector<Lane>& A) const;
    /**
      * Template method that applies the inverse of θ.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a Lane can be a LaneValue for an
      *                 actual evaluation.
      */
    template<class Lane> void inverseTheta(vector<Lane>& A) const;
    /**
      * Template method that applies π.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a Lane can be a LaneValue for an
      *                 actual evaluation.
      */
    template<class Lane> void pi(vector<Lane>& A) const;
    /**
      * Method that applies the π coordinate transformation to a lane position (x,y).
      *
      * @param[in]  x   The input coordinate x.
      * @param[in]  y   The input coordinate y.
      * @param[out] X   The output coordinate x.
      * @param[out] Y   The output coordinate y.
      */
    static void pi(unsigned int x, unsigned int y, unsigned int& X, unsigned int& Y);
    /**
      * Template method that applies the inverse of π.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a Lane can be a LaneValue for an
      *                 actual evaluation.
      */
    template<class Lane> void inversePi(vector<Lane>& A) const;
    /**
      * Method that applies the inverse π coordinate transformation to a lane position (X,Y).
      *
      * @param[in]  X   The input coordinate x.
      * @param[in]  Y   The input coordinate y.
      * @param[out] x   The output coordinate x.
      * @param[out] y   The output coordinate y.
      */
    static void inversePi(unsigned int X, unsigned int Y, unsigned int& x, unsigned int& y);
    /**
      * Template method that applies ρ.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a Lane can be a LaneValue for an
      *                 actual evaluation.
      */
    template<class Lane> void rho(vector<Lane>& A) const;
    /**
      * Method that applies the ρ coordinate transformation to a bit position (x,y,z).
      *
      * @param[in]  x   The input coordinate x.
      * @param[in]  y   The input coordinate y.
      * @param[in]  z   The input coordinate z.
      * @return         The output coordinate z.
      */
    inline unsigned int rho(unsigned int x, unsigned int y, unsigned int z) const
    {
        return (z+rhoOffsets[index(x,y)])%laneSize;
    }
    /**
      * Template method that applies the inverse of ρ.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a Lane can be a LaneValue for an
      *                 actual evaluation.
      */
    template<class Lane> void inverseRho(vector<Lane>& A) const;
    /**
      * Method that applies the inverse ρ coordinate transformation to a bit position (x,y,z).
      *
      * @param[in]  x   The input coordinate x.
      * @param[in]  y   The input coordinate y.
      * @param[in]  z   The input coordinate z.
      * @return         The output coordinate z.
      */
    inline unsigned int inverseRho(unsigned int x, unsigned int y, unsigned int z) const
    {
        return (640+z-rhoOffsets[index(x,y)])%laneSize;
    }
    /**
      * Template method that applies ι, which is its own inverse.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  A       The given state is organized as a vector of 25 lanes.
      *                 The mapping between the coordinates (x, y) and the
      *                 numbering inside the vector is defined by index().
      *                 The parameter type @a Lane can be a LaneValue for an
      *                 actual evaluation.
      * @param  roundIndex  The round index.
      */
    template<class Lane> void iota(vector<Lane>& A, int roundIndex) const;
    /** 
      * Method that retuns the i-th round constant used by ι.
      * 
      * @param  roundIndex  The round index.
      */
    LaneValue getRoundConstant(int roundIndex) const;
    /**
      * Template method that translates a lane along the z-axis.
      * This function is a template to allow both numerical and symbolic values 
      * to be processed (see also KeccakFEquations).
      *
      * @param  L       The given lane.
      * @param offset   The translation offset. It can be any signed
      *                 integer, as it will be reduced modulo laneSize.
      */
    template<class Lane> void ROL(Lane& L, int offset) const;
    /**
      * Method that implementats ROL when the lane is in a 64-bit word LaneValue.
      *
      * @param  L       The given lane.
      * @param offset   The translation offset. It can be any signed
      *                 integer, as it will be reduced modulo laneSize.
      */
    void ROL(LaneValue& L, int offset) const;
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
    void fromBytesToLanes(const UINT8 *in, vector<LaneValue>& out) const;
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
    void fromLanesToBytes(const vector<LaneValue>& in, UINT8 *out) const;
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
private:
    /**
      * Method that initializes the nominal number of rounds according to the
      * specifications.
      */
    void initializeNominalNumberOfRounds();
    /**
      * Method that initializes the round constants according to the
      * specifications.
      */
    void initializeRoundConstants();
    /**
      * Method that initializes the 25 lane translation offsets for ρ according to
      * the specifications.
      */
    void initializeRhoOffsets();
};

template<class Lane>
void KeccakF::forward(vector<Lane>& state) const
{
    for(int i=startRoundIndex; i<startRoundIndex+nrRounds; i++)
        round(state, i);
}

template<class Lane>
void KeccakF::inverse(vector<Lane>& state) const
{
    for(int i=startRoundIndex+nrRounds-1; i>=startRoundIndex; i--)
        inverseRound(state, i);
}

template<class Lane>
void KeccakF::round(vector<Lane>& state, int roundIndex) const
{
    theta(state);
    rho(state);
    pi(state);
    chi(state);
    iota(state, roundIndex);
}

template<class Lane>
void KeccakF::inverseRound(vector<Lane>& state, int roundIndex) const
{
    iota(state, roundIndex);
    inverseChi(state);
    inversePi(state);
    inverseRho(state);
    inverseTheta(state);
}

template<class Lane>
void KeccakF::chi(vector<Lane>& A) const
{
    vector<Lane> C(5);
    for(unsigned int y=0; y<5; y++) { 
        for(unsigned int x=0; x<5; x++)
            C[x] = A[index(x,y)] ^ ((~A[index(x+1,y)]) & A[index(x+2,y)]);
        for(unsigned int x=0; x<5; x++)
            A[index(x,y)] = C[index(x)];
    }
}


template<class Lane>
void KeccakF::inverseChi(vector<Lane>& A) const
{
    for(unsigned int y=0; y<5; y++) { 
        unsigned int length = 5;
        vector<Lane> C(length);
        for(unsigned int x=0; x<length; x++) C[index(x)] = A[index(x,y)];
        for(unsigned int x=0; x<(3*(length-1)/2); x++) {
            unsigned int X = (length-2)*x;
            A[index(X,y)] = C[index(X)] ^ (A[index(X+2,y)] & (~C[index(X+1)]));
        }
    }
}

template<class Lane>
void KeccakF::theta(vector<Lane>& A) const
{
    vector<Lane> C(5);
    for(unsigned int x=0; x<5; x++) {
        C[x] = A[index(x,0)]; 
        for(unsigned int y=1; y<5; y++) 
            C[x] ^= A[index(x,y)];
    }
    vector<Lane> D(5);
    for(unsigned int x=0; x<5; x++) {
        Lane temp = C[index(x+1)];
        ROL(temp, 1);
        D[x] = temp ^ C[index(x-1)];
    }
    for(int x=0; x<5; x++)
        for(unsigned int y=0; y<5; y++)
            A[index(x,y)] ^= D[x];
}

template<class Lane>
void KeccakF::inverseTheta(vector<Lane>& A) const
{
    vector<Lane> C(5);
    for(unsigned int x=0; x<5; x++) {
        C[x] = A[index(x,0)];
        for(unsigned int y=1; y<5; y++){ 
            C[x] ^= A[index(x,y)];
        }
    }
    const LaneValue inversePositions64[5] = {
        0xDE26BC4D789AF134ULL,
        0x09AF135E26BC4D78ULL,
        0xEBC4D789AF135E26ULL,
        0x7135E26BC4D789AFULL,
        0xCD789AF135E26BC4ULL };
    vector<LaneValue> inversePositions(5, 0);
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

template<class Lane>
void KeccakF::pi(vector<Lane>& A) const
{
    vector<Lane> a(A);
    for(unsigned int x=0; x<5; x++) 
    for(unsigned int y=0; y<5; y++) {
        unsigned int X, Y;
        pi(x, y, X, Y);
        A[index(X,Y)] = a[index(x,y)];
    }
}

template<class Lane>
void KeccakF::inversePi(vector<Lane>& A) const
{
    vector<Lane> a(A);
    for(unsigned int X=0; X<5; X++) 
    for(unsigned int Y=0; Y<5; Y++) {
        unsigned int x, y;
        inversePi(X, Y, x, y);
        A[index(x,y)] = a[index(X,Y)];
    }
}

template<class Lane>
void KeccakF::rho(vector<Lane>& A) const
{
    for(unsigned int x=0; x<5; x++) 
    for(unsigned int y=0; y<5; y++)
        ROL(A[index(x,y)], rhoOffsets[index(x,y)]);
}

template<class Lane>
void KeccakF::inverseRho(vector<Lane>& A) const
{
    for(unsigned int x=0; x<5; x++) 
    for(unsigned int y=0; y<5; y++)
        ROL(A[index(x,y)], -rhoOffsets[index(x,y)]);
}

template<class Lane>
void KeccakF::iota(vector<Lane>& A, int roundIndex) const
{
    unsigned int ir = (unsigned int)(((roundIndex % 255) + 255) % 255);
    A[index(0,0)] ^= roundConstants[ir];
}

template<class Lane> 
void KeccakF::ROL(Lane& L, int offset) const
{
    L.ROL(offset, laneSize);
}

/**
  * Class implementing the 7 Keccak-<i>f</i> permutations with a reduced
  * number of rounds, starting from the first nominal round.
  */
class KeccakFfirstRounds : public KeccakF {
public:
    /**
      * The constructor. The width and the number of rounds are
      * given as parameters.
      *
      * @param  aWidth      The width of the Keccak-<i>f</i> permutation.
      *                     It must be one of the valid Keccak-<i>f</i> widths, namely
      *                     25, 50, 100, 200, 400, 800 or 1600.
      * @param  aNrRounds   The desired number of rounds.
      */
    KeccakFfirstRounds(unsigned int aWidth, unsigned int aNrRounds);
    /**
      * The constructor. The width is given as parameter. The number of rounds
      * is nominal.
      *
      * @param  aWidth      The width of the Keccak-<i>f</i> permutation.
      *                     It must be one of the valid Keccak-<i>f</i> widths, namely
      *                     25, 50, 100, 200, 400, 800 or 1600.
      */
    KeccakFfirstRounds(unsigned int aWidth);
};

/**
  * Class implementing the 7 Keccak-<i>f</i> permutations with a reduced
  * number of rounds, ending at the last nominal round.
  */
class KeccakFlastRounds : public KeccakF {
public:
    /**
      * The constructor. The width and the number of rounds are
      * given as parameters.
      *
      * @param  aWidth      The width of the Keccak-<i>f</i> permutation.
      *                     It must be one of the valid Keccak-<i>f</i> widths, namely
      *                     25, 50, 100, 200, 400, 800 or 1600.
      * @param  aNrRounds   The desired number of rounds.
      */
    KeccakFlastRounds(unsigned int aWidth, unsigned int aNrRounds);
    /**
      * The constructor. The width is given as parameter. The number of rounds
      * is nominal.
      *
      * @param  aWidth      The width of the Keccak-<i>f</i> permutation.
      *                     It must be one of the valid Keccak-<i>f</i> widths, namely
      *                     25, 50, 100, 200, 400, 800 or 1600.
      */
    KeccakFlastRounds(unsigned int aWidth);
    string getName() const;
};

/**
  * Class implementing the 7 Keccak-<i>f</i> permutations with a reduced
  * number of rounds, starting at any round index.
  */
class KeccakFanyRounds : public KeccakF {
public:
    /**
      * The constructor. The width and the number of rounds are
      * given as parameters.
      *
      * @param  aWidth      The width of the Keccak-<i>f</i> permutation.
      *                     It must be one of the valid Keccak-<i>f</i> widths, namely
      *                     25, 50, 100, 200, 400, 800 or 1600.
      * @param  aStartRoundIndex    The index of the first round to perform.
      * @param  aNrRounds   The desired number of rounds.
      */
    KeccakFanyRounds(unsigned int aWidth, int aStartRoundIndex, unsigned int aNrRounds);
    /**
      * The constructor. The width is given as parameter. The number of rounds
      * is nominal.
      *
      * @param  aWidth      The width of the Keccak-<i>f</i> permutation.
      *                     It must be one of the valid Keccak-<i>f</i> widths, namely
      *                     25, 50, 100, 200, 400, 800 or 1600.
      */
    KeccakFanyRounds(unsigned int aWidth);
};

#endif
