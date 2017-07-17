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
 * and related or neighboring rights to the source code in this file.
 * http://creativecommons.org/publicdomain/zero/1.0/
 */

#ifndef _BITSTRING_H_
#define _BITSTRING_H_

#include <iostream>
#include <vector>
#include "types.h"

using namespace std;

/**
 * Class implementing a simple bit string
 */
class BitString {
protected:
    unsigned int   vSize;                                            // size in bits -- invariant: v.size() == (vSize+7)/8
    vector<UINT8>  v;                                                // vector storing bytes -- invariant: if (vSize%8), then (v[vSize/8] >> (vSize%8)) == 0
    string *       alias;
    void  truncateLastByte(void);
    void  syncAlias(void);
public:
    BitString();
    ~BitString();
    BitString(unsigned int bit);
    BitString(unsigned int size, UINT8 byte);
    BitString(string &s);                                            // By reference -- BitString will update given string accordingly
    BitString(const string &s);
    BitString(const string &s, unsigned int index, unsigned int size);
    BitString(const BitString &S);
    BitString(const BitString &S, unsigned int index, unsigned int size);
    BitString(const vector<UINT8> &v);
    BitString(const UINT8 *s, unsigned int size);
    string            str() const;
    UINT8 *           array();
    const UINT8 *     array() const;
    unsigned int      size() const;
    static BitString  keypack(const BitString &K, unsigned int size);
    static BitString  substring(const BitString &K, unsigned int index, unsigned int size);
    static BitString  pad10(unsigned int r, unsigned int Mlen);
    static BitString  pad101(unsigned int r, unsigned int Mlen);
    static BitString  zeroes(unsigned int size);
    static BitString  ones(unsigned int size);
    BitString &       truncate(unsigned int size);
    BitString &       overwrite(const BitString &S, unsigned int index);
    BitString &       operator=(const BitString &A);
    friend bool       operator==(const BitString &A, const BitString &B);
    friend BitString  operator||(const BitString &A, unsigned int bit);
    friend BitString  operator||(const BitString &A, const BitString &B);
    friend BitString  operator^(const BitString &A, const BitString &B);
    friend ostream &  operator<<(ostream &os, const BitString &S);
};

bool                  operator==(const BitString &A, const BitString &B);
BitString             operator||(const BitString &A, unsigned int bit);
BitString             operator||(const BitString &A, const BitString &B);
BitString             operator^(const BitString &A, const BitString &B);
ostream &             operator<<(ostream &os, const BitString &S);

/**
 * Class implementing a string of bit strings
 */
class BitStrings
{
	private:
		std::vector<BitString> list;

	public:
		BitStrings();
		BitStrings(const BitString &M);
		size_t size() const;
		const BitString &operator[](size_t i) const;
		BitString &operator[](size_t i);
		BitStrings operator*(const BitString &M) const;
		friend BitStrings operator*(const BitString &M, const BitStrings &B);
};

BitStrings operator*(const BitString &A, const BitString &B);

/**
 * Class implementing a simple block of bits
 */
class Block {
protected:
    BitString &      B;
    const BitString &alias;                                          // Block is mutable if &B == &alias
    unsigned int     index;
    unsigned int     r;
public:
    Block(BitString &S, unsigned int index, unsigned int r);
    Block(const BitString &S, unsigned int index, unsigned int r);
    Block &         operator=(const BitString &S);
    operator BitString() const;
    unsigned int    size() const;
    friend ostream &operator<<(ostream &os, const Block &B);
};

ostream &           operator<<(ostream &os, const Block &B);

/**
 * Class implementing a string of bit blocks
 */
class Blocks {
protected:
    BitString        B;
    BitString &      alias;
    const BitString &const_alias;                                    // Block is mutable if &alias == &const_alias
    unsigned int     r;
public:
    Blocks(unsigned int r);                                          // Use internal BitString for storage
    Blocks(BitString &S, unsigned int r);                            // Use given BitString for storage
    Blocks(const BitString &S, unsigned int r);                      // Idem, but not mutable
    unsigned int    size() const;
    BitString       bits() const;
    Block           operator[](unsigned int i);
    Block           operator[](unsigned int i) const;
    friend ostream &operator<<(ostream &os, const Blocks &B);
};

#endif
