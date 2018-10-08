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
 *
 * and related or neighboring rights to the source code in this file.
 * http://creativecommons.org/publicdomain/zero/1.0/
 */

#include <sstream>
#include <string>
#include <vector>
#include "bitstring.h"

using namespace std;

static UINT8 enc8(unsigned int x)
{
    if (x > 255) {
        stringstream s;
        s << "The integer " << x << " cannot be encoded on 8 bits.";
        throw Exception(s.str());
    }
    else
        return UINT8(x);
}

static void Farfalle_assert(bool condition, const string &synopsis, const char *fct)
{
    if ( !condition ) {
        throw Exception(string(fct) + "(): " + synopsis);
    }
}

#if defined(__GNUC__)
#define assert(cond, msg)  Farfalle_assert(cond, msg, __PRETTY_FUNCTION__)
#else
#define assert(cond, msg)  Farfalle_assert(cond, msg, __FUNCTION__)
#endif

void BitString::truncateLastByte(void)
{
    if ( vSize % 8 ) {
        v[vSize / 8] &= (1 << (vSize % 8)) - 1;                      // zeroize exceeding bits (for operator==())
    }
}

void BitString::syncAlias(void)
{
    if ( alias ) {
        *alias = str();
    }
}

BitString::BitString()
    : vSize(0), v(), alias(NULL)
{}

BitString::~BitString()
{
    syncAlias();                                                     // Our last hope if BitString was modified via array()
}

BitString::BitString(unsigned int bit)
    : vSize(1), v(1, bit), alias(NULL)
{
    assert((0 == bit) || (1 == bit), "bit must be 0 or 1.");
}

BitString::BitString(unsigned int size, UINT8 byte)
    : vSize(size), v((size + 7) / 8, byte), alias(NULL)
{
    truncateLastByte();
}

BitString::BitString(string &s)
    : vSize(s.size() * 8), v(), alias(&s)
{
    v.assign(s.c_str(), s.c_str() + s.size());
}

BitString::BitString(const string &s)
    : vSize(s.size() * 8), v(), alias(NULL)
{
    v.assign(s.c_str(), s.c_str() + s.size());
}

BitString::BitString(const string &s, unsigned int index, unsigned int size)
    : vSize((index >= s.size() * 8) ? 0 : (size + index <= s.size() * 8) ? size : s.size() * 8 - index),
    v(s.begin() + index / 8, s.begin() + index / 8 + ((vSize) + 7) / 8), // vSize must be initialized first! (better enable -Wreorder)
    alias(NULL)
{
    assert((index % 8) == 0, "This implementation only supports index that are multiple of 8.");
    truncateLastByte();
}

BitString::BitString(const BitString &S)
    : vSize(S.vSize), v(S.v), alias(NULL)                            // We don't copy the alias
{}

BitString::BitString(const BitString &S, unsigned int index, unsigned int size)
    : vSize((index >= S.vSize) ? 0 : (size + index <= S.vSize) ? size : S.vSize - index),
    v(S.v.begin() + index / 8, S.v.begin() + index / 8 + ((vSize) + 7) / 8),
    alias(NULL)
{
    assert((index % 8) == 0, "This implementation only supports index that are multiple of 8.");
    truncateLastByte();
}

BitString::BitString(const vector<UINT8> &v)
    : vSize(v.size() * 8), v(v), alias(NULL)
{}

BitString::BitString(const UINT8 *s, unsigned int size)
    : vSize(size), v(s, s + (size + 7) / 8), alias(NULL)
{}

string BitString::str() const
{
    return string(v.begin(), v.end());
}

UINT8 *BitString::array()
{
    assert((vSize % 8) == 0, "Can't get array if BitString length is not a multiple of 8."); // Because caller may modify the array and break the invariant
    return &v[0];
}

const UINT8 *BitString::array() const
{
    return &v[0];
}

unsigned int BitString::size() const
{
    return vSize;
}

BitString BitString::keypack(const BitString &K, unsigned int size)
{
    assert(0 < size,             "size must be positive.");
    assert((size % 8) == 0,      "Keypack length must be a multiple of 8.");
    assert(K.size() + 9 <= size, "The key is too big and does not fit in the key pack.");

    return BitString(8, enc8(size / 8)) || K || BitString::pad10(size - 8, K.size());
}

BitString BitString::substring(const BitString &K, unsigned int index, unsigned int size)
{
    unsigned int prependLength = (8 - index % 8) % 8;
    BitString K2 = BitString::zeroes(prependLength) || K;
    return BitString(K2, prependLength + index, size);
}

BitString BitString::pad10(unsigned int r, unsigned int Mlen)
{
    assert(0 < r, "r must be positive.");
    return BitString(1) || BitString::zeroes(r - 1 - (Mlen % r));
}

BitString BitString::pad101(unsigned int r, unsigned int Mlen)
{
    assert(0 < r, "r must be positive.");
    return BitString(1) || BitString::zeroes((2 * r - 2 - (Mlen % r)) % r) || BitString(1);
}

BitString BitString::zeroes(unsigned int size)
{
    return BitString(size, (UINT8)0);
}

BitString BitString::ones(unsigned int size)
{
    return BitString(size, (UINT8)255);
}

BitString &BitString::truncate(unsigned int size)
{
    if ( size > vSize ) {
        return *this;
    }
    vSize = size;
    v.resize((vSize + 7) / 8);
    truncateLastByte();
    syncAlias();
    return *this;
}

BitString &BitString::overwrite(const BitString &S, unsigned int index)
{
    assert((index % 8) == 0, "This implementation only supports index that are multiple of 8.");

    if ( index + S.vSize > vSize ) {
        vSize = index + S.vSize;
        v.resize((vSize + 7) / 8);
    }

    // Copy all complete bytes, i.e. (S.vSize/8) bytes, leaving (S.vSize%8)<8 bits left
    copy(S.v.begin(), S.v.begin() + (S.vSize / 8), v.begin() + (index / 8));

    // Copy the (S.vSize%8) remaining bits
    if ( S.vSize % 8 ) {
        UINT8  mask = (1 << (S.vSize % 8)) - 1;
        UINT8  src  = S.v[S.vSize / 8];
        UINT8 &dst  = v[(index / 8) + (S.vSize / 8)];
        dst = (dst & ~mask) | (src & mask);
    }

    syncAlias();

    return *this;
}

BitString &BitString::operator=(const BitString &A)
{
    if ( this != &A ) {
        vSize = A.vSize;
        v     = A.v;
    }
    syncAlias();
    return *this;
}

bool operator==(const BitString &A, const BitString &B)
{
    return (A.vSize == B.vSize) && (A.v == B.v);
}

BitString operator||(const BitString &A, unsigned int bit)
{
    assert((0 == bit) || (1 == bit), "bit must be 0 or 1.");

    BitString  Z(A);
    Z.v.resize((Z.vSize + 8) / 8);
    Z.v[Z.vSize / 8] |= bit << (Z.vSize % 8);
    Z.vSize++;

    return Z;
}

BitString operator||(const BitString &A, const BitString &B)
{
    // Allocate C
    BitString  C(A.vSize + B.vSize, 0);

    // Copy A into C
    copy(A.v.begin(), A.v.end(), C.v.begin());

    // Append B to C, starting from index A.vSize -- do it fast if possible, possibly with overflow
    if ((A.vSize % 8) == 0 ) {
        C.overwrite(B, A.vSize);
    }
    else {
        // There are (A.vSize%8) bits in last byte -- append with shift -- TODO: fast copy with 32-bit
        vector<UINT8>::iterator  c     = C.v.begin() + (A.vSize / 8);
        unsigned int             nbits = (A.vSize % 8);
        UINT8                    last  = *c & ((1 << (nbits)) - 1);  // last = -----xxx

        for ( vector<UINT8>::const_iterator b = B.v.begin(); b != B.v.end();) {
            *(c++) = last | (*b << nbits);                           // *c   = -----xxx | xxxxx---
            last   = *(b++) >> (8 - nbits);                          // last = -----xxx
        }

        if ( c != C.v.end()) {
            *(c++) = last;
        }
    }

    return C;
}

BitString operator^(const BitString &A, const BitString &B)
{
    assert(A.vSize == B.vSize, "Cannot xor two BitString of different size.");

    BitString                Z(A.vSize, (UINT8)0);
    vector<UINT8>::iterator  z = Z.v.begin();
    for ( vector<UINT8>::const_iterator a = A.v.begin(), b = B.v.begin(); a != A.v.end(); ++a, ++b, ++z ) {
        *z = *a ^ *b;
    }

    return Z;
}

ostream &operator<<(ostream &os, const BitString &S)
{
    for ( vector<UINT8>::const_iterator i = S.v.begin(); i != S.v.end();) {
        os.width(2);
        os.fill('0');
        os << hex << ((int)UINT8(*i));
        ++i;
        if ( i != S.v.end()) {
            os << " ";
        }
        else {
            os << "(" << ((S.vSize - 1) % 8 + 1) << ")";
        }
    }
    return os;
}

BitStrings::BitStrings()
{
}

BitStrings::BitStrings(const BitString &M)
{
	list.push_back(M);
}

size_t BitStrings::size() const
{
	return list.size();
}

const BitString &BitStrings::operator[](size_t i) const
{
	return list[i];
}

BitString &BitStrings::operator[](size_t i)
{
	return list[i];
}

BitStrings BitStrings::operator*(const BitString &M) const
{
	BitStrings tmp = M;
	tmp.list.insert(tmp.list.end(), list.begin(), list.end());
	return tmp;
}

BitStrings operator*(const BitString &M, const BitStrings &B)
{
	BitStrings tmp = B;
	tmp.list.push_back(M);
	return tmp;
}

BitStrings operator*(const BitString &A, const BitString &B)
{
	return A * BitStrings(B);
}

Block::Block(BitString &S, unsigned int index, unsigned int r)
    : B(S), alias(S), index(index), r(r)
{
    assert(0 < r,             "r must be positive.");
    assert(index <= S.size(), "index must be less than or equal to bit string size.");
}

static BitString  dummy; // Not nice but nobody will notice

Block::Block(const BitString &S, unsigned int index, unsigned int r)
    : B(dummy), alias(S), index(index), r(r)
{
    assert(0 < r,             "r must be positive.");
    assert(index <= S.size(), "index must be less than or equal to bit string size.");
}

Block &Block::operator=(const BitString &S)
{
    assert(S.size() <= r, "String size must be less than or equal to block size.");
    B.overwrite(S, index);

    return *this;
}

Block::operator BitString() const
{
    return BitString(alias, index, r);
}

unsigned int Block::size() const
{
    return alias.size() < (r + index) ? alias.size() - index : r;
}

ostream &operator<<(ostream &os, const Block &B)
{
    return os << "[" << BitString(B) << "]";
}

Blocks::Blocks(unsigned int r)
    : alias(B), const_alias(B), r(r)
{}

Blocks::Blocks(BitString &S, unsigned int r)
    : B(S), alias(S), const_alias(S), r(r)
{}

Blocks::Blocks(const BitString &S, unsigned int r)
    : B(S), alias(B), const_alias(S), r(r)                           // B(S) done to allow using non-const operator[]
{}

unsigned int Blocks::size() const
{
    return const_alias.size() > 0 ? (const_alias.size() + r - 1) / r : 1;
}

BitString Blocks::bits() const
{
    return const_alias;
}

Block Blocks::operator[](unsigned int i)
{
    return Block(alias, i * r, r);
}

Block Blocks::operator[](unsigned int i) const
{
    return Block(const_alias, i * r, r);
}

ostream &operator<<(ostream &os, const Blocks &B)
{
    for ( unsigned int i = 0; i < B.size(); ++i ) {
        os << B[i];
    }
    return os;
}
