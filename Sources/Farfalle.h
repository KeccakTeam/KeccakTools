/*
 * Implementation by the Farfalle and Kravatte Teams, namely, Guido Bertoni,
 * Joan Daemen, Seth Hoffert, Michaël Peeters, Gilles Van Assche and Ronny Van Keer,
 * hereby denoted as "the implementer".
 *
 * For more information, feedback or questions, please refer to our websites:
 * https://keccak.team/farfalle.html
 * https://keccak.team/kravatte.html
 *
 * To the extent possible under law, the implementer has waived all copyright
 * and related or neighboring rights to the source code in this file.
 * http://creativecommons.org/publicdomain/zero/1.0/
 */

#ifndef _FARFALLE_H_
#define _FARFALLE_H_

#include <iostream>
#include <memory>
#include "padding.h"
#include "transformations.h"
#include "types.h"
#include "bitstring.h"

using namespace std;

/**
 * Class implementing an iterable permutation
 */
class BaseIterableTransformation
{
	public:
		const unsigned int            width;
		const unsigned int            rounds;

		BaseIterableTransformation(unsigned int width, unsigned int rounds) : width(width), rounds(rounds) {}

		virtual BitString operator()(const BitString &state) const = 0;
};

template<class T>
class IterableTransformation : public BaseIterableTransformation
{
	protected:
		T                     f;

	public:
		IterableTransformation(unsigned int width) : BaseIterableTransformation(width, 0), f(width) {}
		IterableTransformation(unsigned int width, unsigned int rounds) : BaseIterableTransformation(width, rounds), f(width, rounds) {}

		BitString operator()(const BitString &state) const
		{
			BitString state2 = state;
			f(state2.array());
			return state2;
		}
};

/**
 * Class implementing a rolling function
 */
class BaseRollingFunction
{
	public:
		virtual BitString operator()(const BitString &k, unsigned int i) const = 0;
};

class IdentityRollingFunction : public BaseRollingFunction
{
	public:
		BitString operator()(const BitString &k, unsigned int i) const;
};

typedef Exception  FarException;

/**
 * Class implementing the Farfalle construction
 */
class Farfalle
{
	protected:
		BaseIterableTransformation &p_b;
		BaseIterableTransformation &p_c;
		BaseIterableTransformation &p_d;
		BaseIterableTransformation &p_e;
		BaseRollingFunction        &roll_c;
		BaseRollingFunction        &roll_e;

	public:
		Farfalle(BaseIterableTransformation &p_b, BaseIterableTransformation &p_c, BaseIterableTransformation &p_d, BaseIterableTransformation &p_e, BaseRollingFunction &roll_c, BaseRollingFunction &roll_e);
		BitString     operator()(const BitString &K, const BitStrings &Mseq, unsigned int n, unsigned int q = 0) const;
		unsigned int  width() const;
};

/**
 * Class implementing Farfalle-SAE
 */
class FarfalleSAE
{
	private:
		Farfalle           F;
		const unsigned int t;
		const unsigned int l;
		BitString          K;
		BitStrings         history;
		unsigned int       offset;

	public:
		FarfalleSAE(const Farfalle &F, unsigned int t, unsigned int l, const BitString &K, const BitString &N, BitString &T, bool sender);
		pair<BitString, BitString>  wrap(const BitString &A, const BitString &P);
		BitString                   unwrap(const BitString &A, const BitString &C, const BitString &T);
};

/**
 * Class implementing Farfalle-SIV
 */
class FarfalleSIV
{
	private:
		Farfalle           F;
		const unsigned int t;

	public:
		FarfalleSIV(const Farfalle &F, unsigned int t);
		pair<BitString, BitString>  wrap(const BitString &K, const BitString &A, const BitString &P) const;
		BitString                   unwrap(const BitString &K, const BitString &A, const BitString &C, const BitString &T) const;
};

/**
 * Class implementing Farfalle-WBC
 */
class FarfalleWBC
{
	protected:
		Farfalle           H;
		Farfalle           G;
		const unsigned int l;

		unsigned int split(unsigned int n) const;

	public:
		FarfalleWBC(const Farfalle &H, const Farfalle &G, unsigned int l);
		BitString  encipher(const BitString &K, const BitString &W, const BitString &P) const;
		BitString  decipher(const BitString &K, const BitString &W, const BitString &C) const;
};

/**
 * Class implementing Farfalle-WBC-AE
 */
class FarfalleWBCAE : public FarfalleWBC
{
	private:
		const unsigned int t;

	public:
		FarfalleWBCAE(const Farfalle &H, const Farfalle &G, unsigned int t, unsigned int l);
		BitString  wrap(const BitString &K, const BitString &A, const BitString &P) const;
		BitString  unwrap(const BitString &K, const BitString &A, const BitString &C) const;
};

#endif
