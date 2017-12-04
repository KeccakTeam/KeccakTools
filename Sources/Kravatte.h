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

#ifndef _KRAVATTE_H_
#define _KRAVATTE_H_

#include <iostream>
#include <memory>
#include "padding.h"
#include "transformations.h"
#include "types.h"
#include "bitstring.h"
#include "Keccak-f.h"
#include "Farfalle.h"

using namespace std;

class KravatteCompressionRollingFunction : public BaseRollingFunction
{
	public:
		BitString operator()(const BitString &k, unsigned int i) const;
};

class KravatteExpansionRollingFunction : public BaseRollingFunction
{
	public:
		BitString operator()(const BitString &k, unsigned int i) const;
};

class Kravatte : public Farfalle
{
	public:
		Kravatte();
};

class KravatteSAE : public FarfalleSAE
{
	public:
		KravatteSAE(const BitString &K, const BitString &N, BitString &T, bool sender);
};

class KravatteSIV : public FarfalleSIV
{
	public:
		KravatteSIV();
};

class KravatteWBC : public FarfalleWBC
{
	public:
		KravatteWBC();
};

class KravatteWBCAE : public FarfalleWBCAE
{
	public:
		KravatteWBCAE();
};

#endif
