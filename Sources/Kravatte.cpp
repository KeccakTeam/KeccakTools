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

#include <algorithm>
#include "Kravatte.h"

#define ROL64(a, offset) ((((UINT64)a) << offset) ^ (((UINT64)a) >> (64 - offset)))

/* KravatteCompressionRollingFunction */
BitString KravatteCompressionRollingFunction::operator()(const BitString &k, unsigned int i) const
{
	vector<UINT8> kp(k.array(), k.array() + k.size() / 8);
	UINT64 *lanes = reinterpret_cast<UINT64 *>(kp.data());

	for (unsigned int j = 0; j < i; j++)
	{
		UINT64 x0 = lanes[5 * 4    ];
		UINT64 x1 = lanes[5 * 4 + 1];

		std::rotate(&lanes[5 * 4], &lanes[5 * 4 + 1], &lanes[5 * 5]);
		lanes[5 * 4 + 4] = (ROL64(x0, 7) ^ x1 ^ (x1 >> 3));
	}

	return BitString(kp);
}

/* KravatteExpansionRollingFunction */
BitString KravatteExpansionRollingFunction::operator()(const BitString &k, unsigned int i) const
{
	vector<UINT8> kp(k.array(), k.array() + k.size() / 8);
	UINT64 *lanes = reinterpret_cast<UINT64 *>(kp.data());

	for (unsigned int j = 0; j < i; j++)
	{
		UINT64 x0 = lanes[5 * 3    ];
		UINT64 x1 = lanes[5 * 3 + 1];
		UINT64 x2 = lanes[5 * 3 + 2];

		std::rotate(&lanes[5 * 3], &lanes[5 * 3 + 1], &lanes[5 * 5]);
		lanes[5 * 4 + 4] = (ROL64(x0, 7) ^ ROL64(x1, 18) ^ (x2 & (x1 >> 1)));
	}

	return BitString(kp);
}

/* Kravatte Achouffe instantiation parameters */
namespace KravatteParams
{
	IterableTransformation<KeccakP>    p_b(1600, 6);
	IterableTransformation<KeccakP>    p_c(1600, 6);
	IterableTransformation<KeccakP>    p_d(1600, 6);
	IterableTransformation<KeccakP>    p_e(1600, 6);
	IterableTransformation<Identity>   p_identity(1600);

	KravatteCompressionRollingFunction roll_c;
	KravatteExpansionRollingFunction   roll_e;

	unsigned int                       param_SAE_t = 128;
	unsigned int                       param_SAE_l = 8;

	unsigned int                       param_SIV_t = 256;

	unsigned int                       param_WBC_l = 8;

	unsigned int                       param_WBC_AE_t = 128;
	unsigned int                       param_WBC_AE_l = 8;
};

Farfalle make_Short_Kravatte()
{
	return Farfalle(KravatteParams::p_b, KravatteParams::p_c, KravatteParams::p_identity, KravatteParams::p_e, KravatteParams::roll_c, KravatteParams::roll_e);
}

Farfalle make_Kravatte()
{
	return Farfalle(KravatteParams::p_b, KravatteParams::p_c, KravatteParams::p_d, KravatteParams::p_e, KravatteParams::roll_c, KravatteParams::roll_e);
}

/* Kravatte */
Kravatte::Kravatte()
	: Farfalle(make_Kravatte())
{
}

/* Kravatte-SAE */
KravatteSAE::KravatteSAE(const BitString &K, const BitString &N, BitString &T, bool sender)
	: FarfalleSAE(make_Kravatte(), KravatteParams::param_SAE_t, KravatteParams::param_SAE_l, K, N, T, sender)
{
}

/* Kravatte-SIV */
KravatteSIV::KravatteSIV()
	: FarfalleSIV(make_Kravatte(), KravatteParams::param_SIV_t)
{
}

/* Kravatte-WBC */
KravatteWBC::KravatteWBC()
	: FarfalleWBC(make_Short_Kravatte(), make_Kravatte(), KravatteParams::param_WBC_l)
{
}

/* Kravatte-WBC-AE */
KravatteWBCAE::KravatteWBCAE()
	: FarfalleWBCAE(make_Short_Kravatte(), make_Kravatte(), KravatteParams::param_WBC_AE_t, KravatteParams::param_WBC_AE_l)
{
}
