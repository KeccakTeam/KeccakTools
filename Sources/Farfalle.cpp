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
#include <cmath>
#include <string>
#include <vector>
#include "Farfalle.h"

using namespace std;

static void _assert(bool condition, const string &synopsis, const char *fct)
{
	if (!condition)
	{
		throw Exception(string(fct) + "(): " + synopsis);
	}
}

#if defined(__GNUC__)
#define assert(cond, msg)  _assert(cond, msg, __PRETTY_FUNCTION__)
#else
#define assert(cond, msg)  _assert(cond, msg, __FUNCTION__)
#endif

/* IdentityRollingFunction */
BitString IdentityRollingFunction::operator()(const BitString &k, unsigned int i) const
{
    (void)i;
	return k;
}

/* Farfalle */
Farfalle::Farfalle(BaseIterableTransformation &p_b,
                   BaseIterableTransformation &p_c,
                   BaseIterableTransformation &p_d,
                   BaseIterableTransformation &p_e,
                   BaseRollingFunction        &roll_c,
                   BaseRollingFunction        &roll_e)
	: p_b(p_b), p_c(p_c), p_d(p_d), p_e(p_e), roll_c(roll_c), roll_e(roll_e)
{
	assert((p_b.width % 8) == 0 && (p_c.width % 8) == 0 && (p_d.width % 8) == 0 && (p_e.width % 8) == 0,
		"This implementation only supports permutation width that are multiple of 8."); // Limitation of Transformation class
}

BitString Farfalle::operator()(const BitString &K, const BitStrings &Mseq, unsigned int n, unsigned int q) const
{
	unsigned int b = width();
	if (!(K.size() <= b - 1)) throw Exception("Key length must be less than b bits");
	unsigned int m = Mseq.size();

	BitString Kp = K || BitString::pad10(b, K.size());
	BitString k = p_b(Kp);

	BitString    x = BitString::zeroes(b);
	unsigned int I = 0;

	for (unsigned int j = 0; j <= m - 1; j++)
	{
		unsigned int mu = (Mseq[j].size() + b) / b;
		BitString M = Mseq[j] || BitString::pad10(mu * b, Mseq[j].size());
		Blocks mblocks(M, b);

		for (unsigned int i = I; i <= I + mu - 1; i++)
		{
			x = x ^ p_c(mblocks[i - I] ^ roll_c(k, i));
		}

		I = I + mu + 1;
	}

	BitString kp = roll_c(k, I);
	BitString y = p_d(x);

	Blocks zblocks(width());

	for (unsigned int j = 0; width() * j < n + q; j++)
	{
		zblocks[j] = p_e(roll_e(y, j)) ^ kp;
	}

	BitString Z = BitString::substring(zblocks.bits(), q, n);
	return Z;
}

unsigned int Farfalle::width() const
{
	return p_b.width;
}

/* Farfalle-SAE */
FarfalleSAE::FarfalleSAE(const Farfalle  &F,
                         unsigned int     t,
                         unsigned int     l,
                         const BitString &K,
                         const BitString &N,
                         BitString &T,
                         bool sender)
	: F(F), t(t), l(l), K(K)
{
	offset = l * ((t + l - 1) / l);
	history = N;
	BitString Tp = F(K, history, t);

	if (sender)
	{
		T = Tp;
	}
	else if (!(Tp == T))
	{
		throw Exception("error!");
	}
}

pair<BitString, BitString> FarfalleSAE::wrap(const BitString &A, const BitString &P)
{
	BitString C = P ^ F(K, history, P.size(), offset);

	if (A.size() > 0 || P.size() == 0)
	{
		history = (A || 0) * history;
	}

	if (P.size() > 0)
	{
		history = (C || 1) * history;
	}

	BitString T = F(K, history, t);
	return make_pair(C, T);
}

BitString FarfalleSAE::unwrap(const BitString &A, const BitString &C, const BitString &T)
{
	BitString P = C ^ F(K, history, C.size(), offset);

	if (A.size() > 0 || C.size() == 0)
	{
		history = (A || 0) * history;
	}

	if (C.size() > 0)
	{
		history = (C || 1) * history;
	}

	BitString Tp = F(K, history, t);

	if (Tp == T)
	{
		return P;
	}
	else
	{
		throw Exception("error!");
	}
}

/* Farfalle-SIV */
FarfalleSIV::FarfalleSIV(const Farfalle  &F,
                         unsigned int     t)
	: F(F), t(t)
{
}

pair<BitString, BitString> FarfalleSIV::wrap(const BitString &K, const BitString &A, const BitString &P) const
{
	BitString T = F(K, P * A, t);
	BitString C = P ^ F(K, T * A, P.size());
	return make_pair(C, T);
}

BitString FarfalleSIV::unwrap(const BitString &K, const BitString &A, const BitString &C, const BitString &T) const
{
	BitString P  = C ^ F(K, T * A, C.size());
	BitString Tp = F(K, P * A, t);

	if (Tp == T)
	{
		return P;
	}
	else
	{
		throw Exception("error!");
	}
}

/* Farfalle-WBC */
FarfalleWBC::FarfalleWBC(const Farfalle  &H,
                         const Farfalle  &G,
                         unsigned int     l)
	: H(H), G(G), l(l)
{
}

unsigned int FarfalleWBC::split(unsigned int n) const
{
	unsigned int b = H.width();
	unsigned int n_L;

	if (n <= 2 * b - (l + 2))
	{
		n_L = l * ((n + l) / (2 * l));
	}
	else
	{
		unsigned int q = (n + l + 1 + b) / b;
		unsigned int tx = 1;
		while ((tx << 1) < q) tx <<= 1;
		n_L = (q - tx) * b - l;
	}

	return n_L;
}

BitString FarfalleWBC::encipher(const BitString &K, const BitString &W, const BitString &P) const
{
	unsigned int b = H.width();

	unsigned int n_L = split(P.size());
	unsigned int n_R = P.size() - n_L;
	BitString L = BitString::substring(P, 0, n_L);
	BitString R = BitString::substring(P, n_L, n_R);

	BitString Hval = H(K, (L || 0), min(b, R.size()));
	R = R ^ (Hval || BitString::zeroes(R.size() - Hval.size()));
	L = L ^ G(K, (R || 1) * W, L.size());
	R = R ^ G(K, (L || 0) * W, R.size());
	Hval = H(K, (R || 1), min(b, L.size()));
	L = L ^ (Hval || BitString::zeroes(L.size() - Hval.size()));
	return L || R;
}

BitString FarfalleWBC::decipher(const BitString &K, const BitString &W, const BitString &C) const
{
	unsigned int b = H.width();

	unsigned int n_L = split(C.size());
	unsigned int n_R = C.size() - n_L;
	BitString L = BitString::substring(C, 0, n_L);
	BitString R = BitString::substring(C, n_L, n_R);

	BitString Hval = H(K, (R || 1), min(b, L.size()));
	L = L ^ (Hval || BitString::zeroes(L.size() - Hval.size()));
	R = R ^ G(K, (L || 0) * W, R.size());
	L = L ^ G(K, (R || 1) * W, L.size());
	Hval = H(K, (L || 0), min(b, R.size()));
	R = R ^ (Hval || BitString::zeroes(R.size() - Hval.size()));
	return L || R;
}

/* Farfalle-WBC-AE */
FarfalleWBCAE::FarfalleWBCAE(const Farfalle  &H,
                             const Farfalle  &G,
                             unsigned int     t,
                             unsigned int     l)
	: FarfalleWBC(H, G, l), t(t)
{
}

BitString FarfalleWBCAE::wrap(const BitString &K, const BitString &A, const BitString &P) const
{
	BitString Pp = P || BitString::zeroes(t);
	return encipher(K, A, Pp);
}

BitString FarfalleWBCAE::unwrap(const BitString &K, const BitString &A, const BitString &C) const
{
	unsigned int b = H.width();

	unsigned int n_L = split(C.size());
	unsigned int n_R = C.size() - n_L;
	BitString L = BitString::substring(C, 0, n_L);
	BitString R = BitString::substring(C, n_L, n_R);

	BitString Hval = H(K, (R || 1), min(b, L.size()));
	L = L ^ (Hval || BitString::zeroes(L.size() - Hval.size()));
	R = R ^ G(K, (L || 0) * A, R.size());

	if (R.size() >= b + t)
	{
		if (!(BitString::substring(R, R.size() - t, t) == BitString::zeroes(t))) throw Exception("error!");
		L = L ^ G(K, (R || 1) * A, L.size());
		Hval = H(K, (L || 0), b);
		R = R ^ (Hval || BitString::zeroes(R.size() - Hval.size()));
	}
	else
	{
		L = L ^ G(K, (R || 1) * A, L.size());
		Hval = H(K, (L || 0), min(b, R.size()));
		R = R ^ (Hval || BitString::zeroes(R.size() - Hval.size()));
		if (!(BitString::substring(L || R, C.size() - t, t) == BitString::zeroes(t))) throw Exception("error!");
	}

	BitString Pp = L || R;
	return Pp.truncate(C.size() - t);
}
