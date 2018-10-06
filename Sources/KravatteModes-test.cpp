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

#include <exception>

#include "Kravatte.h"
#include "Keccak.h"

#define OUTPUT
#define VERBOSE_SIV
#define VERBOSE_SAE
#define VERBOSE_WBC
#define VERBOSE_WBC_AE

typedef unsigned char           BitSequence;
typedef size_t                  BitLength;

#define SnP_width               1600
#define SnP_widthInBytes        200
#define dataByteSize            (16*SnP_widthInBytes)
#define ADByteSize              (16*SnP_widthInBytes)
#define keyByteSize             (1*SnP_widthInBytes)
#define nonceByteSize           (2*SnP_widthInBytes)
#define WByteSize               (2*SnP_widthInBytes)
#define dataBitSize             (dataByteSize*8)
#define ADBitSize               (ADByteSize*8)
#define keyBitSize              (keyByteSize*8)
#define nonceBitSize            (nonceByteSize*8)
#define WBitSize                (WByteSize*8)
#define tagLenSANSE             32
#define tagLenSANE              16
#define expansionLenWBCAE       16

#define checksumByteSize        16

#if (defined(OUTPUT) || defined(VERBOSE) || !defined(EMBEDDED))
#include <stdio.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <time.h>

static void Farfalle_assert(int condition)
{
    if (!condition)
    {
        throw std::runtime_error("Invalid");
    }
}

static void randomize( unsigned char* data, unsigned int length)
{
    srand((unsigned int)time(0));
    while (length--)
    {
        *data++ = rand();
    }
}

static void generateSimpleRawMaterial(unsigned char* data, unsigned int length, unsigned char seed1, unsigned int seed2)
{
    unsigned int i;

    for(i=0; i<length; i++) {
        unsigned char iRolled;
        unsigned char byte;
        seed2 = seed2 % 8;
        iRolled = ((unsigned char)i << seed2) | ((unsigned char)i >> (8-seed2));
        byte = seed1 + 161*length - iRolled + i;
        data[i] = byte;
    }
}

/* ------------------------------------------------------------------------- */

static void performTestKravatteSANSE_OneInput(BitLength keyLen, BitLength dataLen, BitLength ADLen, Keccak &rSpongeChecksum)
{
    BitSequence input[dataByteSize];
    BitSequence inputPrime[dataByteSize];
    BitSequence output[dataByteSize];
    BitSequence AD[ADByteSize];
    BitSequence key[keyByteSize];
    unsigned char tag[tagLenSANSE];
    unsigned int seed;
    unsigned int session;

    randomize(key, keyByteSize);
    randomize(input, dataByteSize);
    randomize(inputPrime, dataByteSize);
    randomize(output, dataByteSize);
    randomize(AD, ADByteSize);
    randomize(tag, tagLenSANSE);

    seed = keyLen + dataLen + ADLen;
    seed ^= seed >> 3;
    generateSimpleRawMaterial(key, (keyLen + 7) / 8, 0x4321 - seed, 0x89 + seed);
    if (keyLen & 7)
        key[(keyLen + 7) / 8 - 1] &= (1 << (keyLen & 7)) - 1;
    generateSimpleRawMaterial(input, (dataLen + 7) / 8, 0x6523 - seed, 0x43 + seed);
    if (dataLen & 7)
        input[(dataLen + 7) / 8 - 1] &= (1 << (dataLen & 7)) - 1;
    generateSimpleRawMaterial(AD, (ADLen + 7) / 8, 0x1A29 - seed, 0xC3 + seed);
    if (ADLen & 7)
        AD[(ADLen + 7) / 8 - 1] &= (1 << (ADLen & 7)) - 1;

    #ifdef VERBOSE_SANSE
    printf( "keyLen %5u, dataLen %5u, ADLen %5u (in bits)\n", (unsigned int)keyLen, (unsigned int)dataLen, (unsigned int)ADLen);
    #endif

	KravatteSANSE xpEnc(BitString(key, keyLen));
	KravatteSANSE xpDec(BitString(key, keyLen));

    #ifdef VERBOSE_SANSE
    {
        unsigned int i;
        BitLength len;

        printf("Key of %d bits:", (int)keyLen);
        len = (keyLen + 7) / 8;
        for(i=0; (i<len) && (i<16); i++)
            printf(" %02x", (int)key[i]);
        if (len > 16)
            printf(" ...");
        printf("\n");

        printf("Input of %d bits:", (int)dataLen);
        len = (dataLen + 7) /8;
        for(i=0; (i<len) && (i<16); i++)
            printf(" %02x", (int)input[i]);
        if (len > 16)
            printf(" ...");
        printf("\n");

        printf("AD of %d bits:", (int)ADLen);
        len = (ADLen + 7) / 8;
        for(i=0; (i<len) && (i<16); i++)
            printf(" %02x", (int)AD[i]);
        if (len > 16)
            printf(" ...");
        printf("\n\n");
        fflush(stdout);
    }
    #endif

	for (session = 3; session != 0; --session) {
		const std::pair<BitString, BitString> ct = xpEnc.wrap(BitString(AD, ADLen), BitString(input, dataLen));
		if (ct.first.size() != 0) std::copy(ct.first.array(), ct.first.array() + (ct.first.size() + 7) / 8, output);
		if (ct.second.size() != 0) std::copy(ct.second.array(), ct.second.array() + (ct.second.size() + 7) / 8, tag);

		const BitString p_prime = xpDec.unwrap(BitString(AD, ADLen), ct.first, ct.second);
		if (p_prime.size() != 0) std::copy(p_prime.array(), p_prime.array() + (p_prime.size() + 7) / 8, inputPrime);

        Farfalle_assert(!memcmp(input,inputPrime,(dataLen + 7) / 8));
		rSpongeChecksum.absorb(output, 8 * ((dataLen + 7) / 8));
		rSpongeChecksum.absorb(tag, 8 * tagLenSANSE);
        #ifdef VERBOSE_SANSE
        {
            unsigned int i;
            unsigned int len;

            printf("Output of %d bits:", (int)dataLen);
            len = (dataLen + 7) / 8;
            for(i=0; (i<len) && (i<8); i++)
                printf(" %02x", (int)output[i]);
            if (len > 16)
                printf(" ...");
            if (i < (len - 8))
                i = len - 8;
            for( /* empty */; i<len; i++)
                printf(" %02x", (int)output[i]);
            printf("\n");

            printf("Tag of %d bytes:", (int)tagLenSANSE);
            for(i=0; i<tagLenSANSE; i++)
                printf(" %02x", (int)tag[i]);
            printf("\n");
            fflush(stdout);
            if (session == 1)
                printf("\n");
        }
        #endif
    }
}

static void performTestKravatteSANSE(unsigned char *checksum)
{
    BitLength dataLen, ADLen, keyLen;

    /* Accumulated test vector */
	Keccak spongeChecksum(SnP_width, 0);

    #ifdef OUTPUT
    printf("k ");
    #endif
    dataLen = 128*8;
    ADLen = 64*8;
    for(keyLen=0; keyLen<keyBitSize; keyLen = (keyLen < 2*SnP_width) ? (keyLen+1) : (keyLen+8)) {
        performTestKravatteSANSE_OneInput(keyLen, dataLen, ADLen, spongeChecksum);
    }
    
    #ifdef OUTPUT
    printf("d ");
    #endif
    ADLen = 64*8;
    keyLen = 16*8;
    for(dataLen=0; dataLen<=dataBitSize; dataLen = (dataLen < 2*SnP_width) ? (dataLen+1) : (dataLen+8)) {
        performTestKravatteSANSE_OneInput(keyLen, dataLen, ADLen, spongeChecksum);
    }
    
    #ifdef OUTPUT
    printf("a ");
    #endif
    dataLen = 128*8;
    keyLen = 16*8;
    for(ADLen=0; ADLen<=ADBitSize; ADLen = (ADLen < 2*SnP_width) ? (ADLen+1) : (ADLen+8)) {
        performTestKravatteSANSE_OneInput(keyLen, dataLen, ADLen, spongeChecksum);
    }
    
	spongeChecksum.squeeze(checksum, 8 * checksumByteSize);

    #ifdef VERBOSE_SANSE
    {
        unsigned int i;
        printf("Kravatte-SANSE\n" );
        printf("Checksum: ");
        for(i=0; i<checksumByteSize; i++)
            printf("\\x%02x", (int)checksum[i]);
        printf("\n\n");
    }
    #endif
}

void selfTestKravatteSANSE(const char *expected)
{
    unsigned char checksum[checksumByteSize];

    #if defined(OUTPUT)
    printf("Testing Kravatte-SANSE ");
    fflush(stdout);
    #endif
    performTestKravatteSANSE(checksum);
    Farfalle_assert(memcmp(expected, checksum, checksumByteSize) == 0);
    #if defined(OUTPUT)
    printf(" - OK.\n");
    #endif
}

#ifdef OUTPUT
void writeTestKravatteSANSE_One(FILE *f)
{
    unsigned char checksum[checksumByteSize];
    unsigned int offset;

    printf("Writing Kravatte-SANSE ");
    performTestKravatteSANSE(checksum);
    fprintf(f, "    selfTestKravatteSANSE(\"");
    for(offset=0; offset<checksumByteSize; offset++)
        fprintf(f, "\\x%02x", checksum[offset]);
    fprintf(f, "\");\n");
    printf("\n");
}

void writeTestKravatteSANSE(const char *filename)
{
    FILE *f = fopen(filename, "w");
    Farfalle_assert(f != NULL);
    writeTestKravatteSANSE_One(f);
    fclose(f);
}
#endif

/* ------------------------------------------------------------------------- */

static void performTestKravatteSANE_OneInput(BitLength keyLen, BitLength nonceLen, BitLength dataLen, BitLength ADLen, Keccak &rSpongeChecksum)
{
    BitSequence input[dataByteSize];
    BitSequence inputPrime[dataByteSize];
    BitSequence output[dataByteSize];
    BitSequence AD[ADByteSize];
    BitSequence key[keyByteSize];
    BitSequence nonce[nonceByteSize];
    unsigned char tag[tagLenSANE];
    unsigned char tagInit[tagLenSANE];
    unsigned int seed;
    unsigned int session;

    randomize(key, keyByteSize);
    randomize(nonce, nonceByteSize);
    randomize(input, dataByteSize);
    randomize(inputPrime, dataByteSize);
    randomize(output, dataByteSize);
    randomize(AD, ADByteSize);
    randomize(tag, tagLenSANE);

    seed = keyLen + nonceLen + dataLen + ADLen;
    seed ^= seed >> 3;
    generateSimpleRawMaterial(key, (keyLen + 7) / 8, 0x4371 - seed, 0x59 + seed);
    if (keyLen & 7)
        key[keyLen / 8] &= (1 << (keyLen & 7)) - 1;
    generateSimpleRawMaterial(nonce, (nonceLen + 7) / 8, 0x1327 - seed, 0x84 + seed);
    if (nonceLen & 7)
        nonce[nonceLen / 8] &= (1 << (nonceLen & 7)) - 1;
    generateSimpleRawMaterial(input, (dataLen + 7) / 8, 0x4861 - seed, 0xb1 + seed);
    if (dataLen & 7)
        input[dataLen / 8] &= (1 << (dataLen & 7)) - 1;
    generateSimpleRawMaterial(AD, (ADLen + 7) / 8, 0x243B - seed, 0x17 + seed);
    if (ADLen & 7)
        AD[ADLen / 8] &= (1 << (ADLen & 7)) - 1;

    #ifdef VERBOSE_SANE
    printf( "keyLen %5u, nonceLen %5u, dataLen %5u, ADLen %5u (in bits)\n", (unsigned int)keyLen, (unsigned int)nonceLen, (unsigned int)dataLen, (unsigned int)ADLen);
    #endif

	BitString bits_tagInit;
	KravatteSANE xpEnc(BitString(key, keyLen), BitString(nonce, nonceLen), bits_tagInit, true);
	if (bits_tagInit.size() != 0) std::copy(bits_tagInit.array(), bits_tagInit.array() + (bits_tagInit.size() + 7) / 8, tagInit);

	KravatteSANE xpDec(BitString(key, keyLen), BitString(nonce, nonceLen), bits_tagInit, false);

	rSpongeChecksum.absorb(tagInit, 8 * tagLenSANE);

    #ifdef VERBOSE_SANE
    {
        unsigned int i;
        unsigned int len;

        printf("Key of %d bits:", (int)keyLen);
        len = (keyLen + 7) / 8;
        for(i=0; (i<len) && (i<16); i++)
            printf(" %02x", (int)key[i]);
        if (len > 16)
            printf(" ...");
        printf("\n");

        printf("Nonce of %d bits:", (int)nonceLen);
        len = (nonceLen + 7) / 8;
        for(i=0; (i<len) && (i<16); i++)
            printf(" %02x", (int)nonce[i]);
        if (len > 16)
            printf(" ...");
        printf("\n");

        printf("Input of %d bits:", (int)dataLen);
        len = (dataLen + 7) / 8;
        for(i=0; (i<len) && (i<16); i++)
            printf(" %02x", (int)input[i]);
        if (len > 16)
            printf(" ...");
        printf("\n");

        printf("AD of %d bits:", (int)ADLen);
        len = (ADLen + 7) / 8;
        for(i=0; (i<len) && (i<16); i++)
            printf(" %02x", (int)AD[i]);
        if (len > 16)
            printf(" ...");
        printf("\n");
    }
    #endif

    for (session = 3; session != 0; --session) {
		const std::pair<BitString, BitString> ct = xpEnc.wrap(BitString(AD, ADLen), BitString(input, dataLen));
		if (ct.first.size() != 0) std::copy(ct.first.array(), ct.first.array() + (ct.first.size() + 7) / 8, output);
		if (ct.second.size() != 0) std::copy(ct.second.array(), ct.second.array() + (ct.second.size() + 7) / 8, tag);

		const BitString p_prime = xpDec.unwrap(BitString(AD, ADLen), ct.first, ct.second);
		if (p_prime.size() != 0) std::copy(p_prime.array(), p_prime.array() + (p_prime.size() + 7) / 8, inputPrime);

        Farfalle_assert(!memcmp(input,inputPrime,(dataLen + 7) / 8));
		rSpongeChecksum.absorb(output, 8 * ((dataLen + 7) / 8));
		rSpongeChecksum.absorb(tag, 8 * tagLenSANE);
        #ifdef VERBOSE_SANE
        {
            unsigned int i;
            unsigned int len;

            printf("Output of %d bits:", (int)dataLen);
            len = (dataLen + 7) / 8;
            for(i=0; (i<len) && (i<8); i++)
                printf(" %02x", (int)output[i]);
            if (len > 16)
                printf(" ...");
            if (i < (len - 8))
                i = len - 8;
            for( /* empty */; i<len; i++)
                printf(" %02x", (int)output[i]);
            printf("\n");

            printf("Tag of %d bytes:", (int)tagLenSANE);
            for(i=0; i<tagLenSANE; i++)
                printf(" %02x", (int)tag[i]);
            printf("\n");
            fflush(stdout);
            if (session == 1)
                printf("\n");
        }
        #endif
    }

}


static void performTestKravatteSANE(unsigned char *checksum)
{
    BitLength dataLen, ADLen, keyLen, nonceLen;

    /* Accumulated test vector */
	Keccak spongeChecksum(SnP_width, 0);

    #ifdef OUTPUT
    printf("k ");
    #endif
    dataLen = 128*8;
    ADLen = 64*8;
    nonceLen = 24*8;
    for(keyLen=0; keyLen<keyBitSize; keyLen = (keyLen < 2*SnP_width) ? (keyLen+1) : (keyLen+8)) {
        performTestKravatteSANE_OneInput(keyLen, nonceLen, dataLen, ADLen, spongeChecksum);
    }
    
    #ifdef OUTPUT
    printf("n ");
    #endif
    dataLen = 128*8;
    ADLen = 64*8;
    keyLen = 16*8;
    for(nonceLen=0; nonceLen<=nonceBitSize; nonceLen = (nonceLen < 2*SnP_width) ? (nonceLen+1) : (nonceLen+8)) {
        performTestKravatteSANE_OneInput(keyLen, nonceLen, dataLen, ADLen, spongeChecksum);
    }
    
    #ifdef OUTPUT
    printf("d ");
    #endif
    ADLen = 64*8;
    keyLen = 16*8;
    nonceLen = 24*8;
    for(dataLen=0; dataLen<=dataBitSize; dataLen = (dataLen < 2*SnP_width) ? (dataLen+1) : (dataLen+8)) {
        performTestKravatteSANE_OneInput(keyLen, nonceLen, dataLen, ADLen, spongeChecksum);
    }
    
    #ifdef OUTPUT
    printf("a ");
    #endif
    dataLen = 128*8;
    keyLen = 16*8;
    nonceLen = 24*8;
    for(ADLen=0; ADLen<=ADBitSize; ADLen = (ADLen < 2*SnP_width) ? (ADLen+1) : (ADLen+8)) {
        performTestKravatteSANE_OneInput(keyLen, nonceLen, dataLen, ADLen, spongeChecksum);
    }
    
	spongeChecksum.squeeze(checksum, 8 * checksumByteSize);

    #ifdef VERBOSE_SANE
    {
        unsigned int i;
        printf("Kravatte-SANE\n" );
        printf("Checksum: ");
        for(i=0; i<checksumByteSize; i++)
            printf("\\x%02x", (int)checksum[i]);
        printf("\n\n");
    }
    #endif
}

void selfTestKravatteSANE(const char *expected)
{
    unsigned char checksum[checksumByteSize];

    #if defined(OUTPUT)
    printf("Testing Kravatte-SANE ");
    fflush(stdout);
    #endif
    performTestKravatteSANE(checksum);
    Farfalle_assert(memcmp(expected, checksum, checksumByteSize) == 0);
    #if defined(OUTPUT)
    printf(" - OK.\n");
    #endif
}

#ifdef OUTPUT
void writeTestKravatteSANE_One(FILE *f)
{
    unsigned char checksum[checksumByteSize];
    unsigned int offset;

    printf("Writing Kravatte-SANE ");
    performTestKravatteSANE(checksum);
    fprintf(f, "    selfTestKravatteSANE(\"");
    for(offset=0; offset<checksumByteSize; offset++)
        fprintf(f, "\\x%02x", checksum[offset]);
    fprintf(f, "\");\n");
    printf("\n");
}

void writeTestKravatteSANE(const char *filename)
{
    FILE *f = fopen(filename, "w");
    Farfalle_assert(f != NULL);
    writeTestKravatteSANE_One(f);
    fclose(f);
}
#endif

/* ------------------------------------------------------------------------- */

static void performTestKravatte_WBC_OneInput(BitLength keyLen, BitLength dataLen, BitLength WLen, Keccak &rSpongeChecksum)
{
    BitSequence input[dataByteSize];
    BitSequence inputPrime[dataByteSize];
    BitSequence output[dataByteSize];
    BitSequence key[keyByteSize];
    BitSequence W[WByteSize];
    unsigned int seed;

    randomize(key, keyByteSize);
    randomize(W, WByteSize);
    randomize(input, dataByteSize);
    randomize(inputPrime, dataByteSize);
    randomize(output, dataByteSize);

    seed = keyLen + WLen + dataLen;
    seed ^= seed >> 3;
    generateSimpleRawMaterial(key, (keyLen + 7) / 8, 0x43C1 - seed, 0xB9 + seed);
    if (keyLen & 7)
        key[keyLen / 8] &= (1 << (keyLen & 7)) - 1;
    generateSimpleRawMaterial(W, (WLen + 7) / 8, 0x1727 - seed, 0x34 + seed);
    if (WLen & 7)
        W[WLen / 8] &= (1 << (WLen & 7)) - 1;
    generateSimpleRawMaterial(input, (dataLen + 7) / 8, 0x4165 - seed, 0xA9 + seed);
    if (dataLen & 7)
        input[dataLen / 8] &= (1 << (dataLen & 7)) - 1;

    #ifdef VERBOSE_WBC
    printf( "keyLen %5u, WLen %5u, dataLen %5u (in bits)\n", (unsigned int)keyLen, (unsigned int)WLen, (unsigned int)dataLen);
    #endif

	KravatteWBC kvw;

	// Remember the last byte since the copy() below will blow away the extra bits
	//BitSequence lastByte = (dataLen & 7) ? output[dataLen / 8] : 0;

	const BitString bits_output = kvw.encipher(BitString(key, keyLen), BitString(W, WLen), BitString(input, dataLen));
	if (bits_output.size() != 0) copy(bits_output.array(), bits_output.array() + (bits_output.size() + 7) / 8, output);

	// Restore the excess bits in the last byte
	//if (dataLen & 7)
	//{
	//	BitSequence mask = (1 << (dataLen & 7)) - 1;
	//	output[dataLen / 8] = (output[dataLen / 8] & mask) | (lastByte & ~mask);
	//}

	const BitString p_prime = kvw.decipher(BitString(key, keyLen), BitString(W, WLen), bits_output);
	if (p_prime.size() != 0) copy(p_prime.array(), p_prime.array() + (p_prime.size() + 7) / 8, inputPrime);

    Farfalle_assert(!memcmp(input,inputPrime,(dataLen + 7) / 8));

	rSpongeChecksum.absorb(output, 8 * ((dataLen + 7) / 8));

    #ifdef VERBOSE_WBC
    {
        unsigned int i;
        unsigned int dataByteLen;

        printf("Key of %d bits:", (int)keyLen);
        keyLen += 7;
        keyLen /= 8;
        for(i=0; (i<keyLen) && (i<16); i++)
            printf(" %02x", (int)key[i]);
        if (keyLen > 16)
            printf(" ...");
        printf("\n");

        printf("Tweak of %d bits:", (int)WLen);
        WLen += 7;
        WLen /= 8;
        for(i=0; (i<WLen) && (i<16); i++)
            printf(" %02x", (int)W[i]);
        if (WLen > 16)
            printf(" ...");
        printf("\n");

        printf("Input of %d bits:", (int)dataLen);
        dataByteLen = (dataLen + 7) / 8;
        for(i=0; (i<dataByteLen) && (i<16); i++)
            printf(" %02x", (int)input[i]);
        if (dataByteLen > 16)
            printf(" ...");
        printf("\n");

        printf("Output of %d bits:", (int)dataLen);
        for(i=0; (i<dataByteLen) && (i<8); i++)
            printf(" %02x", (int)output[i]);
        if (dataByteLen > 16)
            printf(" ...");
        if (i < (dataByteLen - 8))
            i = dataByteLen - 8;
        for( /* empty */; i<dataByteLen; i++)
            printf(" %02x", (int)output[i]);
        printf("\n\n");
        fflush(stdout);
    }
    #endif

}


static void performTestKravatte_WBC(unsigned char *checksum)
{
    BitLength dataLen, WLen, keyLen;

    /* Accumulated test vector */
	Keccak spongeChecksum(SnP_width, 0);

    #ifdef OUTPUT
    printf("k ");
    #endif
    dataLen = 128*8;
    WLen = 64*8;
    for(keyLen=0; keyLen<keyBitSize; keyLen = (keyLen < 2*SnP_width) ? (keyLen+1) : (keyLen+8)) {
        performTestKravatte_WBC_OneInput(keyLen, dataLen, WLen, spongeChecksum);
    }
    
    #ifdef OUTPUT
    printf("d ");
    #endif
    WLen = 64*8;
    keyLen = 16*8;
    for(dataLen=0; dataLen<=dataBitSize; dataLen = (dataLen < 2*SnP_width) ? (dataLen+1) : (dataLen+7)) {
        performTestKravatte_WBC_OneInput(keyLen, dataLen, WLen, spongeChecksum);
    }
    
    #ifdef OUTPUT
    printf("w ");
    #endif
    dataLen = 128*8;
    keyLen = 16*8;
    for(WLen=0; WLen<=WBitSize; WLen = (WLen < 2*SnP_width) ? (WLen+1) : (WLen+8)) {
        performTestKravatte_WBC_OneInput(keyLen, dataLen, WLen, spongeChecksum);
    }
    
	spongeChecksum.squeeze(checksum, 8 * checksumByteSize);

    #ifdef VERBOSE_WBC
    {
        unsigned int i;
        printf("Kravatte-WBC\n" );
        printf("Checksum: ");
        for(i=0; i<checksumByteSize; i++)
            printf("\\x%02x", (int)checksum[i]);
        printf("\n\n");
    }
    #endif
}

void selfTestKravatte_WBC(const char *expected)
{
    unsigned char checksum[checksumByteSize];

    printf("Testing Kravatte-WBC ");
    fflush(stdout);
    performTestKravatte_WBC(checksum);
    Farfalle_assert(memcmp(expected, checksum, checksumByteSize) == 0);
    printf(" - OK.\n");
}

#ifdef OUTPUT
void writeTestKravatte_WBC_One(FILE *f)
{
    unsigned char checksum[checksumByteSize];
    unsigned int offset;

    printf("Writing Kravatte-WBC ");
    performTestKravatte_WBC(checksum);
    fprintf(f, "    selfTestKravatte_WBC(\"");
    for(offset=0; offset<checksumByteSize; offset++)
        fprintf(f, "\\x%02x", checksum[offset]);
    fprintf(f, "\");\n");
    printf("\n");
}

void writeTestKravatte_WBC(const char *filename)
{
    FILE *f = fopen(filename, "w");
    Farfalle_assert(f != NULL);

    #if 0
    {
        BitLength n, nl, prevnl;

        prevnl = 0xFFFFFFFF;
        for (n = 0; n <= 8*64*1024; ++n )
        {
            nl = Kravatte_WBC_Split(n);
            if (nl != prevnl)
            {
                printf("n %6u, nl %6u, nr %6u", n, nl, n - nl);
                if (n >= 2*1536)
                    printf(", 2^x %6u", nl / 1536);
                printf("\n");
                prevnl = nl;
            }
        }
    }
    #endif

    writeTestKravatte_WBC_One(f);
    fclose(f);
}
#endif

/* ------------------------------------------------------------------------- */

static void performTestKravatte_WBC_AE_OneInput(BitLength keyLen, BitLength dataLen, BitLength ADLen, Keccak &rSpongeChecksum)
{
    BitLength outputLen = dataLen + 8 * expansionLenWBCAE;
    BitSequence input[dataByteSize];
    BitSequence inputPrime[dataByteSize];
    BitSequence output[dataByteSize];
    BitSequence key[keyByteSize];
    BitSequence AD[ADByteSize];
    unsigned int seed;

    randomize(key, keyByteSize);
    randomize(AD, ADByteSize);
    randomize(input, dataByteSize);
    randomize(inputPrime, dataByteSize);
    randomize(output, dataByteSize);

    seed = keyLen + ADLen + dataLen;
    seed ^= seed >> 3;
    generateSimpleRawMaterial(key, (keyLen + 7) / 8, 0x91FC - seed, 0x5A + seed);
    if (keyLen & 7)
        key[keyLen / 8] &= (1 << (keyLen & 7)) - 1;
    generateSimpleRawMaterial(AD, (ADLen + 7) / 8, 0x8181 - seed, 0x9B + seed);
    if (ADLen & 7)
        AD[ADLen / 8] &= (1 << (ADLen & 7)) - 1;
    generateSimpleRawMaterial(input, (dataLen + 7) / 8, 0x1BF0 - seed, 0xC6 + seed);
    if (dataLen & 7)
        input[dataLen / 8] &= (1 << (dataLen & 7)) - 1;

    #ifdef VERBOSE_WBC_AE
    printf( "keyLen %5u, ADLen %5u, dataLen %5u (in bits)\n", (unsigned int)keyLen, (unsigned int)ADLen, (unsigned int)dataLen);
    #endif

	KravatteWBCAE kvw;

	const BitString bits_output = kvw.wrap(BitString(key, keyLen), BitString(AD, ADLen), BitString(input, dataLen));
	if (bits_output.size() != 0) copy(bits_output.array(), bits_output.array() + (bits_output.size() + 7) / 8, output);

	const BitString p_prime = kvw.unwrap(BitString(key, keyLen), BitString(AD, ADLen), bits_output);
	if (p_prime.size() != 0) copy(p_prime.array(), p_prime.array() + (p_prime.size() + 7) / 8, inputPrime);

    Farfalle_assert(!memcmp(input,inputPrime,(dataLen + 7) / 8));

	rSpongeChecksum.absorb(output, 8 * ((outputLen + 7) / 8));

    #ifdef VERBOSE_WBC_AE
    {
        unsigned int i;
        unsigned int dataByteLen;
        unsigned int outputByteLen;

        printf("Key of %d bits:", (int)keyLen);
        keyLen += 7;
        keyLen /= 8;
        for(i=0; (i<keyLen) && (i<16); i++)
            printf(" %02x", (int)key[i]);
        if (keyLen > 16)
            printf(" ...");
        printf("\n");

        printf("AD of %d bits:", (int)ADLen);
        ADLen += 7;
        ADLen /= 8;
        for(i=0; (i<ADLen) && (i<16); i++)
            printf(" %02x", (int)AD[i]);
        if (ADLen > 16)
            printf(" ...");
        printf("\n");

        printf("Input of %d bits:", (int)dataLen);
        dataByteLen = (dataLen + 7) / 8;
        for(i=0; (i<dataByteLen) && (i<16); i++)
            printf(" %02x", (int)input[i]);
        if (dataByteLen > 16)
            printf(" ...");
        printf("\n");

        printf("Output of %d bits:", (int)outputLen);
        outputByteLen = (outputLen + 7) / 8;
        for(i=0; (i<outputByteLen) && (i<8); i++)
            printf(" %02x", (int)output[i]);
        if (outputByteLen > 16)
            printf(" ...");
        if (i < (outputByteLen - 8))
            i = outputByteLen - 8;
        for( /* empty */; i<outputByteLen; i++)
            printf(" %02x", (int)output[i]);
        printf("\n\n");
        fflush(stdout);
    }
    #endif

}


static void performTestKravatte_WBC_AE(unsigned char *checksum)
{
    BitLength dataLen, ADLen, keyLen;

    /* Accumulated test vector */
	Keccak spongeChecksum(SnP_width, 0);

    #ifdef OUTPUT
    printf("k ");
    #endif
    dataLen = 128*8;
    ADLen = 64*8;
    for(keyLen=0; keyLen<keyBitSize; keyLen = (keyLen < 2*SnP_width) ? (keyLen+1) : (keyLen+8)) {
        performTestKravatte_WBC_AE_OneInput(keyLen, dataLen, ADLen, spongeChecksum);
    }
    
    #ifdef OUTPUT
    printf("d ");
    #endif
    ADLen = 64*8;
    keyLen = 16*8;
	for(dataLen=0; dataLen<=dataBitSize-8*expansionLenWBCAE; dataLen = (dataLen < 2*SnP_width) ? (dataLen+1) : (dataLen+7)) {
        performTestKravatte_WBC_AE_OneInput(keyLen, dataLen, ADLen, spongeChecksum);
    }
    
    #ifdef OUTPUT
    printf("a ");
    #endif
    dataLen = 128*8;
    keyLen = 16*8;
    for(ADLen=0; ADLen<=ADBitSize; ADLen = (ADLen < 2*SnP_width) ? (ADLen+1) : (ADLen+8)) {
        performTestKravatte_WBC_AE_OneInput(keyLen, dataLen, ADLen, spongeChecksum);
    }
    
	spongeChecksum.squeeze(checksum, 8 * checksumByteSize);

    #ifdef VERBOSE_WBC_AE
    {
        unsigned int i;
        printf("Kravatte-WBC-AE\n" );
        printf("Checksum: ");
        for(i=0; i<checksumByteSize; i++)
            printf("\\x%02x", (int)checksum[i]);
        printf("\n\n");
    }
    #endif
}

void selfTestKravatte_WBC_AE(const char *expected)
{
    unsigned char checksum[checksumByteSize];

    printf("Testing Kravatte-WBC-AE ");
    fflush(stdout);
    performTestKravatte_WBC_AE(checksum);
    Farfalle_assert(memcmp(expected, checksum, checksumByteSize) == 0);
    printf(" - OK.\n");
}

#ifdef OUTPUT
void writeTestKravatte_WBC_AE_One(FILE *f)
{
    unsigned char checksum[checksumByteSize];
    unsigned int offset;

    printf("Writing Kravatte-WBC-AE ");
    performTestKravatte_WBC_AE(checksum);
    fprintf(f, "    selfTestKravatte_WBC_AE(\"");
    for(offset=0; offset<checksumByteSize; offset++)
        fprintf(f, "\\x%02x", checksum[offset]);
    fprintf(f, "\");\n");
    printf("\n");
}

void writeTestKravatte_WBC_AE(const char *filename)
{
    FILE *f = fopen(filename, "w");
    Farfalle_assert(f != NULL);

    #if 0
    {
        BitLength n, nl, prevnl;

        prevnl = 0xFFFFFFFF;
        for (n = 0; n <= 8*64*1024; ++n )
        {
            nl = Kravatte_WBC_AE_Split(n);
            if (nl != prevnl)
            {
                printf("n %6u, nl %6u, nr %6u", n, nl, n - nl);
                if (n >= 2*1536)
                    printf(", 2^x %6u", nl / 1536);
                printf("\n");
                prevnl = nl;
            }
        }
    }
    #endif

    writeTestKravatte_WBC_AE_One(f);
    fclose(f);
}
#endif

/* ------------------------------------------------------------------------- */
void testKravatteModes(void)
{
#ifndef KeccakP1600_excluded
#ifdef OUTPUT
//    printKravatteTestVectors();
    writeTestKravatteSANE("Kravatte_SANE.txt");
    writeTestKravatteSANSE("Kravatte_SANSE.txt");
    writeTestKravatte_WBC("Kravatte_WBC.txt");
    writeTestKravatte_WBC_AE("Kravatte_WBC_AE.txt");
#endif

    //selfTestKravatteSANE("\x99\x93\x6b\x93\x37\xcb\x88\xb1\xe0\x97\xf6\x79\x98\xc4\xa6\x75");
    //selfTestKravatteSANSE("\x25\x15\x37\xb4\xc7\x1f\x14\x37\x22\x3d\x59\xa2\x92\xf6\x6d\x82");
    //selfTestKravatte_WBC("\x42\xae\xe6\x1b\xcf\x8f\x06\x34\xc6\xb2\x11\x0a\xf4\xa0\xdd\xc1");
    //selfTestKravatte_WBC_AE("\xda\xa6\x2c\xab\xa7\xe4\xe1\xef\xb5\x4c\x62\x78\x26\xf8\x72\x27");
#endif
}
