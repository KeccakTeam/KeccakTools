/*
 * Implementation by the Keccak, Keyak and Ketje Teams, namely, Guido Bertoni,
 * Joan Daemen, Michaël Peeters, Gilles Van Assche and Ronny Van Keer, hereby
 * denoted as "the implementer".
 *
 * For more information, feedback or questions, please refer to our websites:
 * http://keccak.noekeon.org/
 * http://keyak.noekeon.org/
 * http://ketje.noekeon.org/
 *
 * To the extent possible under law, the implementer has waived all copyright
 * and related or neighboring rights to the source code in this file.
 * http://creativecommons.org/publicdomain/zero/1.0/
 */

#include <fstream>
#include "Ketjev2.h"
#include "Ketjev2-test.h"

static string generateSimpleRawMaterial(unsigned int length, UINT8 seed1, unsigned int seed2)
{
    string  s;

    for ( unsigned int i = 0; i < length; i++ ) {
        UINT8  iRolled = ((UINT8)i << seed2) | ((UINT8)i >> (8 - seed2));
        UINT8  byte    = seed1 + 161 * length - iRolled + i;
        s += char(byte);
    }
    return s;
}

static void displayByteString(ostream &fout, const string &synopsis, const string &data)
{
    fout << synopsis << ":";
    for ( string::const_iterator i = data.begin(); i != data.end(); ++i ) {
        fout << " ";
        fout.width(2);
        fout.fill('0');
        fout << hex << ((int)UINT8(*i));
    }
    fout << endl;
}

static void assert(bool condition, const string &synopsis)
{
    if ( !condition ) {
        throw Exception(synopsis);
    }
}

void testKetje(ostream &fout, Ketje ketje, const string &Texpected)
{
    string        Aglobal;
    unsigned int  width = ketje.getWidth();

    for ( int Klen = (width - 18) / 8; Klen >= 12; Klen -= (Klen > 47) ? 25 : ((Klen > 22) ? 5 : ((Klen > 14) ? 2 : 1))) {
        int  NlenMax = (width - 18) / 8 - Klen;
        for ( int Nlen = ((Klen == 16) ? 0 : NlenMax); Nlen <= NlenMax; Nlen += width <= 400 ? 1 : (width / 200)) {
            string  K = generateSimpleRawMaterial(Klen, 0x12 + Nlen * 8, width);
            string  N = generateSimpleRawMaterial(Nlen, 0x23 + Klen * 8, width);

            cout << "Ketje(" << dec << width << "), key length is " << dec << (Klen * 8) << " bits, nonce length is " << dec <<
            (Nlen * 8) << " bits" << endl;

            Ketje   ketje1(ketje);
            Ketje   ketje2(ketje);
            ketje1.initialize(K, N);
            ketje2.initialize(K, N);

            fout << "***" << endl;
            fout << "initialize with key of " << dec << (Klen * 8) << " bits, nonce of " << dec << (Nlen * 8) << " bits:" << endl;
            displayByteString(fout, "> K (key)",   K);
            displayByteString(fout, "> N (nonce)", N);
            fout << endl;

            for ( unsigned int Alen = 0; Alen <= 50; Alen += 1 + Alen / 3 + (Klen - 12) + Nlen / 32 ) {
                for ( unsigned int Blen = 0;
                      Blen <= 50;
                      Blen += Blen / 2 + 1 + Alen + ((Alen == 0) ? (Klen - 12) : (Nlen / 32 + Klen / 2))) {
                    for ( unsigned int ell = ((Klen == 16) ? 0 : 128); ell <= ((Klen == 16) ? 256 : 128); ell += 64 ) {
                        string  A = generateSimpleRawMaterial(Alen, 0x34 + Blen, 3);
                        string  B = generateSimpleRawMaterial(Blen, 0x45 + Alen, 4);
                        string  Bprime;
                        string  C;
                        string  T;

                        displayByteString(fout, "> A (associated data)", A);
                        displayByteString(fout, "> B (plaintext)",       B);

                        try {
                            C      = ketje1.wrap(A, B, ell, T);

                            displayByteString(fout, "< C (ciphertext)", C);
                            displayByteString(fout, "< T (tag)",        T);

                            Bprime = ketje2.unwrap(A, C, T);

                            assert(B == Bprime, "The plaintexts do not match.");
                        }
                        catch ( ... ) {
                            fout << "Klen: " << dec << Klen << " ";
                            fout << "Nlen: " << dec << Nlen << " ";
                            fout << "Alen: " << dec << Alen << " ";
                            fout << "Blen: " << dec << Blen << endl;
                            displayByteString(fout, "< Bprime (plaintext)", Bprime);
                            throw;
                        }

                        Aglobal += C;
                        Aglobal += T;

                        fout << endl;
                    }
                }
            }
        }
    }

    {
        Ketje   global(ketje);
        string  Tglobal;

        global.initialize(string(""), string(""));
        (void)global.wrap(Aglobal, string(""), 128, Tglobal);
        displayByteString(fout, "+++ Global tag ", Tglobal);
        assert(Tglobal == Texpected, "The global tag does not match.");
    }
}

int tryTestKetje(ostream &fout, Ketje ketje, const string &Texpected)
{
    try {
        testKetje(fout, ketje, Texpected);
    }
    catch ( Exception e ) {
        cout << "Exception: " << e.reason << endl;
        return 1;
    }
    return 0;
}

int testAllKetjev2Instances(void)
{
    int  errors = 0;
    {
        ofstream  fout("KetjeJr.txt");
        errors += tryTestKetje(fout, KetjeJr(), string("\xde\xd0\xa0\x53\xa5\x5c\xbb\x73\xba\x2b\x3a\xb9\x16\x76\x62\xc7"));
    }
    {
        ofstream  fout("KetjeSr.txt");
        errors += tryTestKetje(fout, KetjeSr(), string("\x32\x58\xf4\x15\x47\x45\x84\x64\x3b\x52\x43\x50\x3b\x1a\x25\x41"));
    }
    {
        ofstream  fout("KetjeMn.txt");
        errors += tryTestKetje(fout, KetjeMinor(), string("\x2d\xd3\x9e\xb2\xa3\xfd\x61\xf9\x54\x0c\xf5\x5d\xc8\x37\xb2\x14"));
    }
    {
        ofstream  fout("KetjeMj.txt");
        errors += tryTestKetje(fout, KetjeMajor(), string("\x17\x66\xc5\x9c\xa6\x0e\x3d\x28\xc6\xa8\x87\x72\xe1\x9d\xf4\x95"));
    }

    cout << __FUNCTION__ << ": " << errors << " error(s)." << endl;
    return errors;
}
