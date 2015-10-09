/*
Implementation by the Keccak, Keyak and Ketje Teams, namely, Guido Bertoni,
Joan Daemen, Michaël Peeters, Gilles Van Assche and Ronny Van Keer, hereby
denoted as "the implementer".

For more information, feedback or questions, please refer to our websites:
http://keccak.noekeon.org/
http://keyak.noekeon.org/
http://ketje.noekeon.org/

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/
*/

#include <fstream>
#include "Keyakv2.h"

string generateSimpleRawMaterial(unsigned int length, UINT8 seed1, unsigned int seed2)
{
    string s;
    for(unsigned int i=0; i<length; i++) {
        UINT8 iRolled = ((UINT8)i << seed2) | ((UINT8)i >> (8-seed2));
        UINT8 byte = seed1 + 161*length - iRolled + i;
        s += char(byte);
    }
    return s;
}

void displayByteString(ostream& fout, const string& synopsis, const string& data)
{
    fout << synopsis << ":";
    for(string::const_iterator i=data.begin(); i!=data.end(); ++i) {
        fout << " ";
        fout.width(2); fout.fill('0'); fout << hex << ((int)UINT8(*i));
    }
    fout << endl;
}

void assert(bool condition, const string& synopsis)
{
    if (!condition)
        throw Exception(synopsis);
}

void testKeyakStartEngine(Keyak& global, Keyak& wrap, Keyak& unwrap, ofstream& fout, const string& K, const string& N, bool forgetFlag, bool tagFlag)
{
    stringstream T;
    bool rv;

    fout << "*** " << wrap << "\n";
    fout << "StartEngine(K, N, tagFlag=" << (tagFlag?"true":"false") << ", T, unwrapFlag=false, forgetFlag=" <<  (forgetFlag?"true":"false") << "), with:\n";
    displayByteString(fout, "> K", K);
    displayByteString(fout, "> N", N);
    rv = wrap.StartEngine(K, N, tagFlag, T, false, forgetFlag);
    assert(rv, "wrap.StartEngine() did not return true.");
    if (tagFlag) {
        displayByteString(fout, "< T (tag)", T.str());
        stringstream empty, dummy, TT(T.str());
        global.Wrap(empty, dummy, TT, dummy, false, false);
    }

    T.seekg(0, ios_base::beg);
    rv = unwrap.StartEngine(K, N, tagFlag, T, true, forgetFlag);
    assert(rv, "unwrap.StartEngine() did not return true.");
    fout << endl;
}

void testKeyakWrapUnwrap(Keyak& global, Keyak& wrap, Keyak& unwrap, ofstream& fout, const string& Acontent, const string& Pcontent, bool forgetFlag)
{
    stringstream metadata(Acontent);
    stringstream plaintext(Pcontent);
    bool rv;

    fout << "Wrap(I, O, A, T, unwrapFlag=false, forgetFlag=" <<  (forgetFlag?"true":"false") << "), with:\n";
    displayByteString(fout, "> A (metadata)", metadata.str());
    displayByteString(fout, "> I (plaintext)", plaintext.str());
    stringstream ciphertext;
    stringstream plaintextPrime;
    stringstream tag;
    rv = wrap.Wrap(plaintext, ciphertext, metadata, tag, false, forgetFlag!=0);
    assert(rv, "wrap.Wrap() did not return true.");

    displayByteString(fout, "< O (ciphertext)", ciphertext.str());
    displayByteString(fout, "< T (tag)", tag.str());
    fout << endl;
    {
        stringstream empty, dummy, O(ciphertext.str()), T(tag.str());
        global.Wrap(empty, dummy, O, dummy, false, false);
        global.Wrap(empty, dummy, T, dummy, false, false);
    }

    ciphertext.seekg(0, ios_base::beg);
    metadata.clear(); metadata.seekg(0, ios_base::beg);
    tag.seekg(0, ios_base::beg);
    rv = unwrap.Wrap(ciphertext, plaintextPrime, metadata, tag, true, forgetFlag!=0);
    assert(rv, "unwrap.Wrap() did not return true.");
    assert(plaintext.str() == plaintextPrime.str(), "The plaintexts do not match.");
}

int testKeyak(ofstream& fout, unsigned int b, unsigned int nr, unsigned int Pi, unsigned int c, unsigned int tau, const string& expectedGlobalTag)
{
    Keyak global(b, nr, Pi, c, tau);
    cout << global << endl;
    {
        stringstream dummy;
        global.StartEngine(string(""), string(""), false, dummy, false, false);
    }

    unsigned int Rs = (b == 1600) ? 168 : 68;
    unsigned int Ra = (b == 1600) ? 192 : 96;
    unsigned int W = (b == 1600) ? 8 : 4;

    for(unsigned int Klen=16; Klen<=32; Klen++)
    for(unsigned int Nlen=0; Nlen<=200; Nlen += (Klen == 16) ? 1 : 200)
    for(unsigned int forgetFlag = 0; forgetFlag < 2; ++forgetFlag)
    for(unsigned int tagFlag = 0; tagFlag < 2; ++tagFlag)
    {
        Keyak wrap(b, nr, Pi, c, tau);
        Keyak unwrap(b, nr, Pi, c, tau);

        testKeyakStartEngine(global, wrap, unwrap, fout,
            generateSimpleRawMaterial(Klen, Klen+Nlen+0x12, 3),
            generateSimpleRawMaterial(Nlen, Klen+Nlen+0x45, 6),
            forgetFlag!=0, tagFlag!=0);
        testKeyakWrapUnwrap(global, wrap, unwrap, fout, string("ABC"), string("DEF"), false);
    }

    {
        vector<unsigned int> Alengths;
        Alengths.push_back(0);
        Alengths.push_back(1);
        Alengths.push_back(Pi*(Ra-Rs)-1);
        Alengths.push_back(Pi*(Ra-Rs));
        Alengths.push_back(Pi*(Ra-Rs)+1);
        for(unsigned int forgetFlag = 0; forgetFlag < 2; ++forgetFlag)
        for(unsigned int tagFlag = 0; tagFlag < 2; ++tagFlag)
        for(unsigned int Aleni=0; Aleni<Alengths.size(); Aleni++)
        for(unsigned int Mlen=0; Mlen<=(Rs*Pi+1); Mlen+=(Aleni==0)?1:((Pi+forgetFlag)*(W+tagFlag)+1))
        {
            unsigned int Klen = 16;
            unsigned int Nlen = (b == 1600) ? 150 : 58;
            unsigned int Alen = Alengths[Aleni];
            Keyak wrap(b, nr, Pi, c, tau);
            Keyak unwrap(b, nr, Pi, c, tau);

            testKeyakStartEngine(global, wrap, unwrap, fout,
                generateSimpleRawMaterial(Klen, 0x23+Mlen+Alen, 4),
                generateSimpleRawMaterial(Nlen, 0x56+Mlen+Alen, 7),
                forgetFlag!=0, tagFlag!=0);
            testKeyakWrapUnwrap(global, wrap, unwrap, fout,
                generateSimpleRawMaterial(Alen, 0xAB+Mlen+Alen, 3),
                generateSimpleRawMaterial(Mlen, 0xCD+Mlen+Alen, 4),
                forgetFlag!=0);
            testKeyakWrapUnwrap(global, wrap, unwrap, fout,
                generateSimpleRawMaterial(Alen, 0xCD+Mlen+Alen, 3),
                generateSimpleRawMaterial(Mlen, 0xEF+Mlen+Alen, 4),
                forgetFlag!=0);
        }
    }

    {
        vector<unsigned int> Mlengths;
        Mlengths.push_back(0);
        Mlengths.push_back(1);
        Mlengths.push_back(Pi*Rs-1);
        Mlengths.push_back(Pi*Rs);
        Mlengths.push_back(Pi*Rs+1);
        for(unsigned int forgetFlag = 0; forgetFlag < 2; ++forgetFlag)
        for(unsigned int tagFlag = 0; tagFlag < 2; ++tagFlag)
        for(unsigned int Mleni=0; Mleni<Mlengths.size(); Mleni++)
        for(unsigned int Alen=0; Alen<=(Ra*Pi+1); Alen+=(Mleni==0)?1:((Pi+forgetFlag)*(W+tagFlag)+1))
        {
            unsigned int Klen = 16;
            unsigned int Nlen = (b == 1600) ? 150 : 58;
            unsigned int Mlen = Mlengths[Mleni];
            Keyak wrap(b, nr, Pi, c, tau);
            Keyak unwrap(b, nr, Pi, c, tau);

            testKeyakStartEngine(global, wrap, unwrap, fout,
                generateSimpleRawMaterial(Klen, 0x34+Mlen+Alen, 5),
                generateSimpleRawMaterial(Nlen, 0x45+Mlen+Alen, 6),
                forgetFlag!=0, tagFlag!=0);
            testKeyakWrapUnwrap(global, wrap, unwrap, fout,
                generateSimpleRawMaterial(Alen, 0x01+Mlen+Alen, 5),
                generateSimpleRawMaterial(Mlen, 0x23+Mlen+Alen, 6),
                forgetFlag!=0);
            testKeyakWrapUnwrap(global, wrap, unwrap, fout,
                generateSimpleRawMaterial(Alen, 0x45+Mlen+Alen, 5),
                generateSimpleRawMaterial(Mlen, 0x67+Mlen+Alen, 6),
                forgetFlag!=0);
        }
    }

    {
        int forgetFlag;
        int tagFlag;
        for(forgetFlag=0; forgetFlag<2; forgetFlag++)
        for(tagFlag=0; tagFlag<2; tagFlag++) {
            unsigned int Klen = 16;
            unsigned int Nlen = (b == 1600) ? 150 : 58;
            unsigned int Alen;
            unsigned int Mlen;
            Keyak wrap(b, nr, Pi, c, tau);
            Keyak unwrap(b, nr, Pi, c, tau);

            testKeyakStartEngine(global, wrap, unwrap, fout,
                generateSimpleRawMaterial(Klen, forgetFlag*2+tagFlag, 1),
                generateSimpleRawMaterial(Nlen, forgetFlag*2+tagFlag, 2),
                forgetFlag!=0, tagFlag!=0);

            for(Alen=0; Alen<=(Ra*Pi*2); Alen+=(Alen/3+1))
            for(Mlen=0; Mlen<=(Rs*Pi*2); Mlen+=(Mlen/2+1+Alen))
            {
                testKeyakWrapUnwrap(global, wrap, unwrap, fout,
                    generateSimpleRawMaterial(Alen, 0x34+Mlen+Alen, 3),
                    generateSimpleRawMaterial(Mlen, 0x45+Mlen+Alen, 4),
                    forgetFlag!=0);
            }
        }
    }

    {
        stringstream empty, dummy, T;
        global.Wrap(empty, dummy, empty, T, false, false);
        displayByteString(fout, "+++ Global tag", T.str());
        if (T.str() != expectedGlobalTag) {
            cout << "!!! The global tag does not match." << endl;
            displayByteString(cout, "Expected", expectedGlobalTag);
            displayByteString(cout, "Actual", T.str());
        }
        assert(T.str() == expectedGlobalTag, "The global tag is incorrect.");
    }

    return (0);
}

int testAllKeyakv2Instances()
{
    int errors = 0;
    try {
        { // River Keyak
            ofstream fout("RiverKeyak.txt");
            errors += testKeyak(fout, 800, 12, 1, 256, 128, string("\x6e\xba\x81\x33\x0b\xb8\x5a\x4d\x8d\xb3\x7f\xde\x4d\x67\xcd\x0e", 16));;
        }
        { // Lake Keyak
            ofstream fout("LakeKeyak.txt");
            errors += testKeyak(fout, 1600, 12, 1, 256, 128, string("\x83\x95\xc6\x41\x22\xbb\x43\x04\x32\xd8\xb0\x29\x82\x09\xb7\x36", 16));
        }
        { // Sea Keyak
            ofstream fout("SeaKeyak.txt");
            errors += testKeyak(fout, 1600, 12, 2, 256, 128, string("\xb8\xc0\xe2\x35\x22\xcc\x1d\xe1\x4c\x22\xd0\xb8\xaf\x73\x8e\x33", 16));
        }
        { // Ocean Keyak
            ofstream fout("OceanKeyak.txt");
            errors += testKeyak(fout, 1600, 12, 4, 256, 128, string("\x70\x7c\x06\x47\xf9\xe8\x52\xb6\x00\xee\xd0\xf1\x1c\x66\xe1\x1d", 16));
        }
        { // Lunar Keyak
            ofstream fout("LunarKeyak.txt");
            errors += testKeyak(fout, 1600, 12, 8, 256, 128, string("\xb7\xec\x21\x1d\xc0\x30\xd2\x4d\x66\x70\x44\xc2\xed\x34\x52\x11", 16));
        }
    }
    catch(Exception e) {
        cout << e.reason << endl;
        errors++;
    }
    return errors;
}
