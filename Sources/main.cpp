/*
Tools for the Keccak sponge function family.
Authors: Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche

This code is hereby put in the public domain. It is given as is, 
without any guarantee.

For more information, feedback or questions, please refer to our website:
http://keccak.noekeon.org/
*/

#include <cstdlib>
#include <fstream>
#include <iostream>
#include "Keccak.h"
#include "Keccak-f25LUT.h"
#include "Keccak-fCodeGen.h"
#include "Keccak-fEquations.h"

using namespace std;

void testKeccakF()
{
    KeccakF keccakF(1600);

    UINT8 state[200];
    for(unsigned int i=0; i<200; i++)
        state[i] = 0;
    keccakF(state);
    cout << "Image of the all-zero state:" << endl;
    for(unsigned int i=0; i<200; i++) {
        cout.width(2);
        cout.fill('0');
        cout << hex << (int)state[i] << " ";
    }
    cout << endl;
    cout << "Let's invert this... " << flush;
    keccakF.inverse(state);
    cout << "This should be the all-zero state again:" << endl;
    for(unsigned int i=0; i<200; i++) {
        cout.width(2);
        cout.fill('0');
        cout << hex << (int)state[i] << " ";
    }
    cout << endl;
}

void testKeccak(const char *message, unsigned int length)
{
    Keccak keccak;
    keccak.absorb((UINT8*)message, length);
    keccak.pad();
    UINT8 squeezed[512];
    keccak.squeeze(squeezed, 4096);
    cout << "Message of length " << dec << length << endl;
    for(unsigned int i=0; i<512; i++) {
        cout.width(2);
        cout.fill('0');
        cout << hex << (int)squeezed[i] << " ";
    }
    cout << endl;
}

void testKeccak()
{
    // Test messages from ShortMsgKAT.txt
    const char *message1 = "\x53\x58\x7B\x19"; // last byte aligned on LSB
    unsigned int message1Length = 29;
    const char *message2 = 
        "\x83\xAF\x34\x27\x9C\xCB\x54\x30\xFE\xBE\xC0\x7A\x81\x95\x0D\x30"
        "\xF4\xB6\x6F\x48\x48\x26\xAF\xEE\x74\x56\xF0\x07\x1A\x51\xE1\xBB"
        "\xC5\x55\x70\xB5\xCC\x7E\xC6\xF9\x30\x9C\x17\xBF\x5B\xEF\xDD\x7C"
        "\x6B\xA6\xE9\x68\xCF\x21\x8A\x2B\x34\xBD\x5C\xF9\x27\xAB\x84\x6E"
        "\x38\xA4\x0B\xBD\x81\x75\x9E\x9E\x33\x38\x10\x16\xA7\x55\xF6\x99"
        "\xDF\x35\xD6\x60\x00\x7B\x5E\xAD\xF2\x92\xFE\xEF\xB7\x35\x20\x7E"
        "\xBF\x70\xB5\xBD\x17\x83\x4F\x7B\xFA\x0E\x16\xCB\x21\x9A\xD4\xAF"
        "\x52\x4A\xB1\xEA\x37\x33\x4A\xA6\x64\x35\xE5\xD3\x97\xFC\x0A\x06"
        "\x5C\x41\x1E\xBB\xCE\x32\xC2\x40\xB9\x04\x76\xD3\x07\xCE\x80\x2E"
        "\xC8\x2C\x1C\x49\xBC\x1B\xEC\x48\xC0\x67\x5E\xC2\xA6\xC6\xF3\xED"
        "\x3E\x5B\x74\x1D\x13\x43\x70\x95\x70\x7C\x56\x5E\x10\xD8\xA2\x0B"
        "\x8C\x20\x46\x8F\xF9\x51\x4F\xCF\x31\xB4\x24\x9C\xD8\x2D\xCE\xE5"
        "\x8C\x0A\x2A\xF5\x38\xB2\x91\xA8\x7E\x33\x90\xD7\x37\x19\x1A\x07"
        "\x48\x4A\x5D\x3F\x3F\xB8\xC8\xF1\x5C\xE0\x56\xE5\xE5\xF8\xFE\xBE"
        "\x5E\x1F\xB5\x9D\x67\x40\x98\x0A\xA0\x6C\xA8\xA0\xC2\x0F\x57\x12"
        "\xB4\xCD\xE5\xD0\x32\xE9\x2A\xB8\x9F\x0A\xE1";
    unsigned int message2Length = 2008;

    testKeccak(message1, message1Length);
    testKeccak(message2, message2Length);
}

void generateEquations()
{
    for(unsigned int b=25; b<=1600; b*=2) {
        try {
            KeccakFEquations keccakF(b);
            cout << "Generating equations for " << keccakF << endl;

            string fileName = keccakF.buildFileName("Eq-", ".txt");
            ofstream fout(fileName.c_str());
            fout << "// " << keccakF << endl;
            keccakF.genComponentEquations(fout, "I", "O");
            keccakF.genRoundEquations(fout);
        }
        catch(KeccakException e) {
            cout << e.reason << endl;
        }
    }
}

void generateCode()
{
    {
        KeccakFCodeGen keccakF(1600);

        string fileName = keccakF.buildFileName("", "-64.macros");
        ofstream fout(fileName.c_str());
        keccakF.genMacroFile(fout, true);
    }

    {
        KeccakFCodeGen keccakF(1600);

        string fileName = keccakF.buildFileName("", "-simd64.macros");
        ofstream fout(fileName.c_str());
        keccakF.setOutputMacros(true);
        keccakF.genMacroFile(fout);
    }

    {
        KeccakFCodeGen keccakF(1600);

        string fileName = keccakF.buildFileName("", "-32.macros");
        ofstream fout(fileName.c_str());
        keccakF.setInterleavingFactor(2);
        keccakF.genMacroFile(fout, true);
    }
}

void testKeccakF25LUT()
{
    KeccakF25LUT keccakF25;
    cout << "Lookup table for " << keccakF25 << endl;
    for(SliceValue i=0; i<8; i++) {
        cout << "f(";
        cout.fill('0'); cout.width(7); cout << hex << i;
        cout << ") = ";
        cout.fill('0'); cout.width(7); cout << keccakF25.LUT[i];
        cout << endl;
    }
}

void genKATShortMsg_main();

int main(int argc, char *argv[])
{
    try {
        // Keccak-f[1600] with 33 rounds!
        KeccakF keccakF(1600, 33);
        cout << keccakF << endl;

        // Sponge(r=512, c=1088) using Keccak-f[1600] with 33 rounds!
        Sponge sponge(&keccakF, 1088);
        cout << sponge << endl;

        // Keccak[] = Keccak[r=1024, c=576, d=0]
        Keccak keccak;
        cout << keccak << endl;

        testKeccakF();
        testKeccak();

        //TODO: uncomment the desired function

        //genKATShortMsg_main();

        //generateEquations();

        //generateCode();

        //testKeccakF25LUT();
    }
    catch(SpongeException e) {
        cout << e.reason << endl;
    }
    catch(KeccakException e) {
        cout << e.reason << endl;
    }

    return EXIT_SUCCESS;
}
