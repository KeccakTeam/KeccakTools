/*
The Keccak sponge function, designed by Guido Bertoni, Joan Daemen,
Michaël Peeters and Gilles Van Assche. For more information, feedback or
questions, please refer to our website: http://keccak.noekeon.org/

Implementation by the designers,
hereby denoted as "the implementer".

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/
*/

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include "duplex.h"
#include "Keccak.h"
#include "Keccak-f25LUT.h"
#include "Keccak-fCodeGen.h"
#include "Keccak-fDCEquations.h"
#include "Keccak-fDCLC.h"
#include "Keccak-fEquations.h"
#include "Keccak-fPropagation.h"
#include "Keccak-fTrails.h"

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

void display(const UINT8 *data, unsigned int length)
{
    for(unsigned int i=0; i<length; i++) {
        cout.width(2);
        cout.fill('0');
        cout << hex << (int)data[i] << " ";
    }
    cout << endl;
}

void testKeccakSponge(const char *message, unsigned int length)
{
    Keccak keccak;
    cout << keccak << endl;

    keccak.absorb((UINT8*)message, length);
    UINT8 squeezed[512];
    keccak.squeeze(squeezed, 4096);
    cout << "Message of length " << dec << length << endl;
    display(squeezed, 512);
}

void testKeccakSponge()
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

    testKeccakSponge(message1, message1Length);
    testKeccakSponge(message2, message2Length);
}

void testKeccakDuplex()
{
    KeccakF keccakF(1600);
    MultiRatePadding pad;
    Duplex duplex(&keccakF, &pad, 1026);
    UINT8 output[128];

    cout << duplex << endl;

    duplex.duplexing((const UINT8*)"", 0, output, 1024);
    cout << "First output: "; display(output, 128);
    duplex.duplexing((const UINT8*)"\x00", 1, output, 1024);
    cout << "Second output: "; display(output, 128);
    duplex.duplexing((const UINT8*)"\x03", 2, output, 1024);
    cout << "Third output: "; display(output, 128);
    duplex.duplexing((const UINT8*)"\x06", 3, output, 1024);
    cout << "Fourth output: "; display(output, 128);
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

void testKeccakFDCLC()
{
    KeccakFDCLC keccakFDCLC(200);
    KeccakFPropagation DC(keccakFDCLC, KeccakFPropagation::DC);
    KeccakFPropagation LC(keccakFDCLC, KeccakFPropagation::LC);
    cout << keccakFDCLC << endl;

    ofstream fout("Keccak-f-Chi-DCLC.txt");
    keccakFDCLC.displayAll(fout, &DC, &LC);
}

void displayTrails()
{
    for(unsigned int width=25; width<=200; width*=2) {
        KeccakFDCLC keccakFDCLC(width);
        {
            KeccakFPropagation DC(keccakFDCLC, KeccakFPropagation::DC);
            string fileName = DC.buildFileName("-trails");
            Trail::produceHumanReadableFile(DC, fileName);
        }
        {
            KeccakFPropagation LC(keccakFDCLC, KeccakFPropagation::LC);
            string fileName = LC.buildFileName("-trails");
            Trail::produceHumanReadableFile(LC, fileName);
        }
    }
}

void extendTrailAtTheEnd()
{
    KeccakFDCLC keccakFDCLC(200);
    KeccakFPropagation DC(keccakFDCLC, KeccakFPropagation::DC);
    cout << keccakFDCLC << endl;

    // Trail from file 'KeccakF-200-DC-trails', 4 rounds, last state removed, now 3 rounds
    istringstream sin("8 26 3 15 9 8 0 0 0 0 0 0 849108 1010842 0 1004000 0 0 0 0 0 803000 401000 0 0 0 0 0 0 80010");
    // Read the trail and display it
    Trail trail(sin);
    keccakFDCLC.checkDCTrail(trail); // optional
    trail.display(DC, cout); // for information

    // Create an output file with extensions of the trail at the end
    {
        string fileName = DC.buildFileName("-extensionAtTheEnd-trails");
        {
            ofstream fout(fileName.c_str());
            const vector<SliceValue>& lastStateOfTrail = trail.states.back();
            AffineSpaceOfStates affineSpace = DC.buildStateBase(lastStateOfTrail);
            affineSpace.display(cout); // for information
            for(SlicesAffineSpaceIterator i = affineSpace.getIterator(); !i.isEnd(); ++i) {
                Trail newTrail(trail);
                newTrail.append(*i, DC.getWeight(*i));
                newTrail.save(fout);
            }
        }
        Trail::produceHumanReadableFile(DC, fileName);
    }

    // Create an output file with extensions of the trail at the end,
    // restricting states to be in the kernel
    {
        string fileName = DC.buildFileName("-extensionAtTheEndInTheKernel-trails");
        {
            ofstream fout(fileName.c_str());
            const vector<SliceValue>& lastStateOfTrail = trail.states.back();
            AffineSpaceOfStates affineSpace = DC.buildStateBase(lastStateOfTrail);
            for(SlicesAffineSpaceIterator i = affineSpace.getIteratorWithGivenParity(0); !i.isEnd(); ++i) {
                Trail newTrail(trail);
                newTrail.append(*i, DC.getWeight(*i));
                newTrail.save(fout);
            }
        }
        Trail::produceHumanReadableFile(DC, fileName);
    }
}

void extendTrailAtTheBeginning()
{
    KeccakFDCLC keccakFDCLC(100);
    KeccakFPropagation LC(keccakFDCLC, KeccakFPropagation::LC);
    cout << keccakFDCLC << endl;

    // Trail from file 'KeccakF-100-LC-trails', 6 rounds, first state removed, now 5 rounds
    istringstream sin("4 52 5 8 4 16 14 1c 15a8000 0 0 0 0 4010 0 0 0 0 318c63 9c6318 20004 200002 5800b 80010 800800 1231802 47868 800001");
    // Read the trail and display it
    Trail trail(sin);
    keccakFDCLC.checkLCTrail(trail); // optional
    trail.display(LC, cout); // for information

    // Create an output file with extensions of the trail at the beginning,
    // restricting the first state to have weight up to 16
    {
        string fileName = LC.buildFileName("-extensionAtTheBeginning-trails");
        {
            ofstream fout(fileName.c_str());
            const vector<SliceValue>& firstStateOfTrail = trail.states.front();
            vector<SliceValue> invLambdaOfFirstStateOfTrail;
            // We have to compute λ^-1 on the state to convert it to an ouput state of χ
            LC.reverseLambda(firstStateOfTrail, invLambdaOfFirstStateOfTrail);
            for(ReverseStateIterator i = LC.getReverseStateIterator(invLambdaOfFirstStateOfTrail, 16); !i.isEnd(); ++i) {
                Trail newTrail(trail);
                newTrail.prepend(*i, LC.getWeight(*i));
                newTrail.save(fout);
            }
        }
        Trail::produceHumanReadableFile(LC, fileName);
    }
}

void generateDCTrailEquations()
{
    KeccakFDCEquations keccakFDCEq(50);
    KeccakFPropagation DC(keccakFDCEq, KeccakFPropagation::DC);
    cout << keccakFDCEq << endl;

    // Trail from file 'KeccakF-50-DC-trails', 4 rounds
    istringstream sin("2 1d 4 7 d 5 4 0 849000 84018c a0000 0 3404 4 100000");
    // Read the trail and display it
    Trail trail(sin);
    keccakFDCEq.checkDCTrail(trail); // optional
    trail.display(DC, cout); // for information

    {
        string fileName = string("DC") + keccakFDCEq.getName() + "-equations.txt";
        ofstream fout(fileName.c_str());
        keccakFDCEq.genDCEquations(fout, trail);
    }
}

int main(int argc, char *argv[])
{
    try {
        //TODO: uncomment the desired function
        //testKeccakF();
        //testKeccakSponge();
        //testKeccakDuplex();
        //genKATShortMsg_main();
        //generateEquations();
        //generateCode();
        //testKeccakF25LUT();
        //testKeccakFDCLC();
        //displayTrails();
        //extendTrailAtTheEnd();
        //extendTrailAtTheBeginning();
        //generateDCTrailEquations();
    }
    catch(SpongeException e) {
        cout << e.reason << endl;
    }
    catch(KeccakException e) {
        cout << e.reason << endl;
    }

    return EXIT_SUCCESS;
}
