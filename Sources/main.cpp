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
#include <string.h>
#include "duplex.h"
#include "Keccak.h"
#include "KeccakCrunchyContest.h"
#include "Keccak-f25LUT.h"
#include "Keccak-fCodeGen.h"
#include "Keccak-fDCEquations.h"
#include "Keccak-fDCLC.h"
#include "Keccak-fEquations.h"
#include "Keccak-fPropagation.h"
#include "Keccak-fTrailExtension.h"
#include "Keccak-fTrailExtensionBasedOnParity.h"
#include "Keccak-fTrails.h"
#include "Keccak-fTree.h"
#include "Keyakv2-test.h"
#include "Ketjev2-test.h"
#include "Kravatte.h"
#include "Kravatte-test.h"
#include "KravatteModes-test.h"

using namespace std;

/** Example function that uses the Keccak-f[1600] permutation and its inverse.
  */
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

/** Example function that uses the Keccak[] sponge function.
  */
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

/** Example function that uses the Keccak[] sponge function.
  */
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

/** Example function that uses the Keccak[r=1026, c=574] duplex object.
  */
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

/** Example function that generates the round equations for Keccak-f[25]
  * to Keccak-f[1600].
  */
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
            keccakF.genRoundEquations(fout, 0, keccakF.getNominalNumberOfRounds());
        }
        catch(KeccakException e) {
            cout << e.reason << endl;
        }
    }
}

/** Example function that generate C code to implement Keccak-f[1600].
  */
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

/** Example function that uses the Keccak-f[25] look-up tables.
  */
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

/** Example function that displays DC/LC propagation on rows.
  */
void testKeccakFDCLC()
{
    KeccakFDCLC keccakFDCLC(200);
    KeccakFPropagation DC(keccakFDCLC, KeccakFPropagation::DC);
    KeccakFPropagation LC(keccakFDCLC, KeccakFPropagation::LC);
    cout << keccakFDCLC << endl;

    ofstream fout("Keccak-f-Chi-DCLC.txt");
    keccakFDCLC.displayAll(fout, &DC, &LC);
}

/** Example function that produces the text files
  * in <code>Example trails</code>.
  */
void displayTrails()
{
    for(unsigned int width=25; width<=200; width*=2) {
        KeccakFDCLC keccakFDCLC(width);
        {
            KeccakFPropagation DC(keccakFDCLC, KeccakFPropagation::DC);
            string fileName = DC.buildFileName("-trailcores");
            Trail::produceHumanReadableFile(DC, fileName);
        }
        {
            KeccakFPropagation LC(keccakFDCLC, KeccakFPropagation::LC);
            string fileName = LC.buildFileName("-trailcores");
            Trail::produceHumanReadableFile(LC, fileName);
        }
    }
}

/** Example function that uses the affine base representation
  * to extend a trail forward.
  */
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

/** Example function that extends a trail backward.
  */
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

/** Example function that generates equations in GF(2) for a pair to follow
  * a given differential trail.
  */
void generateDCTrailEquations()
{
    KeccakFDCEquations keccakFDCEq(50);
    KeccakFPropagation DC(keccakFDCEq, KeccakFPropagation::DC);
    cout << keccakFDCEq << endl;

    // Trail core from file 'DCKeccakF-50-trailcores', 4 rounds
    istringstream sin("2 1d 0 c 4 7 d 5 4 3 84018c a0000 0 3404 4 100000 0");
    // Read the trail and display it
    Trail trail(sin);
    DC.specifyFirstStateArbitrarily(trail);
    DC.specifyStateAfterLastChiArbitrarily(trail);
    keccakFDCEq.checkDCTrail(trail); // optional
    trail.display(DC, cout); // for information

    {
        string fileName = string("DC") + keccakFDCEq.getName() + "-equations.txt";
        ofstream fout(fileName.c_str());
        keccakFDCEq.genDCEquations(fout, trail);
    }
}

/** Example function that generates a trail from a pair of inputs.
  * In this example, we use the messages found by I. Dinur, O. Dunkelman
  * and A. Shamir to produce a collision on Keccak[r=1088, c=512] reduced
  * to 4 rounds.
  */
void generateTrailFromDinurDunkelmanShamirCollision()
{
    const UINT8 M1[] =
        "\x32\x1c\xf3\xc4\x6d\xae\x59\x4c\xf4\xf0\x19\x5d\x4b\xe4\xc4\x25"
        "\x32\x30\x85\xd8\xf2\x12\x5e\x8d\xe2\x6e\x6e\xbb\x1e\x3b\xc3\x27"
        "\x58\x10\x09\x6c\xd5\x02\x90\xeb\x6f\xa0\xa4\x3b\xf1\xc7\x0c\x4a"
        "\x51\x5e\xb5\xcc\x83\xd9\x0d\x8d\x43\x08\x0a\x2b\xb0\xd3\x21\x9b"
        "\x75\x90\x67\x53\xd2\xde\x6d\x52\x44\x48\x29\x48\x2c\xed\xf4\x6f"
        "\x15\x2c\xce\x1a\xc7\x1d\x1c\x47\x68\x85\x09\xd4\x39\xf6\xeb\xf1"
        "\x57\xb2\xf7\xea\x87\xae\xfd\x09\xe6\x78\x88\x68\x30\xeb\x75\x48"
        "\x80\x2d\xc3\xc9\xcb\x6f\x9e\x3c\xfa\xbc\x2a\x3c\x7b\x80\xa4\xe6"
        "\xb8\x81\xb2\x2a\xb3\x32\x23";
    const unsigned int M1len = 135;
    UINT8 M2[] =
        "\xf7\x0e\xd3\xa4\x69\x8f\xbb\x80\xdf\x48\xc0\x90\xb9\x13\x72\xeb"
        "\x24\x04\x65\xa6\x3e\xf6\x65\x3a\x81\x88\x26\x8c\x1f\xb8\x51\xb6"
        "\x3c\xfa\xda\xaa\xc3\xa5\x2c\xee\xc2\xea\x78\xdb\x79\xe7\xea\xc8"
        "\x35\x9c\x2f\x44\x87\xe2\x21\x32\x5a\x7a\x01\xb3\x12\x07\x79\x90"
        "\xdc\x8b\x1c\x1b\xa8\x10\x8b\xe0\xca\x25\x9d\x9a\xac\xaa\xe7\x1b"
        "\x9c\x3e\x2f\x4e\xad\x7d\x71\x73\x5a\x01\x66\x55\xb9\xcf\x98\xa1"
        "\xc2\xa8\x1c\x5a\x8a\x34\xe3\xa0\xb1\x0b\x6c\xae\xe4\xf9\x80\x39"
        "\x91\x8b\xfa\xa4\x89\xa9\x81\x6e\xaa\xbc\xa9\x89\xf1\xf1\x2b\xe1"
        "\x95\x95\xef\x30\x45\x8b\x2e";
    const unsigned int M2len = 135;
    {
        UINT8 output[32];
        ReducedRoundKeccak keccak(1088, 512, 0, 4);
        keccak.absorb(M1, M1len*8);
        keccak.squeeze(output, 256);
        for(unsigned int i=0; i<32; i++)
            cout << hex << (int)output[i] << " ";
        cout << endl;
    }
    {
        UINT8 output[32];
        ReducedRoundKeccak keccak(1088, 512, 0, 4);
        keccak.absorb(M2, M2len*8);
        keccak.squeeze(output, 256);
        for(unsigned int i=0; i<32; i++)
            cout << hex << (int)output[i] << " ";
        cout << endl;
    }
    {
        KeccakFDCEquations keccakF(1600);
        KeccakFPropagation DC(keccakF, KeccakFPropagation::DC);
        vector<LaneValue> m1lanes, m2lanes;
        {
            UINT8 temp[200];
            memset(temp, 0, 200);
            memcpy(temp, M1, M1len);
            temp[M1len] = 0x81;
            keccakF.fromBytesToLanes(temp, m1lanes);
        }
        {
            UINT8 temp[200];
            memset(temp, 0, 200);
            memcpy(temp, M2, M2len);
            temp[M2len] = 0x81;
            keccakF.fromBytesToLanes(temp, m2lanes);
        }
        vector<SliceValue> m1slices, m2slices;
        fromLanesToSlices(m1lanes, m1slices, 64);
        fromLanesToSlices(m2lanes, m2slices, 64);
        Trail trail;
        keccakF.buildDCTrailFromPair(m1slices, m2slices, trail, 0, 4);
        {
            ofstream fout("DinurEtAl.trail");
            trail.save(fout);
        }
        Trail::produceHumanReadableFile(DC, "DinurEtAl.trail");
    }
}

/** Example function that takes trails from a file and extends them
  * forward or backward up to a given weight and given number of rounds.
  * @param  DCLC    Whether linear or differential trails are processed.
  * @param  width   The Keccak-f width.
  * @param  inFileName  The name of the file containing trails.
  * @param  nrRounds    The target number of rounds.
  * @param  maxWeight   The maximum weight of the trails to be produced.
  * @param  reverse If true, it does backward extension.
  *                 If false, it does forward extension.
  * @param  allPrefixes If true, the backward extension of trail cores
  *                     looks for all prefixes, not just trail cores.
  * @param  knownSmallWeightStateFileName   Optionally, the name of the file
  *     containing states with small weight to optimize the search. See also
  *     KeccakFTrailExtension::knownSmallWeightStates.
  * @param  maxSmallWeight  Up to which weight the small-weight state file
  *     is complete.
  */
void extendTrails(KeccakFPropagation::DCorLC DCLC, unsigned int width, const string& inFileName, unsigned int nrRounds, int maxWeight, bool reverse, bool allPrefixes=false, const string& knownSmallWeightStateFileName="", int maxSmallWeight=0)
{
    try {
        cout << "Initializing... " << flush;
        KeccakFDCLC keccakF(width);
        cout << endl;
        KeccakFTrailExtension keccakFTE(keccakF, DCLC);
        cout << keccakF << endl;

        if (knownSmallWeightStateFileName != "") {
            keccakFTE.knownSmallWeightStates = new KnownSmallWeightStates(maxSmallWeight);
            keccakFTE.knownSmallWeightStates->loadFromFile(keccakFTE, knownSmallWeightStateFileName);
            cout << "Using '" << knownSmallWeightStateFileName << "'" << endl;
        }

        try {
            TrailFileIterator trailsIn(inFileName, keccakFTE);
            cout << trailsIn << endl;
            string outFileName = inFileName + (reverse ? string("-rev") : string("-dir"));
            ofstream fout(outFileName.c_str());
            TrailSaveToFile trailsOut(fout);
            if (reverse) {
                keccakFTE.showMinimalTrails = true;
                keccakFTE.allPrefixes = allPrefixes;
                keccakFTE.backwardExtendTrails(trailsIn, trailsOut, nrRounds, maxWeight);
            }
            else {
                keccakFTE.showMinimalTrails = true;
                keccakFTE.forwardExtendTrails(trailsIn, trailsOut, nrRounds, maxWeight);
            }
            Trail::produceHumanReadableFile(keccakFTE, outFileName);
        }
        catch(TrailException e) {
            cout << e.reason << endl;
        }
    }
    catch(KeccakException e) {
        cout << e.reason << endl;
    }
}

/** Example function that uses extendTrails().
  */
void extendTrails()
{
    extendTrails(KeccakFPropagation::DC, 1600, "DCKeccakF-1600-FSE2012-3round-trailcores", 6, 75, false);
    extendTrails(KeccakFPropagation::DC, 1600, "DCKeccakF-1600-FSE2012-3round-trailcores", 6, 75, true);
    extendTrails(KeccakFPropagation::LC, 1600, "LCKeccakF-1600-trailcores", 4, 100, false);
}

// This function outputs a file with 2-round trail cores in the kernel with cost below given limit.
// An example function to use it is given below.
void traverseOrbitalTree(KeccakFPropagation::DCorLC DCLC, unsigned int width, unsigned int maxCost, unsigned int alpha, unsigned int beta)
{
    (void)DCLC;
    cout << "Initializing... " << flush;
    KeccakFDCLC keccakFDCLC(width);
    cout << endl;
    KeccakFPropagation keccakProp(keccakFDCLC, KeccakFPropagation::DC);
    cout << keccakFDCLC << endl;
    cout << "Initialized " << flush;
    cout << endl;

    // output file
    stringstream FileName;
    FileName << keccakProp.buildFileName("-TwoRoundTrailCoresInKernel-");
    FileName << "Below-";
    FileName << maxCost;
    string oFileName = FileName.str();
    ofstream fout(oFileName.c_str());
    TrailSaveToFile trailsOut(fout);

    TwoRoundTrailCoreCostFunction costF(alpha, beta);
    OrbitalsSet orbSet(width / 25);
    TwoRoundTrailCoreStack cache(keccakProp);

    OrbitalTreeIterator iterator(orbSet, cache, costF, maxCost);

    for (; !iterator.isEnd(); ++iterator) {
        TwoRoundTrailCore node = *iterator;
        node.save(fout);
    }

    Trail::produceHumanReadableFile(keccakProp, oFileName);

}

// This function outputs a file with 2-round trail cores outside the kernel with cost below given limit.
// An example function to use it is given below
void traverseRunTreeAndOrbitalTree(KeccakFPropagation::DCorLC DCLC, unsigned int width, unsigned int maxCost, unsigned int alpha, unsigned int beta)
{
    (void)DCLC;
    unsigned int laneSize = width / 25;

    cout << "Initializing... " << flush;
    KeccakFDCLC keccakFDCLC(width);
    cout << endl;
    KeccakFPropagation keccakProp(keccakFDCLC, KeccakFPropagation::DC);
    cout << keccakFDCLC << endl;
    cout << "Initialized " << flush;
    cout << endl;

    // output file
    stringstream FileName;
    FileName << keccakProp.buildFileName("-TwoRoundTrailCoresOutsideKernel-");
    FileName << "Below";
    FileName << maxCost;
    string oFileName = FileName.str();
    ofstream fout(oFileName.c_str());
    TrailSaveToFile trailsOut(fout);

    TwoRoundTrailCoreCostBoundFunction costFRun(alpha, beta);
    ColumnsSet colSet(laneSize);
    TwoRoundTrailCoreStack cacheRun(keccakProp);

    RunTreeIterator iteratorRun(colSet, cacheRun, costFRun, maxCost);

    //unsigned int counter = 0;

    for (; !iteratorRun.isEnd(); ++iteratorRun) {
        TwoRoundTrailCore nodeRun = *iteratorRun;
        unsigned int costNodeRun = alpha*nodeRun.w0 + beta*nodeRun.w1;
        bool completeNodeRun = nodeRun.complete;

        if (costNodeRun <= maxCost && completeNodeRun){
            nodeRun.save(fout);
            TwoRoundTrailCoreStack cacheOrb(keccakProp, nodeRun.stateA, nodeRun.stateB, nodeRun.w0, nodeRun.w1, completeNodeRun, nodeRun.zPeriod);
            TwoRoundTrailCoreCostFunction costFOrb(alpha, beta);

            vector<RowValue> C(nodeRun.C), D(nodeRun.D);
            vector<unsigned int> yMin(5 * laneSize, 0);

            for (unsigned int x = 0; x < 5; x++){
                for (unsigned int z = 0; z < laneSize; z++) {
                    bool odd = (getBit(C, x, z) != 0);
                    bool affected = (getBit(D, x, z) != 0);
                    if (affected) {
                        yMin[x + 5 * z] = 5; // no orbitals here
                    }
                    else{
                        if (odd){
                            for (unsigned int y = 0; y < 5; y++){
                                if (getBit(nodeRun.stateA, x, y, z) != 0){
                                    yMin[x + 5 * z] = y + 1;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            // orbital tree with parity-bare trail core at root
            OrbitalsSet orbSet(yMin, laneSize);
            OrbitalTreeIterator iteratorOrb(orbSet, cacheOrb, costFOrb, maxCost);

            for (; !iteratorOrb.isEnd(); ++iteratorOrb) {
                TwoRoundTrailCore nodeOrb = *iteratorOrb;
                nodeOrb.save(fout);
            }
        }
    }

    Trail::produceHumanReadableFile(keccakProp, oFileName);

}


// Function to perform extension in the kernel.
// Example functions to use it are given below.
void extendTrailsInTheKernel(KeccakFPropagation::DCorLC DCLC, unsigned int width, const string& inFileName, int maxWeight, unsigned int nrRounds, bool reverse, bool allPrefixes = false)
{
    try {
        cout << "Initializing... " << flush;
        KeccakFDCLC keccakF(width);
        cout << endl;
        KeccakFTrailExtensionBasedOnParity keccakFTE(keccakF, DCLC);
        cout << keccakF << endl;

        cout << "Extending... " << flush;
        TrailFileIterator trailsIn(inFileName, keccakFTE);
        cout << trailsIn << endl;

        stringstream oFileName;
        oFileName << inFileName + (reverse ? string("-revInKernel") : string("-dirInKernel"));
        oFileName << maxWeight;
        string outFileName = oFileName.str();
        ofstream fout(outFileName.c_str());
        TrailSaveToFile trailsOut(fout);

        if (reverse) {
            keccakFTE.showMinimalTrails = true;
            keccakFTE.allPrefixes = allPrefixes;
            keccakFTE.backwardExtendTrailsInTheKernel(trailsIn, trailsOut, nrRounds, maxWeight);
        }
        else {
            keccakFTE.showMinimalTrails = false;
            keccakFTE.forwardExtendTrailsInTheKernel(trailsIn, trailsOut, nrRounds, maxWeight);
        }
        Trail::produceHumanReadableFile(keccakFTE, outFileName);

    }
    catch (KeccakException e) {
        cout << e.reason << endl;
    }
}

// Function to perform extension outside the kernel.
// Example functions to use it are given below.
void extendTrailsOutsideTheKernel(KeccakFPropagation::DCorLC DCLC, unsigned int width, const string& inFileName, int maxWeight, unsigned int nrRounds, bool reverse, bool allPrefixes = false)
{
    try {
        cout << "Initializing... " << flush;
        KeccakFDCLC keccakF(width);
        cout << endl;
        KeccakFTrailExtensionBasedOnParity keccakFTE(keccakF, DCLC);
        cout << keccakF << endl;

        cout << "Extending... " << flush;
        TrailFileIterator trailsIn(inFileName, keccakFTE);
        cout << trailsIn << endl;

        stringstream oFileName;
        oFileName << inFileName + (reverse ? string("-revOutsideKernel") : string("-dirOutsideKernel"));
        oFileName << maxWeight;
        string outFileName = oFileName.str();
        ofstream fout(outFileName.c_str());
        TrailSaveToFile trailsOut(fout);

        if (reverse) {
            keccakFTE.showMinimalTrails = true;
            keccakFTE.allPrefixes = allPrefixes;
            keccakFTE.backwardExtendTrailsOutsideKernel(trailsIn, trailsOut, nrRounds, maxWeight);
        }
        else {
            keccakFTE.showMinimalTrails = false;
            keccakFTE.forwardExtendTrailsOutsideKernel(trailsIn, trailsOut, nrRounds, maxWeight);
        }
        Trail::produceHumanReadableFile(keccakFTE, outFileName);

    }
    catch (KeccakException e) {
        cout << e.reason << endl;
    }
}


// Example function to generate trail cores in the kernel
void generateTrailCoresInTheKernel()
{
    unsigned int width = 1600;
    unsigned int maxCost = 8;
    unsigned int alpha = 0; // 1 to generate states with w_0 smaller than maxCost
    unsigned int beta = 1;  // 1 to generate states with w_1 smaller than maxCost
    //unsigned int laneSize = width / 25;

    traverseOrbitalTree(KeccakFPropagation::DC, width, maxCost, alpha, beta);

}

// Example function to generate trail cores outside the kernel
void generateTrailCoresOutsideTheKernel()
{
    unsigned int width = 1600;
    unsigned int maxCost = 36;
    unsigned int alpha = 1; // 1 to generate states with w_0 smaller than maxCost
    unsigned int beta = 2;  // 1 to generate states with w_1 smaller than maxCost
    //unsigned int laneSize = width / 25;

    traverseRunTreeAndOrbitalTree(KeccakFPropagation::DC, width, maxCost, alpha, beta);

}


// Example function to perform backward extension in the kernel
void backwardExtendInKernel(){

    unsigned int width = 1600;
    unsigned int nrRounds = 3;
    unsigned int maxWeight = 36;
    string inFileName = "fileName";

    extendTrailsInTheKernel(KeccakFPropagation::DC, width, inFileName, maxWeight, nrRounds, true);

}

// Example function to perform forward extension in the kernel
void forwardExtendInKernel(){

    unsigned int width = 1600;
    unsigned int nrRounds = 3;
    unsigned int maxWeight = 36;
    string inFileName = "fileName";

    extendTrailsInTheKernel(KeccakFPropagation::DC, width, inFileName, maxWeight, nrRounds, false);

}

// Example function to perform forward extension outside the kernel
void forwardExtendOutsideKernel(){

    unsigned int width = 1600;
    unsigned int nrRounds = 3;
    unsigned int maxWeight = 36;
    string inFileName = "fileName";

    extendTrailsOutsideTheKernel(KeccakFPropagation::DC, width, inFileName, maxWeight, nrRounds, false);

}

// Example function to perform backward extension outside the kernel
void backwardExtendOutsideKernel(){

    unsigned int width = 1600;
    unsigned int nrRounds = 3;
    unsigned int maxWeight = 36;
    string inFileName = "fileName";

    extendTrailsOutsideTheKernel(KeccakFPropagation::DC, width, inFileName, maxWeight, nrRounds, true);

}

// this function computes the number of round differentials per weight for a given Keccak width
void weightDistributions(unsigned int width)
{
    unsigned int laneSize = width / 25;
    unsigned int numRows = width / 5;

    cout << "Initializing... " << flush;
    KeccakFDCLC keccakFDCLC(width);
    KeccakFPropagation keccakProp(keccakFDCLC, KeccakFPropagation::DC);
    cout << keccakFDCLC << endl;
    // output file
    stringstream outFileName;
    outFileName << keccakProp.buildFileName("-weightDistributions");
    string oFileName = outFileName.str();
    ofstream fout(oFileName.c_str());

    RowValue rowVal;
    int i;
    unsigned int numRowsPerWeight[nrRowsAndColumns]; // table where numRowsPerWeight[i] contains the number of row inputs with weight i
    double wattab0[20000], wattab1[20000]; // these tables will contain the final number of states with given weight. numbers can grow quite large, hence they are double

    // fill weight table for rows
    for(i = 0; i<nrRowsAndColumns; i++)
        numRowsPerWeight[i] = 0;
    for (rowVal = 0; rowVal<(1 << nrRowsAndColumns); rowVal++)
        numRowsPerWeight[keccakProp.getWeightRow(rowVal)]++;

    for (i = 0; i<20000; i++) wattab0[i] = wattab1[i] = 0;
    wattab1[0] = 1;
    for (unsigned int j = 0; j<numRows; j++){
        for (i = 0; i<int(j*nrRowsAndColumns) + 1; i++) wattab0[i] = wattab1[i]; // copy wattab1 into wattab0
        for (i = 0; i<int(j*nrRowsAndColumns) + 1; i++) wattab1[i] = 0; // clear wattab1
        for (i = 0; i<int(j + 1)*nrRowsAndColumns; i++){
            for (int k = 0; k<nrRowsAndColumns; k++) wattab1[i + k] += wattab0[i] * numRowsPerWeight[k]; // convolution, as many as there are "rows"
        }
    }

    for (i = 0; i<int(nrRowsAndColumns*numRows); i++)
        fout << "w: " << i << " log: " << log(wattab1[i]/laneSize) / log(2) <<  " n : " << wattab1[i]/laneSize << endl;
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
        //verifyChallenges();
        //generateTrailFromDinurDunkelmanShamirCollision();
        //extendTrails();
        //testAllKeyakv2Instances();
        //testAllKetjev2Instances();
        //backwardExtendInKernel();
        //forwardExtendInKernel();
        //backwardExtendOutsideKernel();
        //forwardExtendOutsideKernel();
        //generateTrailCoresOutsideTheKernel();
        //generateTrailCoresInTheKernel();
        //weightDistributions(200);
        //testKravatte();
        //testKravatteModes();
    }
    catch(Exception e) {
        cout << e.reason << endl;
    }

    return EXIT_SUCCESS;
}
