/*
KeccakTools

The Keccak sponge function, designed by Guido Bertoni, Joan Daemen,
Michaël Peeters and Gilles Van Assche. For more information, feedback or
questions, please refer to our website: http://keccak.noekeon.org/

Implementation by the designers,
hereby denoted as "the implementer".

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/
*/

#include <iostream>
#include <sstream>
#include "Keccak-fCodeGen.h"

using namespace std;

KeccakFCodeGen::KeccakFCodeGen(unsigned int aWidth)
    : KeccakF(aWidth), interleavingFactor(1),
    wordSize(laneSize), outputMacros(false), outputSubscripts(false),
    scheduleType(1)
{
}

void KeccakFCodeGen::setInterleavingFactor(unsigned int anInterleavingFactor)
{
    interleavingFactor = anInterleavingFactor;
    wordSize = laneSize/interleavingFactor;
}

void KeccakFCodeGen::setOutputMacros(bool anOutputMacros)
{
    outputMacros = anOutputMacros;
}

void KeccakFCodeGen::setOutputSubscripts(bool anOutputSubscripts)
{
    outputSubscripts = anOutputSubscripts;
}

void KeccakFCodeGen::setScheduleType(unsigned int aScheduleType)
{
    if ((aScheduleType >= 1) && (aScheduleType <= 2))
        scheduleType = aScheduleType;
}

void KeccakFCodeGen::displayRoundConstants()
{
    for(unsigned int i=0; i<roundConstants.size(); ++i) {
        cout << "KeccakF" << dec << width << "RoundConstants[";
        cout.width(2); cout.fill(' ');
        cout << dec << i << "] = ";
        cout.width(16); cout.fill('0');
        cout << hex << roundConstants[i] << endl;
    }
}

void KeccakFCodeGen::displayRhoOffsets(bool moduloWordLength)
{
    int offset = 2;

    cout << "ρ:" << endl;
    cout << "col  |";
    for(unsigned int sx=0; sx<5; sx++) {
        unsigned int x = index(sx-offset);
        cout.width(4); cout.fill(' ');
        cout << dec << x;
    }
    cout << endl;
    cout << "-----+--------------------" << endl;
    for(unsigned int sy=0; sy<5; sy++) {
        unsigned int y = index(5-1-sy-offset);
        cout << dec << "row " << dec << y << "|";
        for(unsigned int sx=0; sx<5; sx++) {
            unsigned int x = index(sx-offset);
            cout.width(4); cout.fill(' ');
            if (moduloWordLength)
                cout << dec << (rhoOffsets[index(x, y)] % laneSize);
            else
                cout << dec << (rhoOffsets[index(x, y)]);
        }
        cout << endl;
    }
    cout << endl;
}

void KeccakFCodeGen::displayPi()
{
    for(unsigned int x=0; x<5; x++)
    for(unsigned int y=0; y<5; y++) {
        cout << dec << "(" << x << "," << y << ") goes to ";
        unsigned int X, Y;
        pi(x, y, X, Y);
        cout << dec << "(" << X << "," << Y << ").";
        cout << endl;
    }
}

string KeccakFCodeGen::buildWordName(const string& prefixSymbol, unsigned int x, unsigned int y, unsigned int z) const
{
    if (outputSubscripts) {
        stringstream temp;
        temp << prefixSymbol << "[" << dec << index(x, y)*interleavingFactor+z << "]";
        return temp.str();
    }
    else
        return buildBitName(laneName(prefixSymbol, x, y), interleavingFactor, z);
}

string KeccakFCodeGen::buildWordName(const string& prefixSymbol, unsigned int x, unsigned int z) const
{
    if (outputSubscripts) {
        stringstream temp;
        temp << prefixSymbol << "[" << dec << x*interleavingFactor+z << "]";
        return temp.str();
    }
    else
        return buildBitName(sheetName(prefixSymbol, x), interleavingFactor, z);
}

string KeccakFCodeGen::buildWordName(const string& prefixSymbol, unsigned int x) const
{
    if (outputSubscripts) {
        stringstream temp;
        temp << prefixSymbol << "[" << dec << x << "]";
        return temp.str();
    }
    else
        return sheetName(prefixSymbol, x);
}

void KeccakFCodeGen::genDeclarations(ostream& fout) const
{
    genDeclarationsLanes(fout, "A");
    genDeclarationsLanes(fout, "B");
    genDeclarationsSheets(fout, "C");
    genDeclarationsSheets(fout, "D");
    genDeclarationsLanes(fout, "E");
    fout << endl;
}

void KeccakFCodeGen::genDeclarationsLanes(ostream& fout, const string& prefixSymbol) const
{
    for(unsigned int y=0; y<5; y++)
    for(unsigned int z=0; z<interleavingFactor; z++) {
        fout << "    ";
        fout << (outputMacros ? "V" : "UINT") << dec << wordSize << " ";
        for(unsigned int x=0; x<5; x++) {
            fout << buildWordName(prefixSymbol, x, y, z);
            if (x < 4)
                fout << ", ";
        }
        fout << "; \\" << endl;
    }
}

void KeccakFCodeGen::genDeclarationsSheets(ostream& fout, const string& prefixSymbol) const
{
    for(unsigned int z=0; z<interleavingFactor; z++) {
        fout << "    ";
        fout << (outputMacros ? "V" : "UINT") << dec << wordSize << " ";
        for(unsigned int x=0; x<5; x++) {
            fout << buildWordName(prefixSymbol, x, z);
            if (x < 4)
                fout << ", ";
        }
        fout << "; \\" << endl;
    }
}

unsigned int KeccakFCodeGen::schedule(unsigned int i) const
{
    if (scheduleType == 1)
        return i;
    else if (scheduleType == 2) {
        const unsigned int theSchedule[10] = { 0, 1, 2, 5, 3, 6, 4, 7, 8, 9 };
        return theSchedule[i];
    }
    else
        return i;
}

void KeccakFCodeGen::genCodePlanePerPlane(ostream& fout, bool prepareTheta,
                                     SliceValue inChiMask, SliceValue outChiMask,
                                     string A, string B, string C,
                                     string D, string E, string header) const
{
    fout << "// --- Code for round";
    if (prepareTheta) fout << ", with prepare-theta";
    if (outChiMask != 0) {
        fout << " (lane complementing pattern '";
        for(unsigned int y=0; y<5; y++)
        for(unsigned int x=0; x<5; x++) {
            if (getRowFromSlice(outChiMask, y) & (1<<x))
                fout << laneName("", x, y);
        }
        fout << "')";
    }
    fout << endl;
    fout << "// --- ";
    if (interleavingFactor > 1)
        fout << "using factor " << interleavingFactor << " interleaving, ";
    fout << dec << laneSize << "-bit lanes mapped to " << wordSize << "-bit words";
    fout << endl;

    if (header.size() > 0)
        fout << header << endl;

    // θ: from C to D
    for(unsigned int x=0; x<5; x++)  {
        fout << "    " << buildWordName(D, x, 0) << " = ";
        fout << strXOR(buildWordName(C, (x+4)%5, 0), strROL(buildWordName(C, (x+1)%5, interleavingFactor-1), 1));
        fout << "; \\" << endl;
        for(unsigned int zeta=1; zeta<interleavingFactor; zeta++) {
            fout << "    " << buildWordName(D, x, zeta) << " = ";
            fout << strXOR(buildWordName(C, (x+4)%5, zeta), buildWordName(C, (x+1)%5, zeta-1));
            fout << "; \\" << endl;
        }
    }

    fout << "\\" << endl;

    for(unsigned int y=0; y<5; y++)
    for(unsigned int zeta=0; zeta<interleavingFactor; zeta++) {
        for(unsigned int i=0; i<10; i++)
        for(unsigned int x=0; x<5; x++) {
            unsigned int xprime, yprime;
            inversePi(x, y, xprime, yprime);
            unsigned int rModS = rhoOffsets[index(xprime, yprime)] % interleavingFactor;
            unsigned int zetaprime = (interleavingFactor + zeta - rModS) % interleavingFactor;
            unsigned int r = (rhoOffsets[index(xprime, yprime)] / interleavingFactor) % wordSize
                 + ((zeta < rModS) ? 1 : 0);
            unsigned j = schedule(i);

            if (j == x) {
                // θ
                fout << "    ";
                fout << strXOReq(buildWordName(A, xprime, yprime, zetaprime), buildWordName(D, xprime, zetaprime));
                fout << "; \\" << endl;

                // ρ then π
                fout << "    " << buildWordName(B, x, y, zeta) << " = ";
                fout << strROL(buildWordName(A, xprime, yprime, zetaprime), r);
                fout << "; \\" << endl;
            }

            if (j == (x+5)) {
                bool M0 = ((((outChiMask ^ inChiMask)>>(x+5*y))&1) == 1);
                bool M1 = (((inChiMask>>((x+1)%5 + 5*y))&1) == 1);
                bool M2 = (((inChiMask>>((x+2)%5 + 5*y))&1) == 1);
                bool LC1 = (M1==M2) && (M0 == M1);
                bool LC2 = (M1==M2) && (M0 != M1);
                bool LOR = ((!M1) && M2) || (M0 && (M1==M2));
                bool LC0 = !LOR==M0;
                // χ
                fout << "    " << buildWordName(E, x, y, zeta) << " = ";
                fout << strXOR(
                    strNOT(buildWordName(B, x, y, zeta), LC0),
                    strANDORnot(buildWordName(B, index(x+1), y, zeta), buildWordName(B, index(x+2), y, zeta), LC1, LC2, LOR));
                fout << "; \\" << endl;

                if ((x == 0) && (y == 0)) {
                    // ι
                    stringstream str;
                    str << "KeccakF" << width << "RoundConstants";
                    if (interleavingFactor > 1)
                        str << "_int" << interleavingFactor << "_" << zeta;
                    str << "[i]";
                    fout << "    ";
                    fout << strXOReq(buildWordName(E, x, y, zeta), strConst(str.str()));
                    fout << "; \\" << endl;
                }

                if (prepareTheta) {
                    // Prepare θ
                    if (y == 0) {
                        fout << "    " << buildWordName(C, x, zeta);
                        fout << " = ";
                        fout << buildWordName(E, x, y, zeta) << "; \\" << endl;
                    }
                    else {
                        fout << "    ";
                        fout << strXOReq(buildWordName(C, x, zeta), buildWordName(E, x, y, zeta));
                        fout << "; \\" << endl;
                    }
                }
            }
        }
        fout << "\\" << endl;
    }
    fout << endl;
}

unsigned int getInPlaceY(unsigned int i, unsigned int x, unsigned int y)
{
    switch(i%4) {
        case 0:
            return y;
            break;
        case 1:
            return (x + 2*y)%5;
            break;
        case 2:
            return (3*x + 4*y)%5;
            break;
        case 3:
            return (2*x + 3*y)%5;
            break;
        default:
            return 0;
    }
}

unsigned int getXplus2Y(unsigned int x, unsigned int y)
{
    return (x + 2*y)%5;
}

void KeccakFCodeGen::genCodeInPlace(ostream& fout,
    bool earlyParity,
    SliceValue inChiMask, SliceValue outChiMask,
    string A, string B, string C, string D, string header) const
{
    if (interleavingFactor > 2)
        throw KeccakException("This routine is only available for interleaving factors 1 and 2.");
    fout << "// --- Code for 4 rounds";
    if (outChiMask != 0) {
        fout << " (lane complementing pattern '";
        for(unsigned int y=0; y<5; y++)
        for(unsigned int x=0; x<5; x++) {
            if (getRowFromSlice(outChiMask, y) & (1<<x))
                fout << laneName("", x, y);
        }
        fout << "')";
    }
    fout << endl;
    fout << "// --- ";
    if (interleavingFactor > 1)
        fout << "using factor " << interleavingFactor << " interleaving, ";
    fout << dec << laneSize << "-bit lanes mapped to " << wordSize << "-bit words";
    fout << endl;

    if (header.size() > 0)
        fout << header << endl;

    vector<unsigned int> O(25, 0);

    for(unsigned int i=0; i<4; i++) {
        if (!earlyParity) {
            for(unsigned int x=0; x<5; x++)
            for(unsigned int zeta=0; zeta<interleavingFactor; zeta++) {
                fout << "    " << buildWordName(C, x, zeta) << " = ";
                fout <<
                    strXOR(buildWordName(A, x, getInPlaceY(i, x, 0), (zeta+O[index(x, getInPlaceY(i, x, 0))])%interleavingFactor),
                    strXOR(buildWordName(A, x, getInPlaceY(i, x, 1), (zeta+O[index(x, getInPlaceY(i, x, 1))])%interleavingFactor),
                    strXOR(buildWordName(A, x, getInPlaceY(i, x, 2), (zeta+O[index(x, getInPlaceY(i, x, 2))])%interleavingFactor),
                    strXOR(buildWordName(A, x, getInPlaceY(i, x, 3), (zeta+O[index(x, getInPlaceY(i, x, 3))])%interleavingFactor),
                           buildWordName(A, x, getInPlaceY(i, x, 4), (zeta+O[index(x, getInPlaceY(i, x, 4))])%interleavingFactor)))));
                fout << "; \\" << endl;
            }
        }

        // θ: from C to D
        for(unsigned int x=0; x<5; x++)  {
            fout << "    " << buildWordName(D, x, 0) << " = ";
            fout << strXOR(buildWordName(C, (x+4)%5, 0), strROL(buildWordName(C, (x+1)%5, interleavingFactor-1), 1));
            fout << "; \\" << endl;
            for(unsigned int zeta=1; zeta<interleavingFactor; zeta++) {
                fout << "    " << buildWordName(D, x, zeta) << " = ";
                fout << strXOR(buildWordName(C, (x+4)%5, zeta), buildWordName(C, (x+1)%5, zeta-1));
                fout << "; \\" << endl;
            }
        }

        fout << "\\" << endl;

        vector<unsigned int> nextO(25, 0);
        for(unsigned int y=0; y<5; y++)
        for(unsigned int zeta=0; zeta<interleavingFactor; zeta++) {
            vector<unsigned int> zetapp(5, 0);
            for(unsigned int x=0; x<5; x++) {
                unsigned int xpp = x;
                unsigned int ypp = getInPlaceY(i+1, x, y);
                unsigned int rModS = rhoOffsets[index(x, getXplus2Y(x, y))] % interleavingFactor;
                unsigned int zetaprime = (interleavingFactor + zeta - rModS) % interleavingFactor;
                unsigned int r = (rhoOffsets[index(x, getXplus2Y(x, y))] / interleavingFactor) % wordSize
                    + ((zeta < rModS) ? 1 : 0);
                zetapp[xpp] = (zetaprime + O[index(xpp, ypp)])%interleavingFactor;

                // θ then ρ then π
                fout << "    " << buildWordName(B, getXplus2Y(x, y)) << " = ";
                fout << strROL((string("(")+strXOR(buildWordName(A, xpp, ypp, zetapp[xpp]),
                    buildWordName(D, x, zetaprime))+")"), r);
                fout << "; \\" << endl;
            }
            for(unsigned int x=0; x<5; x++) {
                unsigned int xpp = x;
                unsigned int ypp = getInPlaceY(i+1, x, y);
                bool M0 = ((((outChiMask ^ inChiMask)>>(x+5*y))&1) == 1);
                bool M1 = (((inChiMask>>((x+1)%5 + 5*y))&1) == 1);
                bool M2 = (((inChiMask>>((x+2)%5 + 5*y))&1) == 1);
                bool LC1 = (M1==M2) && (M0 == M1);
                bool LC2 = (M1==M2) && (M0 != M1);
                bool LOR = ((!M1) && M2) || (M0 && (M1==M2));
                bool LC0 = !LOR==M0;
                // χ
                fout << "    " << buildWordName(A, xpp, ypp, zetapp[xpp]) << " = ";
                fout << strXOR(
                    strNOT(buildWordName(B, x), LC0),
                    strANDORnot(buildWordName(B, index(x+1)), buildWordName(B, index(x+2)), LC1, LC2, LOR));
                fout << "; \\" << endl;

                if ((x == 0) && (y == 0)) {
                    // ι
                    stringstream str;
                    str << "KeccakF" << width << "RoundConstants";
                    if (interleavingFactor > 1)
                        str << "_int" << interleavingFactor << "_" << zeta;
                    str << "[i+" << dec << i << "]";
                    fout << "    ";
                    fout << strXOReq(buildWordName(A, xpp, ypp, zetapp[xpp]), strConst(str.str()));
                    fout << "; \\" << endl;
                }
                if (earlyParity) {
                    // Prepare θ
                    if (y == 0) {
                        fout << "    " << buildWordName(C, x, zeta);
                        fout << " = ";
                        fout << buildWordName(A, xpp, ypp, zetapp[xpp]) << "; \\" << endl;
                    }
                    else {
                        fout << "    ";
                        fout << strXOReq(buildWordName(C, x, zeta), buildWordName(A, xpp, ypp, zetapp[xpp]));
                        fout << "; \\" << endl;
                    }
                }
            }
            fout << "\\" << endl;
            if (zeta == 0) {
                for(unsigned int x=0; x<5; x++)
                    nextO[index(x, getInPlaceY(i+1, x, y))] = zetapp[x];
            }
        }
        O = nextO;
    }
    fout << endl;
}

void KeccakFCodeGen::genCodeForPrepareTheta(ostream& fout, string A, string C) const
{
    for(unsigned int x=0; x<5; x++)
    for(unsigned int z=0; z<interleavingFactor; z++) {
        fout << "    " << buildWordName(C, x, z) << " = ";
        fout <<
            strXOR(buildWordName(A, x, 0, z),
            strXOR(buildWordName(A, x, 1, z),
            strXOR(buildWordName(A, x, 2, z),
            strXOR(buildWordName(A, x, 3, z), buildWordName(A, x, 4, z)))));
        fout << "; \\" << endl;
    }
    fout << endl;
}

void KeccakFCodeGen::genRoundConstants(ostream& fout) const
{
    vector<vector<LaneValue> > interleavedRC;

    for(vector<LaneValue>::const_iterator i=roundConstants.begin(); i!=roundConstants.end(); ++i) {
        vector<LaneValue> iRC(interleavingFactor, 0);
        for(unsigned int z=0; z<laneSize; z++)
            if (((*i) & ((LaneValue)1 << z)) != 0)
                iRC[z%interleavingFactor] |= ((LaneValue)1 << (z/interleavingFactor));
        interleavedRC.push_back(iRC);
    }

    for(unsigned int z=0; z<interleavingFactor; z++) {
        fout << "const UINT" << dec << wordSize << " ";
        fout << "KeccakF" << width << "RoundConstants";
        if (interleavingFactor > 1)
            fout << "_int" << interleavingFactor << "_" << z;
        fout << "[" << dec << interleavedRC.size() << "] = {" << endl;
        for(unsigned int i=0; i<interleavedRC.size(); i++) {
            if (i > 0)
                fout << "," << endl;
            fout << "    0x";
            fout.fill('0'); fout.width((wordSize+3)/4);
            fout << hex << interleavedRC[i][z];
            if (wordSize == 64)
                fout << "ULL";
            else if (wordSize == 32)
                fout << "UL";
        }
        fout << " };" << endl;
        fout << endl;
    }
}

void KeccakFCodeGen::genCopyFromStateAndXor(ostream& fout, unsigned int bitsToXor, string A, string state, string input) const
{
    unsigned int i=0;
    for(unsigned int y=0; y<5; y++)
    for(unsigned int x=0; x<5; x++)
    for(unsigned int z=0; z<interleavingFactor; z++) {
        if (outputMacros) {
            fout << "    " << buildWordName(A, x, y, z);
            fout << " = ";
            if (i*wordSize < bitsToXor)
                fout << "XOR" << dec << wordSize << "(";
            fout << "LOAD" << dec << wordSize << "(" << state << "[";
            fout.width(2); fout.fill(' '); fout << dec << i << "]";
            fout << ")";
            if (i*wordSize < bitsToXor) {
                fout << ", ";
                fout << "LOAD" << dec << wordSize << "(" << input << "[";
                fout.width(2); fout.fill(' '); fout << dec << i << "]";
                fout << ")";
                fout << ")";
            }
        }
        else {
            fout << "    " << buildWordName(A, x, y, z);
            fout << " = " << state << "[";
            fout.width(2); fout.fill(' '); fout << dec << i << "]";
            if (i*wordSize < bitsToXor) {
                fout << "^" << input << "[";
                fout.width(2); fout.fill(' '); fout << dec << i << "]";
            }
        }
        fout << "; \\" << endl;
        i++;
    }
    fout << endl;
}

void KeccakFCodeGen::genCopyToState(ostream& fout, string A, string state, string input) const
{
    (void)input;
    unsigned int i=0;
    for(unsigned int y=0; y<5; y++)
    for(unsigned int x=0; x<5; x++)
    for(unsigned int z=0; z<interleavingFactor; z++) {
        if (outputMacros) {
            fout << "    STORE" << dec << wordSize << "(" << state << "[";
            fout.width(2); fout.fill(' '); fout << dec << i << "], ";
            fout << buildWordName(A, x, y, z);
            fout << ")";
        }
        else {
            fout << "    ";
            fout << state << "[";
            fout.width(2); fout.fill(' '); fout << dec << i << "] = ";
            fout << buildWordName(A, x, y, z);
        }
        fout << "; \\" << endl;
        i++;
    }
    fout << endl;
}

void KeccakFCodeGen::genCopyStateVariables(ostream& fout, string X, string Y) const
{
    for(unsigned int y=0; y<5; y++)
    for(unsigned int x=0; x<5; x++)
    for(unsigned int z=0; z<interleavingFactor; z++) {
        fout << "    ";
        fout << buildWordName(X, x, y, z);
        fout << " = ";
        fout << buildWordName(Y, x, y, z);
        fout << "; \\" << endl;
    }
    fout << endl;
}

void KeccakFCodeGen::genMacroFile(ostream& fout, bool laneComplementing) const
{
    fout << "/*" << endl;
    fout << "Code automatically generated by KeccakTools!" << endl;
    fout << endl;
    fout << "The Keccak sponge function, designed by Guido Bertoni, Joan Daemen," << endl;
    fout << "Micha\xC3\xABl Peeters and Gilles Van Assche. For more information, feedback or" << endl;
    fout << "questions, please refer to our website: http://keccak.noekeon.org/" << endl;
    fout << endl;
    fout << "Implementation by the designers," << endl;
    fout << "hereby denoted as \"the implementer\"." << endl;
    fout << endl;
    fout << "To the extent possible under law, the implementer has waived all copyright" << endl;
    fout << "and related or neighboring rights to the source code in this file." << endl;
    fout << "http://creativecommons.org/publicdomain/zero/1.0/" << endl;
    fout << "*/" << endl;
    fout << endl;
    fout << "#define declareABCDE \\" << endl;
    genDeclarations(fout);
    fout << "#define prepareTheta \\" << endl;
    genCodeForPrepareTheta(fout);
    if (laneComplementing) {
        const SliceValue inChiMask = 0x9d14ad;
        const SliceValue outChiMask = 0x121106; // see "Keccak implementation overview", Section "The lane complementing transform"
        fout << "#ifdef UseBebigokimisa" << endl;
        genCodePlanePerPlane(fout, true, inChiMask, outChiMask,
            "A##", "B", "C", "D", "E##",
            "#define thetaRhoPiChiIotaPrepareTheta(i, A, E) \\");
        genCodePlanePerPlane(fout, false, inChiMask, outChiMask,
            "A##", "B", "C", "D", "E##",
            "#define thetaRhoPiChiIota(i, A, E) \\");
        fout << "#else // UseBebigokimisa" << endl;
    }
    genCodePlanePerPlane(fout, true, 0, 0,
        "A##", "B", "C", "D", "E##",
        "#define thetaRhoPiChiIotaPrepareTheta(i, A, E) \\");
    genCodePlanePerPlane(fout, false, 0, 0,
        "A##", "B", "C", "D", "E##",
        "#define thetaRhoPiChiIota(i, A, E) \\");
    if (laneComplementing)
        fout << "#endif // UseBebigokimisa" << endl << endl;
    genRoundConstants(fout);
    if (width == 1600) {
        fout << "#define copyFromStateAndXor576bits(X, state, input) \\" << endl;
        genCopyFromStateAndXor(fout, 576);
        fout << "#define copyFromStateAndXor832bits(X, state, input) \\" << endl;
        genCopyFromStateAndXor(fout, 832);
        fout << "#define copyFromStateAndXor1024bits(X, state, input) \\" << endl;
        genCopyFromStateAndXor(fout, 1024);
        fout << "#define copyFromStateAndXor1088bits(X, state, input) \\" << endl;
        genCopyFromStateAndXor(fout, 1088);
        fout << "#define copyFromStateAndXor1152bits(X, state, input) \\" << endl;
        genCopyFromStateAndXor(fout, 1152);
        fout << "#define copyFromStateAndXor1344bits(X, state, input) \\" << endl;
        genCopyFromStateAndXor(fout, 1344);
    }
    fout << "#define copyFromState(X, state) \\" << endl;
    genCopyFromStateAndXor(fout, 0);
    fout << "#define copyToState(state, X) \\" << endl;
    genCopyToState(fout);
    fout << "#define copyStateVariables(X, Y) \\" << endl;
    genCopyStateVariables(fout);
}

string KeccakFCodeGen::strROL(const string& symbol, unsigned int amount) const
{
    stringstream str;

    if (amount > 0)
        str << "ROL" << dec << wordSize << "(";
    str << symbol;
    if (amount > 0) {
        str << ", " << dec;
        str << amount << ")";
    }
    return str.str();
}

string KeccakFCodeGen::strXOR(const string& A, const string& B) const
{
    stringstream str;

    if (outputMacros) {
        str << "XOR" << dec << wordSize << "(";
        str << A << ", " << B << ")";
    }
    else
        str << A << "^" << B;
    return str.str();
}

string KeccakFCodeGen::strXOReq(const string& A, const string& B) const
{
    stringstream str;

    if (outputMacros) {
        str << "XOReq" << dec << wordSize << "(";
        str << A << ", " << B << ")";
    }
    else
        str << A << " ^= " << B;
    return str.str();
}

string KeccakFCodeGen::strANDORnot(const string& A, const string& B, bool LC1, bool LC2, bool LOR) const
{
    stringstream str;

    if (outputMacros) {
        str << (LOR ? "OR" : "AND") << (LC1 ? "n" : "u") << (LC2 ? "n" : "u")
            << dec << wordSize << "(";
        str << A << ", " << B;
        str << ")";
    }
    else {
        str << "(";
        str << (LC1 ? "(~" : "  ") << A << (LC1 ? ")" : " ");
        str << (LOR ? "|" : "&");
        str << (LC2 ? "(~" : "  ") << B << (LC2 ? ")" : " ");
        str << ")";
    }
    return str.str();
}

string KeccakFCodeGen::strNOT(const string& A, bool complement) const
{
    stringstream str;

    if (outputMacros) {
        if (complement)
            str << "NOT" << dec << wordSize << "(";
        str << A;
        if (complement)
            str << ")";
    }
    else {
        if (complement)
            str << "(~";
        else
            str << "  ";
        str << A;
        if (complement)
            str << ")";
        else
            str << " ";
    }
    return str.str();
}

string KeccakFCodeGen::strConst(const string& A) const
{
    if (outputMacros) {
        stringstream str;

        str << "CONST" << dec << wordSize << "(";
        str << A;
        str << ")";
        return str.str();
    }
    else
        return A;
}

string KeccakFCodeGen::getName() const
{
    stringstream a;
    a << "KeccakF-" << dec << width;
    return a.str();
}
