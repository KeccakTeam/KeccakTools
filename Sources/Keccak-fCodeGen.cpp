/*
Tools for the Keccak sponge function family.
Authors: Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche

This code is hereby put in the public domain. It is given as is, 
without any guarantee.

For more information, feedback or questions, please refer to our website:
http://keccak.noekeon.org/
*/

#include <iostream>
#include "Keccak-fCodeGen.h"

using namespace std;

KeccakFCodeGen::KeccakFCodeGen(unsigned int aWidth, unsigned int aNrRounds)
    : KeccakF(aWidth, aNrRounds)
{
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

    cout << "Rho:" << endl;
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

string KeccakFCodeGen::buildWordName(const string& prefixSymbol, unsigned int x, unsigned int y, unsigned int z, unsigned int interleavingFactor)
{
    return buildBitName(laneName(prefixSymbol, x, y), interleavingFactor, z);
}

string KeccakFCodeGen::buildWordName(const string& prefixSymbol, unsigned int x, unsigned int z, unsigned int interleavingFactor)
{
    return buildBitName(sheetName(prefixSymbol, x), interleavingFactor, z);
}

void KeccakFCodeGen::genDeclarations(ostream& fout, unsigned int interleavingFactor) const
{
    fout << "// --- Declarations" << endl;
    genDeclarationsLanes(fout, "A", interleavingFactor);
    genDeclarationsLanes(fout, "B", interleavingFactor);
    genDeclarationsSheets(fout, "C", interleavingFactor);
    genDeclarationsSheets(fout, "D", interleavingFactor);
    genDeclarationsLanes(fout, "E", interleavingFactor);
}

void KeccakFCodeGen::genDeclarationsLanes(ostream& fout, const string& prefixSymbol, unsigned int interleavingFactor) const
{
    for(unsigned int y=0; y<5; y++)
    for(unsigned int z=0; z<interleavingFactor; z++) {
        fout << "    UINT" << dec << (laneSize/interleavingFactor) << " ";
        for(unsigned int x=0; x<5; x++) {
            fout << buildWordName(prefixSymbol, x, y, z, interleavingFactor);
            if (x < 4)
                fout << ", ";
        }
        fout << "; \\" << endl;
    }
}

void KeccakFCodeGen::genDeclarationsSheets(ostream& fout, const string& prefixSymbol, unsigned int interleavingFactor) const
{
    for(unsigned int z=0; z<interleavingFactor; z++) {
        fout << "    UINT" << dec << (laneSize/interleavingFactor) << " ";
        for(unsigned int x=0; x<5; x++) {
            fout << buildWordName(prefixSymbol, x, z, interleavingFactor);
            if (x < 4)
                fout << ", ";
        }
        fout << "; \\" << endl;
    }
}

unsigned int schedule(unsigned int i)
{
    switch(i) {
    case 3:
        return 5;
        break;
    case 4:
        return  3;
        break;
    case 5:
        return  6;
        break;
    case 6:
        return  4;
        break;
    default:
        return i;
    }
}

void KeccakFCodeGen::genCodeForRound(ostream& fout, bool prepareTheta, 
                                     SliceValue inChiMask, SliceValue outChiMask, 
                                     unsigned int interleavingFactor) const
{
    fout << "// --- Theta Rho Pi Chi";
    if (prepareTheta) fout << " Iota Prepare-theta";
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
    fout << dec << laneSize << "-bit lanes mapped to " << laneSize/interleavingFactor << "-bit words";
    fout << endl;
    for(unsigned int Y=0; Y<5; Y++)
    for(unsigned int Z=0; Z<interleavingFactor; Z++) {
        for(unsigned int i=0; i<10; i++)
        for(unsigned int X=0; X<5; X++) {
            unsigned int x, y;
            inversePi(X, Y, x, y);
            unsigned int r = rhoOffsets[index(x, y)] % laneSize;
            unsigned int rMod = r%interleavingFactor;
            unsigned int rDiv = r/interleavingFactor;
            unsigned int z = (interleavingFactor + Z - rMod)%interleavingFactor;
            unsigned int rInt = rDiv + (z + rMod)/interleavingFactor;
            unsigned j = schedule(i);
            const unsigned int iF = interleavingFactor;

            if (j == X) {
                // Theta
                fout << "    " << buildWordName("A", x, y, z, iF) << " ^= ";
                fout << buildWordName("D", x, z, iF) << "; \\" << endl;

                // Rho Pi
                fout << "    " << buildWordName("B", X, Y, Z, iF) << " = ";
                if (rInt > 0)
                    fout << "ROL" << dec << laneSize/interleavingFactor << "(";
                fout << buildWordName("A", x, y, z, iF);
                if (rInt > 0) {
                    fout << ", " << dec;
                    fout.width(2); fout.fill(' ');
                    fout << rInt << ")";
                }
                fout << "; \\" << endl;
            }

            if (j == (X+5)) {
                bool M0 = ((((outChiMask ^ inChiMask)>>(X+5*Y))&1) == 1);
                bool M1 = (((inChiMask>>((X+1)%5 + 5*Y))&1) == 1);
                bool M2 = (((inChiMask>>((X+2)%5 + 5*Y))&1) == 1);
                bool LC1 = (M1==M2) && (M0 == M1);
                bool LC2 = (M1==M2) && (M0 != M1);
                bool LOR = ((!M1) && M2) || (M0 && (M1==M2));
                bool LC0 = !LOR==M0; 
                // Chi
                fout << "    " << buildWordName("E", X, Y, Z, iF) << " = ";
                fout << (LC0 ? "~" : " ") << buildWordName("B", X, Y, Z, iF);
                fout << " ^ (";
                fout << (LC1 ? "(~" : "  ")
                    << buildWordName("B", index(X+1), Y, Z, iF)
                    << (LC1 ? ")" : " ");
                fout << (LOR ? "|" : "&");
                fout << (LC2 ? "(~" : "  ")
                    << buildWordName("B", index(X+2), Y, Z, iF)
                    << (LC2 ? ")" : " ");
                fout << "); \\" << endl;

                if ((X == 0) && (Y == 0)) {
                    // Iota
                    fout << "    " << buildWordName("E", X, Y, Z, iF) << " ^= ";
                    fout << "KeccakF" << width << "RoundConstants";
                    if (interleavingFactor > 1)
                        fout << "_int" << interleavingFactor << "_" << Z;
                    fout << "[i]; \\" << endl;
                }

                if (prepareTheta) {
                    // Prepare Theta
                    fout << "    " << buildWordName("C", X, Z, iF);
                    if (Y == 0)
                        fout << " = ";
                    else
                        fout << " ^= ";
                    fout << buildWordName("E", X, Y, Z, iF) << "; \\" << endl;
                }
            }
        }
        fout << "\\" << endl;
    }
}
