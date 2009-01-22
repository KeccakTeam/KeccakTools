/*
Tools for the Keccak sponge function family.
Authors: Guido Bertoni, Joan Daemen, Michaël Peeters and Gilles Van Assche

This code is hereby put in the public domain. It is given as is, 
without any guarantee.

For more information, feedback or questions, please refer to our website:
http://keccak.noekeon.org/
*/

#ifndef _KECCAKFCODEGEN_H_
#define _KECCAKFCODEGEN_H_

#include "Keccak-f.h"
#include "Keccak-fParts.h"

using namespace std;

/**
  * Class implementing code generation for the Keccak-f permutations.
  */
class KeccakFCodeGen : public KeccakF {
public:
    /**
      * The constructor. See KeccakF() for more details.
      */
    KeccakFCodeGen(unsigned int aWidth, unsigned int aNrRounds = 0);
    /**
      * Method that displays the round constants.
      */
    void displayRoundConstants();
    /**
      * Method that displays the translation offsets for ρ.
      */
    void displayRhoOffsets(bool moduloWordLength);
    /**
      * Method that displays the lane moves for π.
      */
    void displayPi();
    /**
      * Method that generates declarations for the C code produced by 
      * genCodeForRound().
      *
      * @param  fout    The output stream where the code is generated.
      * @param  interleavingFactor  The interleaving factor, see genCodeForRound().
      */
    void genDeclarations(ostream& fout, unsigned int interleavingFactor=1) const;
    /**
      * Method that generates C code to compute one round.
      * The produced code assumes that the state is stored in the variables
      * starting with letter A. It also assumes that the variables starting with
      * D are the lanes that will be XORed into the state to perform θ.
      * The generated code then XORs the D's into the A's for θ, moves
      * the lanes and rotates them into the B's to perform ρ and π.
      * The evaluation of χ is done from the B's into the variables starting
      * with E. As χ is computed, the generated code optionally also computes
      * sheet parities into 5 variables starting with C.
      * TODO: generate code for round constants
      * TODO: generate code from the C's to the D's.
      * The generated code also assumes that the lanes are complemented
      * according to patterns @a inChiMask (after the linear steps, before χ)
      * and @a outChiMask (after χ, before θ).
      *
      * @param  fout    The output stream where the code is generated.
      * @param prepareTheta A Boolean value telling whether the sheet
      *                 parities are to be computed into the C's.
      * @param  inChiMask   The lane complementing pattern at the input of χ
      *                 (or after ρ and π).
      * @param  outChiMask  The lane complementing pattern at the output of χ
      *                 (or before θ).
      * @param  interleavingFactor  The interleaving factor, i.e., the ratio
      *                 between the lane size and the target word size.
      *                 For instance, to generate 32-bit interleaved code
      *                 for Keccak-f[1600], @a interleavingFactor must be
      *                 set to 2 (=64/32). By default, the interleavingFactor
      *                 is 1, meaning no interleaving. The interleavingFactor
      *                 must divide the lane size.
      */
    void genCodeForRound(ostream& fout, bool prepareTheta, SliceValue inChiMask=0, SliceValue outChiMask=0, unsigned int interleavingFactor=1) const;
    static string buildWordName(const string& prefixSymbol, unsigned int x, unsigned int y, unsigned int z, unsigned int interleavingFactor);
    static string buildWordName(const string& prefixSymbol, unsigned int x, unsigned int z, unsigned int interleavingFactor);
protected:
    void genDeclarationsLanes(ostream& fout, const string& prefixSymbol, unsigned int interleavingFactor) const;
    void genDeclarationsSheets(ostream& fout, const string& prefixSymbol, unsigned int interleavingFactor) const;
};

#endif
