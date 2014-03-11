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

#ifndef _KECCAKFEQUATIONS_H_
#define _KECCAKFEQUATIONS_H_

#include <iostream>
#include <vector>
#include "Keccak-f.h"

using namespace std;

class SymbolicLane;

/**
  * Class implementing equation generation for the Keccak-<i>f</i> permutations.
  */
class KeccakFEquations : public KeccakF {
public:
    /**
      * The constructor. See KeccakF() for more details.
      */
    KeccakFEquations(unsigned int aWidth);
    /**
      * Method that generate equations for the rounds of the chosen Keccak-<i>f</i>
      * instance. The variables starting with letter A are the input of round #0.
      * The B's are the output of round #0 and the input of round #1.
      * The C's are the output of round #1 and the input of round #2, etc.
      *
      * @param  fout    The stream to which the equations are sent.
      * @param  startRoundIndex The index of the first round to perform.
      * @param  nrRounds    The number of rounds to perform.
      * @param  forSage A Boolean value telling whether the syntax of Sage
      *                 should be followed.
      */
    void genRoundEquations(ostream& fout, int startRoundIndex, unsigned nrRounds, bool forSage=false) const;
    /**
      * Method that generate equations for the mappings θ, ρ, π, χ and ι 
      * of the chosen Keccak-<i>f</i> instance.
      *
      * @param  fout    The stream to which the equations are sent.
      * @param  prefixInput     The prefix of the input variables.
      * @param  prefixOutput    The prefix of the output variables.
      */
    void genComponentEquations(ostream& fout, const string& prefixInput, const string& prefixOutput) const;
    /**
      * Method that, given the input to the chosen Keccak-<i>f</i> permutation,
      * generates initialization equations to set the absolute value 
      * just before χ at each round.
      * This can be used in combination with KeccakFDCEquations::buildDCTrailFromPair()
      * and KeccakFDCEquations::genDCEquations() to test that the conditions
      * produced by the latter method are indeed satisfied by a known pair value.
      *
      * @param  fout    The stream to which the equations are sent.
      * @param  input   The input state.
      * @param  prefix  The prefix of the variables.
      * @param  startRoundIndex The index of the first round to perform.
      * @param  nrRounds    The number of rounds to perform.
      */
    void genAbsoluteValuesBeforeChi(ostream& fout, const vector<LaneValue>& input, const string& prefix, int startRoundIndex, unsigned nrRounds) const;
protected:
    /**
      * Internal method to generate the equations from symbolic lanes. 
      * The left hand side of the equations are the variables with prefix @a prefixOutput.
      * The right hand side of the equations are the symbolic bits 
      * in the symbolic lanes in @a state.
      *
      * @param  fout    The stream to which the equations are sent.
      * @param  state   The symbolic input state (right hand side of equations).
      * @param  prefixOutput    The prefix of the variables.
      * @param  forSage A Boolean value telling whether the syntax of Sage
      *                 should be followed.
      */
    void genEquations(ostream& fout, const vector<SymbolicLane>& state, const string& prefixOutput, bool forSage=false) const;
    /**
      * Internal method that initializes the symbolic bits of a symbolic state 
      * with variables using the given prefix.
      *
      * @param  state   The symbolic state to initialize.
      * @param  prefix  The prefix of the variables.
      */
    void initializeState(vector<SymbolicLane>& state, const string& prefix) const;
public:
    /**
      * Method that initializes the symbolic bits of a symbolic state
      * with variables using the given prefix.
      *
      * @param  state   The symbolic state to initialize.
      * @param  prefix  The prefix of the variables.
      * @param  laneSize    The number of bits per lane.
      */
    static void initializeState(vector<SymbolicLane>& state, const string& prefix, unsigned int laneSize);
};

/**
  * Class implementing a symbolic bit in GF(2).
  * The value is a string that can be modified with the methods of this class.
  */
class SymbolicBit
{
public:
    string value;
    bool containsAddition;
public:
    SymbolicBit();
    SymbolicBit(bool aValue);
    SymbolicBit(const string& aValue, bool aContainsAddition=false);
    void complement();
    void add(const SymbolicBit& a);
    void multiply(const SymbolicBit& a);
};

/**
  * Class implementing a vector of symbolic bits to represent a symbolic lane.
  */
class SymbolicLane
{
public:
    vector<SymbolicBit> values;
public:
    SymbolicLane();
    SymbolicLane(LaneValue aValues);
    SymbolicLane(unsigned int laneLength, const string& prefixSymbol);
    void ROL(int offset, unsigned int laneLength);
    friend SymbolicLane operator~(const SymbolicLane& lane);
    friend SymbolicLane operator^(const SymbolicLane& a, LaneValue b);
    friend SymbolicLane operator^(const SymbolicLane& a, const SymbolicLane& b);
    friend SymbolicLane operator&(const SymbolicLane& a, const SymbolicLane& b);
    SymbolicLane& operator^=(const SymbolicLane& b);
};

#endif
