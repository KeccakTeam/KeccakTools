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

#ifndef _KECCAKFPARITYBOUNDS_H_
#define _KECCAKFPARITYBOUNDS_H_

#include <vector>
#include "Keccak-fParts.h"
#include "Keccak-fPropagation.h"

/** Given the parity @a C and the θ-effect @a D, this function computes
  * a lower bound on the total number of active rows before and after λ.
  * This follows the algorithm in the paper "Differential propagation of Keccak"
  * at Fast Software Encryption 2012.
  * @param   DCorLC The propagation context of the trail,
  *                 as a reference to a KeccakFPropagation object.
  * @param  C   The parity as a vector or row values.
  * @param  D   The θ-effect as a vector or row values.
  * @return A lower bound on the number of total active rows a state
  *     has before and after λ if it has the given parity and θ-effect.
  */
unsigned int getLowerBoundTotalActiveRows(const KeccakFPropagation& DCorLC,
    const vector<RowValue>& C, const vector<RowValue>& D);

/** Given the total Hamming weight of a state before and after λ, this
  * function computes a lower bound on the total weight
  * before and after λ, i.e., the minimum reverse weight before λ
  * plus the propagation weight after λ.
  * @param   DCorLC The propagation context of the trail,
  *                 as a reference to a KeccakFPropagation object.
  * @param  totalHW The total Hamming weight before and after λ.
  * @return A lower bound on the minimum reverse weight before λ
  *     plus the propagation weight after λ of any state having the given
  *     total Hamming weight.
  */
unsigned int getBoundOfTotalWeightGivenTotalHammingWeight(const KeccakFPropagation& DCorLC, unsigned int totalHW);

/** Class representing a parity run, with its starting coordinate and its length.
  */
class Run
{
public:
    /** The @a t coordinate of the starting point. */
    unsigned int tStart;
    /** The length of the run. */
    unsigned int length;
public:
    /** This constructor creates an empty run. */
    Run() : tStart(0), length(0) {}
    /** This constructor creates a run.
      * @param  aTStart The @a t coordinate of the starting point.
      * @param  aLength The length of the run.
      */
    Run(unsigned int aTStart, unsigned int aLength) : tStart(aTStart), length(aLength) {}
    /** Function that returns a string describing the run. */
    string display() const;
};

/** Class representing the parity of a state as a collection of distinct runs.
  */
class ParityAsRuns
{
public:
    /** The collection of runs. */
    std::vector<Run> runs;
public:
    /** This constructor initializes to the all-zero parity (no runs). */
    ParityAsRuns() {}
    /** Method to convert the runs into a parity expressed as vector of
      * row values and its associated θ-effect.
      * @param   DCorLC The propagation context ,
      *                 as a reference to a KeccakFPropagation object.
      * @param  C   The parity to be stored as a vector or row values.
      * @param  D   The θ-effect to be stored as a vector or row values.
      */
    void toParityAndParityEffect(const KeccakFPropagation& DCorLC, vector<RowValue>& C, vector<RowValue>& D) const;
    /** Function that returns a string describing the runs. */
    string display() const;
    /** This method returns a lower bound
      * on the total Hamming weight before and after λ
      * for any state having this parity.
      * This is 10 per affected column and 2 per unaffected odd column.
      * @param   DCorLC The propagation context ,
      *                 as a reference to a KeccakFPropagation object.
      *  @return    The lower bound on the total Hamming weight.
      */
    unsigned int getLowerBoundTotalHammingWeight(const KeccakFPropagation& DCorLC) const;
    /** This method returns a lower bound
      * on the total number of active rows before and after λ
      * for any state having this parity.
      * See also ::getLowerBoundTotalActiveRows().
      * @param   DCorLC The propagation context ,
      *                 as a reference to a KeccakFPropagation object.
      *  @return    The lower bound on the total number of active rows.
      */
    unsigned int getLowerBoundTotalActiveRows(const KeccakFPropagation& DCorLC) const;
    /** This method is like getLowerBoundTotalActiveRows(),
      * except that the bound does not take into account the contribution
      * of unaffected odd columns.
      */
    unsigned int getLowerBoundTotalActiveRowsUsingOnlyAC(const KeccakFPropagation& DCorLC) const;
};

/** This function looks for all parities (up to translation in z) such that
  * the lower bound provided by ::getLowerBoundTotalActiveRows() is not
  * higher than a given target.
  * The target is given as total weight
  * (the minimum reverse weight before λ
  * plus the propagation weight after λ), which is lower bounded
  * by 2 times the total number of active rows.
  * The search is done as explained in the paper "Differential propagation of Keccak"
  * at Fast Software Encryption 2012.
  * @param   DCorLC The propagation context ,
  *                 as a reference to a KeccakFPropagation object.
  * @param  out The output file where to store the found parities.
  * @param  targetWeight    The target total weight.
  * @param  verbose If true, the funtion displays information on the standard output.
  */
void lookForRunsBelowTargetWeight(const KeccakFPropagation& DCorLC, ostream& out,  unsigned int targetWeight, bool verbose = false);

#endif
