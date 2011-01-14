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

#ifndef _KECCAKFTRAILS_H_
#define _KECCAKFTRAILS_H_

#include <fstream>
#include <iostream>
#include "Keccak-fParts.h"

class KeccakFPropagation;

/**
 * Exception that can be thrown by the class Trail.
 */
class TrailException {
public:
    /** A string expressing the reason for the exception. */
    string reason;
    /**
     * The constructor.
     */
    TrailException() : reason() {}
    /**
     * The constructor.
     * @param   aReason A string giving the reason of the exception.
     */
    TrailException(const string& aReason) : reason(aReason) {}
};

/** This class implements a container for a differential or linear trail. 
  * A trail makes sense only in the context of a KeccakFPropagation object.
  * The main attribute of the class is the sequence of state values before χ. 
  * If S<sub>i</sub> = @a states[i], the trail is: 
  * S<sub>0</sub> χλ S<sub>1</sub> χλ S<sub>2</sub> χλ ... 
  * χλ S<sub><i>n</i>-1</sub>,
  * with <i>n</i> = @a states.size().
  */
class Trail {
public:
    /** This attribute contains the list of states round after round, before χ.
      */
    vector<vector<SliceValue> > states;
    /** This attribute contains the propagation weights of the states 
      * in @a states. So, @a weights has the same size 
      * as @a states.
      */
    vector<unsigned int> weights;
    /** This attribute contains the sum of the weights[i].
      */
    unsigned int totalWeight;
public:
    /** This constructor creates an empty trail.
    */
    Trail();
    /** This constructor loads a trail from an input stream.
      * @param   fin    The input stream to read the trail from.
      */
    Trail(istream& fin);
    /** This constructor initializes a trail by copying the
      * trail given in parameter.
      * @param   other  The original trail to copy.
      */
    Trail(const Trail& other)
        : states(other.states), weights(other.weights),
        totalWeight(other.totalWeight) {}
    /** This method appends a state to the end of @a states, 
      * with its corresponding propagation weight.
      * @param   state  The state to add.
      * @param   weight The propagation weight.
      */
    void append(const vector<SliceValue>& state, unsigned int weight);
    /** This method appends another trail to the current trail.
      * @param   otherTrail The trail to append.
      */
    void append(const Trail& otherTrail);
    /** This method inserts a state at the beginning of @a states, 
      * with its corresponding propagation weight.
      * @param   state  The state to add.
      * @param   weight The propagation weight.
      */
    void prepend(const vector<SliceValue>& state, unsigned int weight);
    /** This method displays the trail for in a human-readable form.
      * @param   DCorLC The propagation context of the trail, 
      *                 as a reference to a KeccakFPropagation object.
      * @param  fout    The stream to display to.
      */
    void display(const KeccakFPropagation& DCorLC, ostream& fout) const;
    /** This methods loads the trail from a stream (e.g., file).
      * @param   fin    The input stream to read the trail from.
      */
    void load(istream& fin);
    /** This methods outputs the trail to save it in, e.g., a file.
      * @param  fout    The stream to save the trail to.
      */
    void save(ostream& fout) const;
    /** This function reads all the trails in a file, checks their consistency
      * and then produces a report.
      * The report is output in a file with the same file name plus ".txt".
      * See also KeccakFPropagation::displayTrailsAndCheck().
      * @param   DCorLC     The propagation context of the trails, 
      *                     as a reference to a KeccakFPropagation object.
      * @param   fileName   The name of the file containing the trails.
      * @param   verbose    If true, the function will display the name of
      *                     the file written to cout.
      * @param   maxWeight  As in KeccakFPropagation::displayTrailsAndCheck().
      * @return The number of trails read and checked.
      */
    static UINT64 produceHumanReadableFile(const KeccakFPropagation& DCorLC, const string& fileName,
        bool verbose = true, unsigned int maxWeight = 0);
};

#endif
