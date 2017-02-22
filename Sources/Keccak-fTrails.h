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

typedef Exception TrailException;

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
    /** This attribute tells whether the first state (states[0]) is specified.
      *  If false, the trail is actually a trail core. In this case, weights[0]
      * is the minimum reverse weight of λ<sup>-1</sup>(states[1]).
      */
    bool firstStateSpecified;
    /** This attribute contains the list of states round after round, before χ.
      */
    vector<vector<SliceValue> > states;
    /** This attribute tells whether the state after the last χ is specified.
      *  If true, the state after the last χ is specified in stateAfterLastChi.
      */
    bool stateAfterLastChiSpecified;
    /**
      *  If stateAfterLastChiSpecified==true, this attribute contains the state after the last χ.
      */
    vector<SliceValue> stateAfterLastChi;
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
        : firstStateSpecified(other.firstStateSpecified),
        states(other.states),
        stateAfterLastChiSpecified(other.stateAfterLastChiSpecified),
        stateAfterLastChi(other.stateAfterLastChi),
        weights(other.weights),
        totalWeight(other.totalWeight) {}
    /** This method returns the number of rounds the trail represents.
      *  @return    The number of rounds.
      */
    unsigned int getNumberOfRounds() const;
    /** This method empties the trail.
     */
    void clear();
    /** For a trail core, this sets the minimum reverse weight
      * of the first state.
      * @param   weight The minimum reverse weight.
      */
    void setFirstStateReverseMinimumWeight(unsigned int weight);
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

/** This base class represents a filter on trails, to be used with the class TrailIterator
 * or descendants.
 */
class TrailFilter
{
public:
    /** This method tells whether to keep (true) or discard (false) the given trail.
      * @param   DCorLC     The propagation context of the trails,
      *                     as a reference to a KeccakFPropagation object.
      * @param  trail   The trail to consider.
      * @return The result of the filter: whether to keep (true) or discard (false) the given trail.
      */
    virtual bool filter(const KeccakFPropagation& DCorLC, const Trail& trail) const = 0;
};

/** This class implements a filter on trails, combining one or more given filters
 * with a logical 'and'.
 */
class TrailFilterAND : public TrailFilter
{
public:
    /** The set of filters. */
    vector<TrailFilter *> filters;
public:
    /** The default constructor. No filter is added to filters be default.*/
    TrailFilterAND();
    /**  The constructor specifying two filters to combine.
     *  @param  filter1 A pointer to the first filter.
     *  @param  filter2 A pointer to the second filter.
     */
    TrailFilterAND(TrailFilter *filter1, TrailFilter *filter2);
    /** This method returns the logical 'and' of the filters in the filters attribute.
     * See TrailFilter::filter() for more details.
     */
    virtual bool filter(const KeccakFPropagation& DCorLC, const Trail& trail) const;
};

/** This base class represents an iterator on a set of trails.
 */
class TrailIterator {
protected:
    const KeccakFPropagation& DCorLC;
    TrailFilterAND filterAND;
    TrailFilter *filter;
public:
    /** The default constructor.
      * @param   DCorLC     The propagation context of the trails,
      *                     as a reference to a KeccakFPropagation object.
      */
    TrailIterator(const KeccakFPropagation& aDCorLC)
        : DCorLC(aDCorLC), filter(0) {}
    /** The constructor with one filter.
      * @param   DCorLC     The propagation context of the trails,
      *                     as a reference to a KeccakFPropagation object.
      * @param  aFilter     A pointer to the filter.
      */
    TrailIterator(const KeccakFPropagation& aDCorLC, TrailFilter *aFilter)
        : DCorLC(aDCorLC), filter(aFilter) {}
    /** The constructor with two filters.
      * @param   DCorLC     The propagation context of the trails,
      *                     as a reference to a KeccakFPropagation object.
      * @param  aFilter1    A pointer to the first filter.
      * @param  aFilter2    A pointer to the second filter.
      */
    TrailIterator(const KeccakFPropagation& aDCorLC, TrailFilter *aFilter1, TrailFilter *aFilter2)
        : DCorLC(aDCorLC), filterAND(aFilter1, aFilter2), filter(&filterAND) {}
    /** This method indicates whether the iterator has reached the end of the set of trails.
     *  @return True if there are no more trails to consider.
     */
    virtual bool isEnd() = 0;
    /** This method indicates whether the set of trails is empty.
     *  @return True if the set is empty.
     */
    virtual bool isEmpty() = 0;
    /** This method indicates whether the iterator knows the cardinality of the set.
     *  @return True if the cardinality is known.
     */
    virtual bool isBounded() = 0;
    /** This method returns the index of the current trail encountered,
     * starting from 0 for the first trail.
     *  @return The index of the current trail.
     */
    virtual UINT64 getIndex() = 0;
    /** If the cardinality of the set is known, this method returns it.
     *  @return The cardinality of the set.
     */
    virtual UINT64 getCount() = 0;
    /** The '++' operator moves the iterator to the next trail in the set.
     */
    virtual void operator++() = 0;
    /** The '*' operator gives a constant reference to the current trail.
     *  @return A constant reference to the current trail.
     */
    virtual const Trail& operator*() = 0;
};

/** This class implements an iterator on a set of trails read from a file.
 */
class TrailFileIterator : public TrailIterator {
protected:
    ifstream fin;
    string fileName;
    bool prefetch;
    UINT64 i, count, unfilteredCount;
    bool end;
    Trail current;
public:
    /** The constructor of the iterator.
      * @param  aFileName   The name of the file to read from.
      * @param   DCorLC     The propagation context of the trails,
      *                     as a reference to a KeccakFPropagation object.
      * @param  aPrefetch   Whether the file has to be first read to determine
      *                     the number of trails in the file.
      */
    TrailFileIterator(const string& aFileName, const KeccakFPropagation& aDCorLC,
        bool aPrefetch = true);
    /** The constructor of the iterator, with one filter.
      * @param  aFileName   The name of the file to read from.
      * @param   DCorLC     The propagation context of the trails,
      *                     as a reference to a KeccakFPropagation object.
      * @param  aFilter A pointer to the filter.
      * @param  aPrefetch   Whether the file has to be first read to determine
      *                     the number of trails in the file.
      */
    TrailFileIterator(const string& aFileName, const KeccakFPropagation& aDCorLC,
        TrailFilter *aFilter, bool aPrefetch = true);
    /** The constructor of the iterator, with two filters.
      * @param  aFileName   The name of the file to read from.
      * @param   DCorLC     The propagation context of the trails,
      *                     as a reference to a KeccakFPropagation object.
      * @param  aFilter1    A pointer to the first filter.
      * @param  aFilter2    A pointer to the second filter.
      * @param  aPrefetch   Whether the file has to be first read to determine
      *                     the number of trails in the file.
      */
    TrailFileIterator(const string& aFileName, const KeccakFPropagation& aDCorLC,
        TrailFilter *aFilter1, TrailFilter *aFilter2, bool aPrefetch = true);
    /** This method displays information about the file being read.
      * @param  fout    The stream to display to.
      */
    void display(ostream& fout) const;
    /** See TrailIterator::isEnd(). */
    virtual bool isEnd();
    /** See TrailIterator::isEmpty(). */
    virtual bool isEmpty();
    /** See TrailIterator::operator++(). */
    virtual void operator++();
    /** See TrailIterator::operator*(). */
    virtual const Trail& operator*();
    /** See TrailIterator::isBounded(). */
    virtual bool isBounded();
    /** See TrailIterator::getIndex(). */
    virtual UINT64 getIndex();
    /** See TrailIterator::getCount(). */
    virtual UINT64 getCount();
    /** This method returns the number of trails in the file before applying the filter(s).
      * @return The number of trails in the file before applying the filter(s).
      */
    UINT64 getUnfilteredCount() const;
    friend ostream& operator<<(ostream& a, const TrailFileIterator& tfi);
protected:
    void initialize();
    void next();
};

/** This base class represents the output of trails, which can be
  * for instance saved or further processed.
  */
class TrailFetcher {
public:
    /** Method to call when a trail is produced.
      * @param  trail   The trail to save or to process.
      */
    virtual void fetchTrail(const Trail& trail) = 0;
};

/** This class implements a TrailFetcher and saves the trails in a file.
  */
class TrailSaveToFile : public TrailFetcher {
protected:
    ostream& fout;
public:
    /** The constructor.
      * @param  aFout   The output stream to save the trails to.
      */
    TrailSaveToFile(ostream& aFout);
    /** See TrailFetcher::fetchTrail().*/
    void fetchTrail(const Trail& trail);
};

#endif
