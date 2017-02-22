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

#ifndef _KECCAKFTRAILEXTENSION_H_
#define _KECCAKFTRAILEXTENSION_H_

#include <map>
#include <vector>
#include "Keccak-fPropagation.h"
#include "progress.h"

using namespace std;

/** Class that maintains a list of minimum weights for each number of rounds
  * below which there is no need to look for trails.
  * This is typically used during trail extension to avoid looking for trails
  * with weight lower than the lower bound, or to exclude a subspace
  * already covered.
  */
class LowWeightExclusion {
protected:
    /** The explicitly excluded weights per number of rounds. */
    map<unsigned int, int> excludedWeight;
    /** The interpolated minimum weights per number of rounds.
      * Note that minWeight[nrRounds-1] contains the minimum weight
      * for nrRounds rounds.
      */
    vector<int> minWeight;
public:
    /** The constructor. */
    LowWeightExclusion();
    /** This method tells to exclude trails with weight below the given value
      * for the given number of rounds.
      * @param  nrRounds    The number of rounds.
      * @param  weight  The weight below which trails are to be excluded.
      */
    void excludeBelowWeight(unsigned int nrRounds, int weight);
    /** For a given number of rounds, this function returns the minimum weight
      * to consider.
      * @param  nrRounds    The number of rounds.
      * @return The minimum weight to consider.
      */
    int getMinWeight(unsigned int nrRounds);
    friend ostream& operator<<(ostream& out, const LowWeightExclusion& lwe);
protected:
    void computeExcludedLowWeight(unsigned int upToNrRounds);
};

/** Class that contains all states C after χ such that D=λ(C) has low weight.
  * This allows to quickly check whether a state B can be connected
  * to such a state C through χ.
  */
class KnownSmallWeightStates {
protected:
    /** Vector that contains all states C per weight of D=λ(C). The first index level
      * is for the weight of D=λ(C). All states are stored up to translation in z.
      */
    vector<vector<vector<SliceValue> > > statesAfterChiPerWeight;
    /** This attribute tells up to which weight the set of states in
     * statesAfterChiPerWeight is complete.
     */
    int maxCompleteWeight;
public:
    /** Constructor that initializes an empty set of states.
      * @param  aMaxCompleteWeight  The intended weight up to which the set is complete.
      */
    KnownSmallWeightStates(int aMaxCompleteWeight);
    /** Constructor that initializes a set of states from a file.
      * @param  DCorLC The propagation context,
      *                 as a reference to a KeccakFPropagation object.
      * @param  fileName    The name of the file to extract states from.
      *     See loadFromFile() for more information.
      * @param  aMaxCompleteWeight  The intended weight up to which the set is complete.
      */
    KnownSmallWeightStates(const KeccakFPropagation& DCorLC, const string& fileName, int aMaxCompleteWeight);
    /** Method that loads states from a file.
      * @param  DCorLC The propagation context,
      *                 as a reference to a KeccakFPropagation object.
      * @param  fileName    The name of the file to extract states from.
      *     The file has to contain trails (or trail cores or trail prefixes)
      *     and all states with weight less than maxCompleteWeight will
      *     be fetched.
      *     The states fetched are interpreted as D=λ(C).
      */
    void loadFromFile(const KeccakFPropagation& DCorLC, const string& fileName);
    /** Method that returns maxCompleteWeight.
      * @return The intended weight up to which the set is complete.
      */
    int getMaxCompleteWeight() const;
    /** Method that attempts at finding states D that are compatible
      * through λ after χ with the given state B.
      * @param  DCorLC The propagation context,
      *                 as a reference to a KeccakFPropagation object.
      * @param  inputState  The input state B.
      * @param  maxWeightOut    The maximum weight of D that needs
      *     to be considered.
      * @param  compatibleStates    The list of states D that are compatible
      *     with B.
      */
    void connect(const KeccakFPropagation& DCorLC, const vector<SliceValue>& inputState,
        int maxWeightOut, vector<vector<SliceValue> >& compatibleStates) const;
    /** Method that stores all the states in a file.
      * The states are saved as D=λ(C) in 1-round trail prefixes.
      * @param  DCorLC The propagation context,
      *                 as a reference to a KeccakFPropagation object.
      * @param  fileName    The name of the file to save to.
      */
    void saveToFile(const KeccakFPropagation& DCorLC, const string& fileName) const;
protected:
    void addState(const KeccakFPropagation& DCorLC, const vector<SliceValue>& state);
    void connect(const KeccakFPropagation& DCorLC, const vector<SliceValue>& inputState,
        const vector<SliceValue>& candidate, vector<vector<SliceValue> >& compatibleStates) const;
};

/** This class provides trail extension services.
  */
class KeccakFTrailExtension : public KeccakFPropagation
{
public:
    /** If true, the best trails found will be output, even if their weight
      * is higher than the requested maximum.
      */
    bool showMinimalTrails;
    /** If true, the backward extension of trail cores looks for all prefixes,
      * not just trail cores.
      */
    bool allPrefixes;
    /** This LowWeightExclusion object specifies search areas to exclude,
      * primarily because bounds are known. For instance, this expresses
      * that any 2-round trail in Keccak-f has at least weight 8,
      * and this knowledge can limit the weight of intermediate states
      * that are to be considered.
      */
    LowWeightExclusion knownBounds;
    /** This optional KnownSmallWeightStates object pointer provides
      * a list of states with low weight. This is used to optimize forward
      * trail extension.
      */
    KnownSmallWeightStates *knownSmallWeightStates;
protected:
    vector<int> minWeightSoFar;
    ProgressMeter progress;
public:
    /** The constructor. See KeccakFPropagation::KeccakFPropagation(). */
    KeccakFTrailExtension(const KeccakFDCLC& aParent, KeccakFPropagation::DCorLC aDCorLC);
    /** The destructor.
      * This frees the memory taken by @a knownSmallWeightStates.
      */
    virtual ~KeccakFTrailExtension();
    /** Starting from a given trail (prefix or core), this method
      * appends states to it, and systematically looks
      * for all trails with @a nrRounds rounds
      * up to total weight @a maxTotalWeight
      * that have the given trail as prefix.
      * @param  trail   The starting trail core or trail prefix.
      * @param  trailsOut   Where to output the found trails.
      * @param  nrRounds    The target number of rounds.
      * @param  maxTotalWeight  The maximum total weight to consider.
      */
    void forwardExtendTrail(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight);
    /** This function is like forwardExtendTrail(), except that it processes
      * all the trails from @a trailsIn.
      * @param  trailsIn    The starting trail cores or trail prefixes.
      * @param  trailsOut   Where to output the found trails.
      * @param  nrRounds    The target number of rounds.
      * @param  maxTotalWeight  The maximum total weight to consider.
      */
    void forwardExtendTrails(TrailIterator& trailsIn, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight);
    /** Starting from a given trail (prefix or core), this method
      * prepends states to it, and systematically looks
      * for all trails with @a nrRounds rounds
      * up to total weight @a maxTotalWeight
      * that have the given trail as suffix.
      * If the given trail is a trail core, then the function looks for all trail cores,
      * unless @a allPrefixes is set to true, in which case all trail prefixes
      * are output.
      * Otherwise if the given trail is a trail prefix, then the function looks for all trail prefixes.
      * @param  trail   The starting trail core or trail prefix.
      * @param  trailsOut   Where to output the found trails.
      * @param  nrRounds    The target number of rounds.
      * @param  maxTotalWeight  The maximum total weight to consider.
      */
    void backwardExtendTrail(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight);
    /** This function is like backwardExtendTrails(), except that it processes
      * all the trails from @a trailsIn.
      * @param  trailsIn    The starting trail cores or trail prefixes.
      * @param  trailsOut   Where to output the found trails.
      * @param  nrRounds    The target number of rounds.
      * @param  maxTotalWeight  The maximum total weight to consider.
      */
    void backwardExtendTrails(TrailIterator& trailsIn, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight);



protected:
    void recurseForwardExtendTrail(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight);
    void recurseBackwardExtendTrail(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight, bool allPrefixes);
    bool isLessThanMinWeightSoFar(unsigned int nrRounds, int weight);
};

#endif
