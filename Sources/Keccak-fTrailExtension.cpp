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

#include <sstream>
#include "Keccak-fTrailExtension.h"
#include "translationsymmetry.h"

LowWeightExclusion::LowWeightExclusion()
{
}

void LowWeightExclusion::excludeBelowWeight(unsigned int nrRounds, int weight)
{
    excludedWeight[nrRounds] = weight;
    minWeight.clear();
}

int LowWeightExclusion::getMinWeight(unsigned int nrRounds)
{
    if (nrRounds == 0)
        return 0;
    if (nrRounds > minWeight.size())
        computeExcludedLowWeight(nrRounds);
    return minWeight[nrRounds-1];
}

void LowWeightExclusion::computeExcludedLowWeight(unsigned int upToNrRounds)
{
    minWeight.clear();
    for(unsigned int nrRounds=1; nrRounds<=upToNrRounds; nrRounds++) {
        map<unsigned int, int>::iterator i = excludedWeight.find(nrRounds);
        if (i != excludedWeight.end()) {
            minWeight.push_back(i->second);
        }
        else {
            int max = 0;
            for(unsigned int n1=1; n1<=(nrRounds-1); n1++) {
                unsigned int n2 = nrRounds - n1;
                int sum = minWeight[n1-1] + minWeight[n2-1];
                if (sum > max) max = sum;
            }
            minWeight.push_back(max);
        }
    }
}

ostream& operator<<(ostream& out, const LowWeightExclusion& lwe)
{
    for(unsigned int nrRounds=1; nrRounds<=lwe.minWeight.size(); nrRounds++) {
        map<unsigned int, int>::const_iterator i = lwe.excludedWeight.find(nrRounds);
        bool extrapolated = (i == lwe.excludedWeight.end());
        out.width(2); out.fill(' '); out << dec << nrRounds << " rounds: ";
        out.width(3); out.fill(' '); out << dec << lwe.minWeight[nrRounds-1] << " ";
        if (extrapolated)
            out << "+";
        out << endl;
    }
    return out;
}

KnownSmallWeightStates::KnownSmallWeightStates(const KeccakFPropagation& DCorLC, const string& fileName, int aMaxCompleteWeight)
: maxCompleteWeight(aMaxCompleteWeight)
{
    statesAfterChiPerWeight.resize(maxCompleteWeight+1);
    loadFromFile(DCorLC, fileName);
}

KnownSmallWeightStates::KnownSmallWeightStates(int aMaxCompleteWeight)
: maxCompleteWeight(aMaxCompleteWeight)
{
    statesAfterChiPerWeight.resize(maxCompleteWeight+1);
}

int KnownSmallWeightStates::getMaxCompleteWeight() const
{
    return maxCompleteWeight;
}

void KnownSmallWeightStates::connect(const KeccakFPropagation& DCorLC, const vector<SliceValue>& inputState,
    int maxWeightOut, vector<vector<SliceValue> >& compatibleStates) const
{
    unsigned int inputNrRows = getNrActiveRows(inputState);
    for(unsigned int weight=2; (weight<=(unsigned int)maxWeightOut) && (weight<statesAfterChiPerWeight.size()); weight++)
        for(unsigned int i=0; i<statesAfterChiPerWeight[weight].size(); i++)
            if (getNrActiveRows(statesAfterChiPerWeight[weight][i]) == inputNrRows)
                connect(DCorLC, inputState, statesAfterChiPerWeight[weight][i], compatibleStates);
}

void KnownSmallWeightStates::connect(const KeccakFPropagation& DCorLC, const vector<SliceValue>& inputState,
    const vector<SliceValue>& candidate, vector<vector<SliceValue> >& compatibleStates) const
{
    for(unsigned int z=0; z<DCorLC.laneSize; z++) {
        vector<SliceValue> candidateZ(DCorLC.laneSize);
        for(unsigned int iz=0; iz<DCorLC.laneSize; iz++)
            candidateZ[iz] = candidate[(iz+z)%DCorLC.laneSize];
        if (DCorLC.isChiCompatible(inputState, candidateZ)) {
            vector<SliceValue> candidateZbeforeChi(DCorLC.laneSize);
            DCorLC.directLambda(candidateZ, candidateZbeforeChi);
            compatibleStates.push_back(candidateZbeforeChi);
        }
    }
}

void KnownSmallWeightStates::loadFromFile(const KeccakFPropagation& DCorLC, const string& fileName)
{
    TrailFileIterator fin(fileName, DCorLC);
    for( ; !fin.isEnd(); ++fin) {
        const Trail& trail = *fin;
        for(unsigned int i=(trail.firstStateSpecified ? 0 : 1); i<trail.weights.size(); i++)
            if (trail.weights[i] <= (unsigned int)maxCompleteWeight)
                addState(DCorLC, trail.states[i]);
    }
}

void KnownSmallWeightStates::saveToFile(const KeccakFPropagation& DCorLC, const string& fileName) const
{
    ofstream fout(fileName.c_str());
    for(unsigned int weight=0; weight<statesAfterChiPerWeight.size(); weight++) {
        for(unsigned int i=0; i<statesAfterChiPerWeight[weight].size(); i++) {
            vector<SliceValue> stateBeforeChi;
            DCorLC.directLambda(statesAfterChiPerWeight[weight][i], stateBeforeChi);
            vector<SliceValue> stateBeforeChiMinZ;
            getSymmetricMinimum(stateBeforeChi, stateBeforeChiMinZ);
            Trail trail;
            trail.append(stateBeforeChiMinZ, weight);
            trail.save(fout);
        }
    }
}

void KnownSmallWeightStates::addState(const KeccakFPropagation& DCorLC, const vector<SliceValue>& state)
{
    unsigned int weight = DCorLC.getWeight(state);
    if (weight > (unsigned int)maxCompleteWeight) return;
    vector<SliceValue> stateAfterChi;
    DCorLC.reverseLambda(state, stateAfterChi);
    statesAfterChiPerWeight[weight].push_back(stateAfterChi);
}

KeccakFTrailExtension::KeccakFTrailExtension(const KeccakFDCLC& aParent, KeccakFPropagation::DCorLC aDCorLC)
    : KeccakFPropagation(aParent, aDCorLC),
        showMinimalTrails(false), allPrefixes(false),
        knownSmallWeightStates(0)
{
    knownBounds.excludeBelowWeight(1, 2);
    knownBounds.excludeBelowWeight(2, 8);
    if (parent.getWidth() == 100) {
        if (aDCorLC == KeccakFPropagation::DC) {
            knownBounds.excludeBelowWeight(3, 19);
            knownBounds.excludeBelowWeight(4, 30);
        }
        else {
            knownBounds.excludeBelowWeight(3, 20);
            knownBounds.excludeBelowWeight(4, 38);
        }
    }
    else if (parent.getWidth() == 200) {
        if (aDCorLC == KeccakFPropagation::DC) {
            knownBounds.excludeBelowWeight(3, 20);
            knownBounds.excludeBelowWeight(4, 46);
        }
        else {
            knownBounds.excludeBelowWeight(3, 20);
            knownBounds.excludeBelowWeight(4, 46);
        }
    }
    else if (parent.getWidth() == 1600) {
        if (aDCorLC == KeccakFPropagation::DC) {
            knownBounds.excludeBelowWeight(3, 32);
        }
    }
}

KeccakFTrailExtension::~KeccakFTrailExtension()
{
    if (knownSmallWeightStates)
        delete knownSmallWeightStates;
}

bool KeccakFTrailExtension::isLessThanMinWeightSoFar(unsigned int nrRounds, int weight)
{
    if (nrRounds >= minWeightSoFar.size())
        minWeightSoFar.resize(nrRounds+1, -1);
    if ((minWeightSoFar[nrRounds] < 0) || (weight < minWeightSoFar[nrRounds])) {
        minWeightSoFar[nrRounds] = weight;
        return true;
    }
    else
        return false;
}

void KeccakFTrailExtension::forwardExtendTrails(TrailIterator& trailsIn, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight)
{
    progress.stack("File", trailsIn.getCount());
    for( ; !trailsIn.isEnd(); ++trailsIn) {
        forwardExtendTrail(*trailsIn, trailsOut, nrRounds, maxTotalWeight);
        ++progress;
    }
    progress.unstack();
}

void KeccakFTrailExtension::forwardExtendTrail(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight)
{
    if (trail.stateAfterLastChiSpecified)
        throw KeccakException("KeccakFTrailExtension::forwardExtendTrail() can work only with trail cores or trail prefixes.");
    recurseForwardExtendTrail(trail, trailsOut, nrRounds, maxTotalWeight);
}

void KeccakFTrailExtension::recurseForwardExtendTrail(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight)
{
    int baseWeight = trail.totalWeight;
    int baseNrRounds  = trail.getNumberOfRounds();
    unsigned int curNrRounds = baseNrRounds + 1;
    int curWeight = trail.weights.back();
    int maxWeightOut = maxTotalWeight - baseWeight
        - knownBounds.getMinWeight(nrRounds-baseNrRounds-1);
    if (maxWeightOut < knownBounds.getMinWeight(1))
        return;
    string synopsis;
    {
        stringstream str;
        str << "Weight " << dec << curWeight << " towards round " << dec << curNrRounds;
        str << " (limiting weight to " << dec << maxWeightOut << ")";
        synopsis = str.str();
    }

    const int minWeightInLookingForSmallWeightStates = 16;
    if ((curWeight >= minWeightInLookingForSmallWeightStates) && (knownSmallWeightStates != 0)
            && (maxWeightOut <= knownSmallWeightStates->getMaxCompleteWeight())) {
        vector<vector<SliceValue> > compatibleStates;
        knownSmallWeightStates->connect(*this, trail.states.back(), maxWeightOut, compatibleStates);
        progress.stack(synopsis + " [known small-weight states]", compatibleStates.size());
        for(vector<vector<SliceValue> >::const_iterator i=compatibleStates.begin(); i!=compatibleStates.end(); ++i) {
            int weightOut = getWeight(*i);
            int curWeight = baseWeight + weightOut;
            if (curNrRounds == nrRounds) {
                bool minTrail = showMinimalTrails && isLessThanMinWeightSoFar(curNrRounds, curWeight);
                if (minTrail)
                    cout << "! " << dec << curNrRounds << "-round trail of weight " << dec << curWeight << " found" << endl;
                if ((curWeight <= maxTotalWeight) || minTrail) {
                    Trail newTrail(trail);
                    newTrail.append((*i), weightOut);
                    trailsOut.fetchTrail(newTrail);
                }
            }
            else {
                if (weightOut <= maxWeightOut) {
                    Trail newTrail(trail);
                    newTrail.append((*i), weightOut);
                    recurseForwardExtendTrail(newTrail, trailsOut, nrRounds, maxTotalWeight);
                }
            }
            ++progress;
        }
        progress.unstack();
    }
    else {
        AffineSpaceOfStates base = buildStateBase(trail.states.back());
        SlicesAffineSpaceIterator i(base.originalGenerators, base.offset);
        progress.stack(synopsis + " [affine base]", i.getCount());
        for(; !i.isEnd(); ++i) {
            int weightOut = getWeight(*i);
            int curWeight = baseWeight + weightOut;
            if (curNrRounds == nrRounds) {
                bool minTrail = showMinimalTrails && isLessThanMinWeightSoFar(curNrRounds, curWeight);
                if (minTrail)
                    cout << "! " << dec << curNrRounds << "-round trail of weight " << dec << curWeight << " found" << endl;
                if ((curWeight <= maxTotalWeight) || minTrail) {
                    Trail newTrail(trail);
                    newTrail.append((*i), weightOut);
                    trailsOut.fetchTrail(newTrail);
                }
            }
            else {
                if (weightOut <= maxWeightOut) {
                    Trail newTrail(trail);
                    newTrail.append((*i), weightOut);
                    recurseForwardExtendTrail(newTrail, trailsOut, nrRounds, maxTotalWeight);
                }
            }
            ++progress;
        }
        progress.unstack();
    }
}

void KeccakFTrailExtension::backwardExtendTrails(TrailIterator& trailsIn, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight)
{
    progress.stack("File", trailsIn.getCount());
    for( ; !trailsIn.isEnd(); ++trailsIn) {
        backwardExtendTrail(*trailsIn, trailsOut, nrRounds, maxTotalWeight);
        ++progress;
    }
    progress.unstack();
}

void KeccakFTrailExtension::backwardExtendTrail(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight)
{
    bool isPrefix = trail.firstStateSpecified;
    if (isPrefix) {
        recurseBackwardExtendTrail(trail, trailsOut, nrRounds, maxTotalWeight, true);
    }
    else {
        Trail trimmedTrailPrefix; // cut wrev(a0)
        for(unsigned int i=1; i<trail.states.size(); i++)
            trimmedTrailPrefix.append(trail.states[i], trail.weights[i]);
        recurseBackwardExtendTrail(trimmedTrailPrefix, trailsOut, nrRounds, maxTotalWeight, allPrefixes);
    }
}

void KeccakFTrailExtension::recurseBackwardExtendTrail(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight, bool allPrefixes)
{
    if (!allPrefixes && (nrRounds == (trail.getNumberOfRounds()+1))) {
        int baseWeight = trail.totalWeight;
        vector<SliceValue> stateAfterChi;
        reverseLambda(trail.states[0], stateAfterChi);
        int curMinReverseWeight = getMinReverseWeight(stateAfterChi);
        int curWeight = baseWeight + curMinReverseWeight;
        bool minTrail = showMinimalTrails && isLessThanMinWeightSoFar(nrRounds, curWeight);
        if (minTrail)
            cout << "! " << dec << nrRounds << "-round trail of weight " << dec << curWeight << " found" << endl;
        if ((curWeight <= maxTotalWeight) || minTrail) {
            Trail newTrail;
            newTrail.setFirstStateReverseMinimumWeight(curMinReverseWeight);
            newTrail.append(trail);
            trailsOut.fetchTrail(newTrail);
        }
    }
    else {
        int baseWeight = trail.totalWeight;
        int baseNrRounds  = trail.getNumberOfRounds();
        int maxWeightOut = maxTotalWeight - baseWeight
            - knownBounds.getMinWeight(nrRounds-baseNrRounds-1);
        if (maxWeightOut < knownBounds.getMinWeight(1))
            return;
        vector<SliceValue> stateAfterChi;
        reverseLambda(trail.states[0], stateAfterChi);
        ReverseStateIterator i(stateAfterChi, *this, maxWeightOut);
        if (i.isEmpty())
            return;
        unsigned int curNrRounds = baseNrRounds + 1;
        {
            stringstream str;
            str << dec << getNrActiveRows(stateAfterChi) << " active rows towards round -" << dec << curNrRounds;
            str << " (limiting weight to " << dec << maxWeightOut << ")";
            progress.stack(str.str());
        }
        for(; !i.isEnd(); ++i) {
            int weightOut = getWeight(*i);
            int curWeight = baseWeight + weightOut;
            if (curNrRounds == nrRounds) {
                bool minTrail = showMinimalTrails && isLessThanMinWeightSoFar(nrRounds, curWeight);
                if (minTrail)
                    cout << "! " << dec << nrRounds << "-round trail of weight " << dec << curWeight << " found" << endl;
                if ((curWeight <= maxTotalWeight) || minTrail) {
                    Trail newTrail(trail);
                    newTrail.prepend((*i), weightOut);
                    trailsOut.fetchTrail(newTrail);
                }
            }
            else {
                int minPrevWeight = getMinReverseWeightAfterLambda(*i);
                if ((curWeight + minPrevWeight + knownBounds.getMinWeight(nrRounds-curNrRounds-1)) <= maxTotalWeight) {
                    Trail newTrail(trail);
                    newTrail.prepend((*i), weightOut);
                    recurseBackwardExtendTrail(newTrail, trailsOut, nrRounds, maxTotalWeight, allPrefixes);
                }
            }
            ++progress;
        }
        progress.unstack();
    }
}

