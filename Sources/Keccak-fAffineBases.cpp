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

#include <algorithm>
#include <iostream>
#include <math.h>
#include <set>
#include "Keccak-fAffineBases.h"
#include "Keccak-fDisplay.h"

using namespace std;


// -------------------------------------------------------------
//
// AffineSpaceOfRows
//
// -------------------------------------------------------------

void AffineSpaceOfRows::display(ostream& fout) const
{
    fout.fill('0'); fout.width(2); fout << hex << (int)offset;
    fout << " + <";
    for(unsigned int i=0; i<generators.size(); i++) {
        fout.fill('0'); fout.width(2); fout << hex << (int)generators[i];
        if (i < (generators.size()-1))
            fout << ", ";
    }
    fout << ">" << endl;
}

// -------------------------------------------------------------
//
// AffineSpaceOfSlices
//
// -------------------------------------------------------------

AffineSpaceOfSlices::AffineSpaceOfSlices(vector<SliceValue>& aGenerators, vector<RowValue>& aGeneratorParities, SliceValue aOffset, RowValue aOffsetParity)
    : offset(aOffset), 
    offsetParity(aOffsetParity)
{
    setGenerators(aGenerators, aGeneratorParities);
}

void AffineSpaceOfSlices::setGenerators(vector<SliceValue>& aGenerators, vector<RowValue>& aGeneratorParities)
{
    // Copy the generators into originalGenerators;
    originalGenerators = aGenerators;

    // Upper-triangularize the parities
    for(unsigned int x=0; x<nrRowsAndColumns; x++) {
        // Look for a generator with a parity 1 at position x
        RowValue selectX = 1 << x;
        bool found = false;
        SliceValue foundSlice;
        RowValue foundParity;
        for(unsigned int i=0; i<aGenerators.size(); i++) {
            if ((aGeneratorParities[i] & selectX) != 0) {
                foundSlice = aGenerators[i];
                offsetGenerators.push_back(foundSlice);
                foundParity = aGeneratorParities[i];
                offsetParities.push_back(foundParity);
                found = true;
                break;
            }
        }
        // If found, cancel all parities at position x for the other aGenerators
        if (found) {
            for(unsigned int i=0; i<aGenerators.size(); i++) {
                if ((aGeneratorParities[i] & selectX) != 0) {
                    aGenerators[i] ^= foundSlice;
                    aGeneratorParities[i] ^= foundParity;
                }
            }
        }
    }
    // The remaining aGenerators have zero parity
    for(unsigned int i=0; i<aGenerators.size(); i++)
        if (aGenerators[i] != 0)
            kernelGenerators.push_back(aGenerators[i]);
}

void AffineSpaceOfSlices::display(ostream& fout) const
{
    fout << "Offset = " << endl;
    displaySlice(fout, offset);
    fout << endl;

    if (originalGenerators.size() == 0)
        fout << "No generators" << endl;
    else {
        fout << dec << originalGenerators.size() << " generators:" << endl;
        for(unsigned int i=0; i<originalGenerators.size(); i++) {
            displaySlice(fout, originalGenerators[i]);
            fout << endl;
        }
        if (offsetGenerators.size() == 0)
            fout << "No parity-offset generators" << endl;
        else {
            fout << dec << offsetGenerators.size() << " parity-offset generators:" << endl;
            for(unsigned int i=0; i<offsetGenerators.size(); i++) {
                displaySlice(fout, offsetGenerators[i]);
                fout << endl;
            }
        }
        if (kernelGenerators.size() == 0)
            fout << "No parity-kernel generators" << endl;
        else {
            fout << dec << kernelGenerators.size() << " parity-kernel generators:" << endl;
            for(unsigned int i=0; i<kernelGenerators.size(); i++) {
                displaySlice(fout, kernelGenerators[i]);
                fout << endl;
            }
        }
    }
}

bool AffineSpaceOfSlices::getOffsetWithGivenParity(RowValue parity, SliceValue& output) const
{
    output = offset;
    RowValue outputParity = offsetParity;
    RowValue correctionParity = parity^offsetParity;

    unsigned int i = 0;

    for(unsigned int x=0; x<nrRowsAndColumns; x++) {
        RowValue mask = (1<<(x+1))-1;
        if ((correctionParity & (1<<x)) != 0) {
            while((i<offsetParities.size()) && ((offsetParities[i] & mask) != (1<<x)))
                i++;
            if (i<offsetParities.size()) {
                output ^= offsetGenerators[i];
                correctionParity ^= offsetParities[i];
            }
            else
                return false;
        }
    }
    return (correctionParity == 0);
}

// -------------------------------------------------------------
//
// AffineSpaceOfStates
//
// -------------------------------------------------------------

AffineSpaceOfStates::AffineSpaceOfStates(unsigned int aLaneSize, vector<vector<SliceValue> >& aGenerators, vector<PackedParity>& aGeneratorParities, const vector<SliceValue>& aOffset, PackedParity aOffsetParity)
    : laneSize(aLaneSize), offset(aOffset), offsetParityPacked(aOffsetParity), packed(true)
{
    setGenerators(aGenerators, aGeneratorParities);
}

AffineSpaceOfStates::AffineSpaceOfStates(unsigned int aLaneSize, vector<vector<SliceValue> >& aGenerators, vector<vector<RowValue> >& aGeneratorParities, const vector<SliceValue>& aOffset, const vector<RowValue>& aOffsetParity)
    : laneSize(aLaneSize), offset(aOffset), offsetParity(aOffsetParity), packed(false)
{
    setGenerators(aGenerators, aGeneratorParities);
}

void AffineSpaceOfStates::setGenerators(vector<vector<SliceValue> >& aGenerators, vector<PackedParity>& aGeneratorParities)
{
    // Copy the generators into originalGenerators;
    originalGenerators = aGenerators;
    for(unsigned int i=0; i<aGeneratorParities.size(); i++) {
        vector<RowValue> parities;
        for(unsigned int z=0; z<laneSize; z++)
            parities.push_back(getParityFromPackedParity(aGeneratorParities[i], z));
        originalParities.push_back(parities);
    }

    // Upper-triangularize the parities
    for(unsigned int xz=0; xz<(nrRowsAndColumns*laneSize); xz++) {
        // Look for a generator with a parity 1 at position xz
        PackedParity selectXZ = (PackedParity)1 << xz;
        bool found = false;
        vector<SliceValue> foundState;
        PackedParity foundParity;
        for(unsigned int i=0; i<aGenerators.size(); i++) {
            if ((aGeneratorParities[i] & selectXZ) != 0) {
                foundState = aGenerators[i];
                offsetGenerators.push_back(foundState);
                foundParity = aGeneratorParities[i];
                offsetParitiesPacked.push_back(foundParity);
                found = true;
                break;
            }
        }
        // If found, cancel all parities at position x for the other aGenerators
        if (found) {
            for(unsigned int i=0; i<aGenerators.size(); i++) {
                if ((aGeneratorParities[i] & selectXZ) != 0) {
                    for(unsigned int z=0; z<laneSize; z++)
                        aGenerators[i][z] ^= foundState[z];
                    aGeneratorParities[i] ^= foundParity;
                }
            }
        }
    }
    // The remaining aGenerators have zero parity
    for(unsigned int i=0; i<aGenerators.size(); i++) {
        bool zero = true;
        for(unsigned int j=0; j<aGenerators[i].size(); j++)
            if (aGenerators[i][j] != 0) {
                zero = false;
                break;
            }
        if (!zero)
            kernelGenerators.push_back(aGenerators[i]);
    }
}

void AffineSpaceOfStates::setGenerators(vector<vector<SliceValue> >& aGenerators, vector<vector<RowValue> >& aGeneratorParities)
{
    // Copy the generators into originalGenerators;
    originalGenerators = aGenerators;
    originalParities = aGeneratorParities;

    // Upper-triangularize the parities
    for(unsigned int z=0; z<laneSize; z++)
    for(unsigned int x=0; x<nrRowsAndColumns; x++) {
        // Look for a generator with a parity 1 at position (x,z)
        RowValue selectX = 1 << x;
        bool found = false;
        vector<SliceValue> foundState;
        vector<RowValue> foundParity;
        for(unsigned int i=0; i<aGenerators.size(); i++) {
            if ((aGeneratorParities[i][z] & selectX) != 0) {
                foundState = aGenerators[i];
                offsetGenerators.push_back(foundState);
                foundParity = aGeneratorParities[i];
                offsetParities.push_back(foundParity);
                found = true;
                break;
            }
        }
        // If found, cancel all parities at position (x,z) for the other aGenerators
        if (found) {
            for(unsigned int i=0; i<aGenerators.size(); i++) {
                if ((aGeneratorParities[i][z] & selectX) != 0) {
                    for(unsigned int jz=0; jz<laneSize; jz++)
                        aGenerators[i][jz] ^= foundState[jz];
                    for(unsigned int jz=0; jz<laneSize; jz++)
                        aGeneratorParities[i][jz] ^= foundParity[jz];
                }
            }
        }
    }
    // The remaining aGenerators have zero parity
    for(unsigned int i=0; i<aGenerators.size(); i++) {
        bool zero = true;
        for(unsigned int j=0; j<aGenerators[i].size(); j++)
            if (aGenerators[i][j] != 0) {
                zero = false;
                break;
            }
        if (!zero)
            kernelGenerators.push_back(aGenerators[i]);
    }
}

void AffineSpaceOfStates::display(ostream& fout) const
{
    fout << "Offset = " << endl;
    displayState(fout, offset);
    fout << endl;

    if (originalGenerators.size() == 0)
        fout << "No generators" << endl;
    else {
        fout << dec << originalGenerators.size() << " generators:" << endl;
        for(unsigned int i=0; i<originalGenerators.size(); i++) {
            displayState(fout, originalGenerators[i]);
            fout << endl;
        }
        if (offsetGenerators.size() == 0)
            fout << "No parity-offset generators" << endl;
        else {
            fout << dec << offsetGenerators.size() << " parity-offset generators:" << endl;
            for(unsigned int i=0; i<offsetGenerators.size(); i++) {
                displayState(fout, offsetGenerators[i]);
                fout << endl;
            }
        }
        if (kernelGenerators.size() == 0)
            fout << "No parity-kernel generators" << endl;
        else {
            fout << dec << kernelGenerators.size() << " parity-kernel generators:" << endl;
            for(unsigned int i=0; i<kernelGenerators.size(); i++) {
                displayState(fout, kernelGenerators[i]);
                fout << endl;
            }
        }
    }
}

bool AffineSpaceOfStates::getOffsetWithGivenParity(PackedParity parity, vector<SliceValue>& output) const
{
    if (!packed) {
        vector<RowValue> parityUnpacked;
        for(unsigned int z=0; z<laneSize; z++)
            parityUnpacked.push_back(getParityFromPackedParity(parity, z));
        return getOffsetWithGivenParity(parityUnpacked, output);
    }
    else {
        output = offset;
        PackedParity outputParity = offsetParityPacked;
        PackedParity correctionParity = parity^offsetParityPacked;

        unsigned int i = 0;

        for(unsigned int xz=0; xz<(nrRowsAndColumns*laneSize); xz++) {
            PackedParity mask = ((PackedParity)1<<(xz+1))-1;
            PackedParity selectXZ = (PackedParity)1 << xz;
            if ((correctionParity & selectXZ) != 0) {
                while((i<offsetParitiesPacked.size()) && ((offsetParitiesPacked[i] & mask) != selectXZ))
                    i++;
                if (i<offsetParitiesPacked.size()) {
                    for(unsigned int z=0; z<laneSize; z++)
                        output[z] ^= offsetGenerators[i][z];
                    correctionParity ^= offsetParitiesPacked[i];
                }
                else
                    return false;
            }
        }
        return (correctionParity == 0);
    }
}

bool oneAndZeroesBefore(const vector<RowValue>& parity, RowValue maskX, RowValue selectX, unsigned int z)
{
    for(unsigned int iz=0; iz<z; iz++)
        if (parity[iz] != 0)
            return false;
    return (parity[z] & maskX) == selectX;
}

bool AffineSpaceOfStates::getOffsetWithGivenParity(const vector<RowValue>& parity, vector<SliceValue>& output) const
{
    if (packed)
        throw KeccakException("AffineBaseOfState initialized with PackedParity, not accessible without PackedParity.");
    output = offset;
    vector<RowValue> outputParity = offsetParity;
    vector<RowValue> correctionParity(parity);
    for(unsigned int z=0; z<laneSize; z++)
        correctionParity[z] ^= offsetParity[z];

    unsigned int i = 0;

    for(unsigned int z=0; z<laneSize; z++)
    for(unsigned int x=0; x<nrRowsAndColumns; x++) {
        RowValue maskX = (1<<(x+1))-1;
        RowValue selectX = 1 << x;
        if ((correctionParity[z] & selectX) != 0) {
            while((i<offsetParities.size()) && (!oneAndZeroesBefore(offsetParities[i], maskX, selectX, z)))
                i++;
            if (i<offsetParities.size()) {
                for(unsigned int jz=0; jz<laneSize; jz++)
                    output[jz] ^= offsetGenerators[i][jz];
                for(unsigned int jz=0; jz<laneSize; jz++)
                    correctionParity[jz] ^= offsetParities[i][jz];
            }
            else
                return false;
        }
    }
    for(unsigned int jz=0; jz<correctionParity.size(); jz++)
        if (correctionParity[jz] != 0)
            return false;
    return true;
}

SlicesAffineSpaceIterator AffineSpaceOfStates::getIteratorWithGivenParity(PackedParity parity) const
{
    vector<SliceValue> offset;
    
    if (getOffsetWithGivenParity(parity, offset))
        return SlicesAffineSpaceIterator(kernelGenerators, offset);
    else
        return SlicesAffineSpaceIterator();
}

SlicesAffineSpaceIterator AffineSpaceOfStates::getIteratorWithGivenParity(const vector<RowValue>& parity) const
{
    vector<SliceValue> offset;
    
    if (getOffsetWithGivenParity(parity, offset))
        return SlicesAffineSpaceIterator(kernelGenerators, offset);
    else
        return SlicesAffineSpaceIterator();
}

SlicesAffineSpaceIterator AffineSpaceOfStates::getIteratorInKernel() const
{
    vector<RowValue> parity(laneSize, 0);
    vector<SliceValue> offset;
    
    if (getOffsetWithGivenParity(parity, offset))
        return SlicesAffineSpaceIterator(kernelGenerators, offset);
    else
        return SlicesAffineSpaceIterator();
}

SlicesAffineSpaceIterator AffineSpaceOfStates::getIterator() const
{
    return SlicesAffineSpaceIterator(originalGenerators, offset);
}
