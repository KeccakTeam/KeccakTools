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

#ifndef _KECCAKFAFFINEBASES_H_
#define _KECCAKFAFFINEBASES_H_

#include <string>
#include "Keccak-fDCLC.h"
#include "Keccak-fParity.h"

using namespace std;

/** Like ListOfRowPatterns, this class lists the row patterns (differences or linear masks)
  * for a given row input pattern, but in the form of an affine space.
  * The members of the list are determined by offset plus any linear combination
  * of the generators.
  */
class AffineSpaceOfRows {
public:
    /** The generators. */
    vector<RowValue> generators;
    /** The offset. */
    RowValue offset;
public:
    /** The constructor creates an affine space with a single element: 0. */
    AffineSpaceOfRows() { offset = 0; }
    /** This method adds a generator to the set of generators.
      * @param   generator  The generator to add
      */
    void addGenerator(RowValue generator) { generators.push_back(generator); }
    /** This method sets the offset of the affine space.
      * @param   anOffset   The new offset value.
      */
    void setOffset(RowValue anOffset) { offset = anOffset; }
    /** This method displays the offset and generators.
      * @param   fout       The stream to display to.
      */
    void display(ostream& fout) const;
    /**
      * This method returns the number of generators.
      * After they are set, the number of generators happens to be equal
      * to the propagation weight of the input pattern
      * whose output space is represented by this affine space.
      * @return    The number of generators (propagation weight).
      */
    unsigned int getWeight() const { return generators.size(); }
};

/** This class expresses an affine space of slice values.
  * The members of the affine space are determined by the offset
  * plus any linear combination of the generators.
  *
  * The parity of a slice is the sum of all its rows.
  * In order to easily generate all the slice values with a given parity,
  * the generators are grouped in two sets.
  * - The parity-kernel generators (in kernelGenerators) are the generators with parity zero.
  * - The parity-offset generators (in offsetGenerators) are the generators with non-zero parity.
  *
  * The generators in offsetGenerators have linearly independent parities.
  * Hence, to generate all the slice values with a given requested parity, one has to
  * - first, take the parity of the slice offset
  * - then, combine it with some parity-offset generators to obtain the requested parity (if possible);
  *   when possible, this is uniquely determined thanks to the linear independence
  * - finally, combine it with the vector space of parity-kernel generators (which will not alter the parity)
  */
class AffineSpaceOfSlices {
public:
    /** The set of generators, before separation into parity-kernel and
      * parity-offset sets.
      */
    vector<SliceValue> originalGenerators;
    /** The set of parity-kernel generators of the affine space.
      * This generates a subspace with all-zero parity pattern.
      */
    vector<SliceValue> kernelGenerators;
    /** The set of parity-offset generators of the affine space.
      * This can be used to generate an offset for each of the parity patterns.
      * The parity-offset generators are organized such that their parity
      * makes an upper-triangular matrix. In other words, whenever
      * a generator has parity bit #n at 1 and parity bits #m&lt;n at 0,
      * the next generators necessarily have parity bits #m&lt;=n at 0.
      */
    vector<SliceValue> offsetGenerators;
    /** This vector contains the parities of offsetGenerators, i.e.,
      * offsetParities[i] contains the parity of offsetGenerators[i].
      */
    vector<RowValue> offsetParities;
    /** The offset of the affine space.
      */
    SliceValue offset;
    /** The parity of the offset of the affine space.
      */
    RowValue offsetParity;
public:
    /** Constructor of AffineSpaceOfSlices.
      * Warning: the passed vectors aGenerators and aGeneratorParities are modified.
      */
    AffineSpaceOfSlices(vector<SliceValue>& aGenerators, vector<RowValue>& aGeneratorParities, SliceValue aOffset, RowValue aOffsetParity);
    /** This method returns an offset slice (in argument output) with a given parity pattern.
      * If such a pattern could not be built, false is returned.
      */
    bool getOffsetWithGivenParity(RowValue parity, SliceValue& output) const;
    /** This method displays information about the object.
      */
    void display(ostream& fout) const;
private:
    void setGenerators(vector<SliceValue>& aGenerators, vector<RowValue>& aGeneratorParities);
};

/** This class implements an iterator over the affine space generated by the given
  * base and offset.
  */
template<class T>
class AffineSpaceIterator {
private:
    vector<vector<T> > *emptyBase;
    const vector<vector<T> > *base;
    vector<T> current;
    UINT64 i, end;
public:
    /** This constructor creates an empty iterator. */
    AffineSpaceIterator();
    /** This constructor initializes the affine space iterator with a given generator base
      * and a given offset.
      * @param   aBase      The generator base, as a reference to the set (vector) of generators.
      *                     Each generator is a vector of elements of template type @a T.
      * @param   aOffset    The offset, as a reference to a vector of elements of template type @a T.
      */
    AffineSpaceIterator(const vector<vector<T> >& aBase, const vector<T>& aOffset);
    /** The destructor. */
    ~AffineSpaceIterator();
    /** This method tells whether the last element of the affine space
      * has been reached.
      * @return True iff the last element of the affine space has been reached.
      */
    bool isEnd() const;
    /** This method moves the iterator to the next element in the affine space. */
    void operator++();
    /** This method returns a constant reference to the current element in the affine space.
      * @return The current element in the affine space.
      */
    const vector<T>& operator*() const;
    /** This method displays the offset and generators.
      * @param   fout       The stream to display to.
      */
    void display(ostream& fout) const;
    /** This method returns the number of elements in the affine space.
      * @return The number of elements in the affine space.
      */
    UINT64 getCount() const;
};

template<class T>
AffineSpaceIterator<T>::AffineSpaceIterator()
    : emptyBase(new vector<vector<T> >), base(emptyBase), current(), i(0), end(0)
{
}

template<class T>
AffineSpaceIterator<T>::AffineSpaceIterator(const vector<vector<T> >& aBase, const vector<T>& aOffset)
    : emptyBase(0), base(&aBase), current(aOffset), i(0), end((UINT64)1<<base->size())
{
}

template<class T>
AffineSpaceIterator<T>::~AffineSpaceIterator()
{
    if (emptyBase)
        delete emptyBase;
}

template<class T>
bool AffineSpaceIterator<T>::isEnd() const
{
    return i >= end;
}

template<class T>
void AffineSpaceIterator<T>::operator++()
{
    if (i < (end-1)) {
        unsigned int index = 0;
        while((i & ((UINT64)1<<index)) != 0)
            index++;
        for(unsigned int z=0; z<current.size(); z++)
            current[z] ^= (*base)[index][z];
    }
    i++;
}

template<class T>
const vector<T>& AffineSpaceIterator<T>::operator*() const
{
    return current;
}

template<class T>
void AffineSpaceIterator<T>::display(ostream& fout) const
{
    fout << "Offset: " << endl;
    for(unsigned int i=0; i<current.size(); i++)
        fout << current[i] << " ";
    fout << endl;
    fout << endl;
    fout << "Base: " << endl;
    for(unsigned int i=0; i<base->size(); i++) {
        fout << dec << i << ". ";
        for(unsigned int j=0; j<(*base)[i].size(); j++)
            fout << (*base)[i][j] << " ";
        fout << endl;
    }
}

template<class T>
UINT64 AffineSpaceIterator<T>::getCount() const
{
    return end;
}

/** This class implements an iterator over the affine space generated by the given
  * base and offset. Both are states expressed as vector of slices.
  */
typedef AffineSpaceIterator<SliceValue> SlicesAffineSpaceIterator;


/** This class expresses an affine space of states.
  * The members of the affine space are determined by the offset
  * plus any linear combination of the generators.
  * The generators are split into two sets:
  * - one set of generators that have a zero parity (called "parity-kernel"), and
  * - one set of generators that have a non-zero parity (called "parity-offset").
  * This allows one to address a subspace of states with a given parity.
  */
class AffineSpaceOfStates {
public:
    /** The set of generators, before separation into parity-kernel and
      * parity-offset sets.
      */
    vector<vector<SliceValue> > originalGenerators;
    /** The set of parities, before separation into parity-kernel and
      * parity-offset sets, i.e., originalParities[i] contains the parity
      * of originalGenerators[i].
      */
    vector<vector<RowValue> > originalParities;
    /** The set of parity-kernel generators of the affine space.
      */
    vector<vector<SliceValue> > kernelGenerators;
    /** The set of parity-offset generators of the affine space.
      * This can be used to generate an offset for a given parity.
      * The parity-offset generators are organized such that their parity
      * makes an upper-triangular matrix. In other words, whenever
      * a generator has parity bit <i>n</i> set to 1
      * and parity bits <i>m</i>&lt;<i>n</i> set to 0,
      * the next generators necessarily have parity bits <i>m</i>&lt;=<i>n</i> set to 0.
      */
    vector<vector<SliceValue> > offsetGenerators;
    /** If packed is true, this vector contains the parities of offsetGenerators, i.e.,
      * offsetParitiesPacked[i] contains the parity of offsetGenerators[i].
      */
    vector<PackedParity> offsetParitiesPacked;
    /** If packed is false, this vector contains the parities of offsetGenerators, i.e.,
      * offsetParities[i] contains the parity of offsetGenerators[i].
      */
    vector<vector<RowValue> > offsetParities;
    /** The offset of the affine space.
      */
    vector<SliceValue> offset;
    /** If packed is true, the parity of the offset of the affine space.
      */
    PackedParity offsetParityPacked;
    /** If packed is false, the parity of the offset of the affine space.
      */
    vector<RowValue> offsetParity;
    /** This attribute indicates whether we are using the packed parities or not.
      */
    bool packed;
protected:
    /** The lane size. */
    unsigned int laneSize;
public:
    /** This constructor initializes the different attributes from the given generators,
      * the offset and their parities.
      * This function is to be called only if the number of slices is low enough
      * so that PackedParity can contain all the parities.
      * @note The passed vectors aGenerators and aGeneratorParities are modified in the process.
      * @param   aLaneSize          The lane size.
      * @param   aGenerators        The set of generators of the affine space,
      *                             each given as a vector of slices.
      * @param   aGeneratorParities The corresponding parities, each given in a PackedParity type.
      * @param   aOffset            The offset of the affine space, as a vector of slices.
      * @param   aOffsetParity      The parity of the offset, given in a PackedParity type.
      */
    AffineSpaceOfStates(unsigned int aLaneSize, vector<vector<SliceValue> >& aGenerators, vector<PackedParity>& aGeneratorParities, const vector<SliceValue>& aOffset, PackedParity aOffsetParity);
    /** This constructor initializes the different attributes from the given generators,
      * the offset and their parities.
      * @note   The passed vectors aGenerators and aGeneratorParities are modified in the process.
      * @param   aLaneSize          The lane size.
      * @param   aGenerators        The set of generators of the affine space,
      *                             each given as a vector of slices.
      * @param   aGeneratorParities The corresponding parities, each given as a vector of rows.
      * @param   aOffset            The offset of the affine space, as a vector of slices.
      * @param   aOffsetParity      The parity of the offset, given as a vector of rows.
      */
    AffineSpaceOfStates(unsigned int aLaneSize, vector<vector<SliceValue> >& aGenerators, vector<vector<RowValue> >& aGeneratorParities, const vector<SliceValue>& aOffset, const vector<RowValue>& aOffsetParity);
    /** This method returns a state value (in argument @a output) with a given parity.
      * From the offset and the parity-offset generators of the affine space, this method
      * computes an element that has the given parity. (Note that other elements with the same parity
      * can be generated by adding any linear combination of parity-kernel generators.)
      * If the given parity cannot be reached, false is returned.
      * @param   parity     The requested parity.
      * @param   output     The state returned by the method with the requested parity, if possible.
      * @return The method returns true iff a state can be found in the affine space with the requested parity.
      */
    bool getOffsetWithGivenParity(PackedParity parity, vector<SliceValue>& output) const;
    /** This method returns a state value (in argument @a output) with a given parity.
      * From the offset and the parity-offset generators of the affine space, this method
      * computes an element that has the given parity. (Note that other elements with the same parity
      * can be generated by adding any linear combination of parity-kernel generators.)
      * If the given parity cannot be reached, false is returned.
      * @param   parity     The requested parity.
      * @param   output     The state returned by the method with the requested parity, if possible.
      * @return The method returns true iff a state can be found in the affine space with the requested parity.
      */
    bool getOffsetWithGivenParity(const vector<RowValue>& parity, vector<SliceValue>& output) const;
    /** This method returns an iterator to all states in the affine space.
      * @return The iterator to the whole affine space.
      */
    SlicesAffineSpaceIterator getIterator() const;
    /** This method returns an iterator to the affine space restricted to a given parity.
      * In the returned iterator, its offset is computed as in method getOffsetWithGivenParity()
      * and its generators are the parity-kernel generators.
      * @param   parity     The requested parity.
      * @return An iterator to all states in the affine space, restricted to the requested parity.
      */
    SlicesAffineSpaceIterator getIteratorWithGivenParity(PackedParity parity) const;
    /** This method returns an iterator to the affine space restricted to a given parity.
      * In the returned iterator, its offset is computed as in method getOffsetWithGivenParity()
      * and its generators are the parity-kernel generators.
      * @param   parity     The requested parity.
      * @return An iterator to all states in the affine space, restricted to the requested parity.
      */
    SlicesAffineSpaceIterator getIteratorWithGivenParity(const vector<RowValue>& parity) const;
    /** Like getIteratorWithGivenParity() with the requested parity being zero.
      * @return An iterator to all states in the affine space with zero parity (i.e., in the kernel).
      */
    SlicesAffineSpaceIterator getIteratorInKernel() const;
    /** This method displays the offset and generators.
      * @param   fout       The stream to display to.
      */
    void display(ostream& fout) const;
private:
    void setGenerators(vector<vector<SliceValue> >& aGenerators, vector<vector<RowValue> >& aGeneratorParities);
    void setGenerators(vector<vector<SliceValue> >& aGenerators, vector<PackedParity>& aGeneratorParities);
};


#endif
