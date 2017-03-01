/*
KeccakTools

This code implements the techniques described in FSE2017 paper 
"New techniques for trail bounds and application to differential trails in Keccak"
by Silvia Mella, Joan Daemen and Gilles Van Assche.
http://tosc.iacr.org/
http://eprint.iacr.org/2017/181

Implementation by Silvia Mella, hereby denoted as "the implementer".

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/
*/

#ifndef _KECCAKFTRAILEXTENSIONBASEDONPARITY_H_
#define _KECCAKFTRAILEXTENSIONBASEDONPARITY_H_

#include <algorithm>
#include <math.h>
#include <string>
#include "Keccak-f.h"
#include "Keccak-fAffineBases.h"
#include "Keccak-fDCLC.h"
#include "Keccak-fParts.h"
#include "Keccak-fPropagation.h"
#include "Keccak-fTrailExtension.h"

using namespace std;

/** This method builds the χ<sup>-1</sup>-envelope space of a given input state.
* This space cotains the set of states compatible with the given input state through χ<sup>-1</sup>.
* @param   stateAfterChi  The state after χ to propagate, given as a vector of slices.
* @return The affine space as a AffineSpaceOfStates object.
*/
AffineSpaceOfStates buildBasisBeforeChiGivenPatternAfterChi(const vector<SliceValue>& stateAfterChi);

/** This method returns the image of a given input affine space through ρ<sup>-1</sup> and π<sup>-1</sup>.
* @param   keccakFTE  The keccak propagation context.
* @param   basisBeforeChi  The affine space after ρ and π, given as a AffineSpaceOfStates object.
* @return The affine space before ρ and π as a AffineSpaceOfStates object.
*/
AffineSpaceOfStates getdBasisAfterThetaGivenPatternBeforeChi(const KeccakFTrailExtension& keccakFTE, const AffineSpaceOfStates& basisBeforeChi);

/** This method returns a basis for the intersection of @basis and the kernel.
* @param   basis  The affine space to intersect with the kernel.
* @return  The intersection space as a AffineSpaceOfStates object.
*/
AffineSpaceOfStates buildBasisIntersectionWithKernel(const AffineSpaceOfStates& basis);

/** This method checks if @basis intersects the kernel.
* @param   basis  The affine space to intersect with the kernel.
* @return True if the intersection is non-empty. False otherwise.
*/
bool intersectionWithKernel(const AffineSpaceOfStates& basis);

/** This method sets a basis for the output row difference patterns through χ for each input row difference pattern.
* @param   keccakFTE  The keccak propagation context.
* @param   aBasisPerInput The vector of bases. The index of the vector corresponds to the input row difference pattern.
*/
void setBasisPerInput(const KeccakFTrailExtension& keccakFTE, vector<AffineSpaceOfRows>& aBasisPerInput);

/** This method builds the affine space of states compatible with a given input state through χ.
* @param affinePerInput The output row difference patterns through χ.
* @param   stateBeforeChi  The state before χ to propagate, given as a vector of slices.
* @return The affine space as a AffineSpaceOfStates object.
*/
AffineSpaceOfStates buildBasisAfterChiGivenPatternBeforeChi(vector<AffineSpaceOfRows>& affinePerInput, const vector<SliceValue>& stateBeforeChi);

/** This method computes the parity pattern compatible with a given parity pattern through theta.
* @param A The input parity pattern as a vector of LaneValue.
* @param laneSize The lane size of the context.
*/
void inverseThetaOnParity(vector<LaneValue>& A, unsigned int& aLaneSize);

/** This method computes the parity pattern compatible with a given parity pattern through theta.
* @param A The input parity pattern as a vector of RowValue.
* @param laneSize The lane size of the context.
*/
void inverseThetaOnParity(vector<RowValue>& A, unsigned int& aLaneSize);

/** This method counts, for each column position, the number of basis vectors with active bit in that column.
* @param basis The basis vectors.
* @return The number of basis vectors with active bits for each column position.
*/
vector<vector<unsigned int> > getnNrBasisVectorsPerColumn(AffineSpaceOfStates& basis);

/** This method counts, for each column position, the number of basis vectors with active bit in that column excluding orbitals.
* @param basis The basis vectors.
* @return The number of basis vectors with active bits for each column position.
*/
vector<vector<unsigned int> > getnNrBasisVectorsPerColumnNoOrbitals(AffineSpaceOfStates& basis);

/** This method computes the space of all the possible values of a row in the parity plane, given a basis for that row.
* @param offset The offset for the row.
* @param basis The basis for the row.
* @return The set of possible row values.
*/
vector<RowValue> getSetOfRowValues(RowValue& offset, vector<RowValue>& basis);

/** This method computes the space of all the possible values for each row in the parity plane, given a basis for the parity plane.
* @param basis The basis.
* @return The space of possible row values.
*/
vector< vector<RowValue> > getRowValuesFromBasis(AffineSpaceOfStates& basis);

/** This method looks for the slice to start from for backward extension outside the kernel.
* That is the slice with maximum number of empty slices before it.
* @param basisAfterTheta The basis of the affine space after theta.
* @return The z-coordinate of the starting slice.
*/
unsigned int getStartingSlice(AffineSpaceOfStates& basisAfterTheta);


/** This base class represents an iterator on a set of parity patterns.
*/
class parityIterator {
protected:
	/* The propagation context as a reference to a KeccakFPropagation object. */
	const KeccakFPropagation& DCorLC;
	/* The lane size. */
	unsigned int laneSize;
public:
	/** The default constructor.
	* @param DCorLC The propagation context of the parity, as a reference to a KeccakFPropagation object.
	* @param laneSize The lane size pf the context
	*/
	parityIterator(const KeccakFPropagation& aDCorLC)
		: DCorLC(aDCorLC), laneSize(aDCorLC.laneSize){}

	/** This method indicates whether the iterator has reached the end of the set of parities.
	*  @return True if there are no more parities to consider.
	*/
	virtual bool isEnd() = 0;
	/** This method indicates whether the set of parities is empty.
	*  @return True if the set is empty.
	*/
	virtual bool isEmpty() = 0;
	/** This method indicates whether the iterator knows the cardinality of the set.
	*  @return True if the cardinality is known.
	*/
	virtual bool isBounded() = 0;
	/** This method returns the index of the current parity pattern encountered,
	* starting from 0 for the first parity pattern.
	*  @return The index of the current parity pattern.
	*/
	virtual UINT64 getIndex() = 0;
	/** If the cardinality of the set is known, this method returns it.
	*  @return The cardinality of the set.
	*/
	virtual UINT64 getCount() = 0;
	/** The '++' operator moves the iterator to the next parity pattern in the set.
	*/
	virtual void operator++() = 0;
	/** The '*' operator gives a constant reference to the current parity pattern.
	*  @return A constant reference to the current parity pattern.
	*/
	virtual const vector<RowValue>& operator*() = 0;
};


/** This class represents an iterator on a set of parity patterns compatible with a given state pattern through χ<sup>-1</sup>.
*/
class parityBackwardIterator : public parityIterator
{
public:
	/** The slice which the backward extension starts from. */
	unsigned int start;
	/** The parity of the offset of the basis after theta. */
	vector<RowValue> offsetParity;
	/** The offset of the basis vectors after theta. */
	vector<SliceValue> offset;
	/** The space of all possible values for each row of the parity plane. */
	vector< vector<RowValue> > rowsValues;
	/** This vector contains a reference to the current value of each row of the parity plane.  */
	vector<int> indexes;
	/** The number of basis vectors with active bit for each column of the parity plane. */
	vector<vector<unsigned int> > nrBasisVectors;
	/** The parity pattern before theta. */
	vector<RowValue> a;
	/** The parity pattern after theta. */
	vector<RowValue> b;
	/** The working z-ccordinate. */
	unsigned int current;
	/** The lower bound on the current parity pattern. */
	vector<unsigned int> lowerBound;
	/** The current parity pattern. */
	vector<RowValue> parity;
	/** The maximum weight to generate states. */
	int maxWeight;
    /** Attribute that indicates whether the iterator has reached the end. */
	bool end;
	/** Attribute that indicates whether the iterator has been initialized. */
	bool initialized;
	/** Attribute that indicates whether the set of parities is empty. */
	bool empty;
	/** Number of the current iteration. */
	UINT64 index;
public:
	/** The default constructor.
	* @param DCorLC The propagation context of the parity, as a reference to a KeccakFPropagation object.
	* @param aOffset The offset of the basis after theta.
	* @param aOffsetParity The parity of the offset of the basis after theta.
	* @param aValuesX The set of possible values of the rows in the parity plane.
	* @param aNrBasisVectors  The number of basis vectors with active bit for each column of the parity plane.
	* @param aStart The z-coordinate of the starting slice.
	* @param aGuess The guess on the value of the starting slice.
	* @param aMaxWeight The maximum total weight to consider.
	*/
	parityBackwardIterator(const KeccakFPropagation& aDCorLC, vector<SliceValue> &aOffset, vector<RowValue> &aOffsetParity,
		vector<vector<RowValue> > &aValuesX, vector<vector<unsigned int> > &aNrBasisVectors, unsigned int& aStart,
		RowValue& aGuess, unsigned int aMaxWeight);

	/** This method initializes the iterator, by calling the method first. */
	void initialize();

	/** This method method generates the first valid parity pattern.
	* @return true if a first valid pattern exists, false otherwise.
	*/
	bool first();

	/** This method generates the next valid parity pattern.
	* @return true if a next valid pattern exists, false otherwise.
	*/
	bool next();

	/** This method iterates a given row value and pushes it to the parity plane.
	* @param workingRow The current row value to iterate.
	* @return true if a value exists, false otherwise.
	*/
	bool iterateRowAndPush(RowValue& workingRow);

	/** This method computes the next possible value of a given row.
	* @param workingRow The current row value.
	* @return true if a value exists, false otherwise.
	*/
	bool successorOf(RowValue& workingRow);

	/** This method pushes a given row to the parity plane.
	* @param workingRow The row value to push.
	*/
	void pushRow(RowValue& workingRow);

	/** This method pops the highest row value from the parity pattern and returns it.
	* @param workingRow The highest row value.
	* @return true if a row has been popped, false otherwise.
	*/
	bool popRow(RowValue& workingRow);

	/** This method computes the minimum number of active rows before theta in a slice.
	* @param z The z-coordinate of the slice.
	* @return The minimum number of active rows.
	*/
	unsigned int getNrActiveRowsBeforeTheta(unsigned int& z);

	/** This method computes the number of possible orbitals in a column.
	* @param x The x-coordinate of the column.
	* @param z The z-coordinate of the column.
	* @param activeBasis The basis vectors.
	* @return The number of orbitals.
	*/
	unsigned int getNrOrbitals(unsigned int& x, unsigned int& z, unsigned int& activeBasis);

	/** This method checks if lowerBound is smaller than maximum weight.
	* @return true if the lower bound is maller than maxWeight, false otherwise.
	*/
	bool checkWeight();

	/** This method checks if the current parity pattern is compatible with the guess on the starting slice.
	* @return true if the pattern is compatible, false otherwise.
	*/
	bool checkCompatibility();

	/** This method sets the current parity pattern. */
	void getParity();
	/** See parityIterator::isEnd(). */
	bool isEnd();
	/** See parityIterator::isEmpty(). */
	bool isEmpty();
	/** See parityIterator::isBounded().*/
	bool isBounded();
	/** See parityIterator::getIndex(). */
	UINT64 getIndex();
	/** See parityIterator::getCount(). */
	UINT64 getCount();
	/** See parityIterator::operator++(). */
	void operator++();
	/** See parityIterator::operator*(). */
	const vector<RowValue>& operator*();

};



/** This class represents an iterator on the possible values of a state after chi given the state before chi.*/
class stateForwardIterator
{
protected:
	/**  The propagation context of the parity as a reference to a KeccakFPropagation object.*/
	const KeccakFPropagation& DCorLC;
	/** The lane size. */
	unsigned int laneSize;
public:
	/** The basis of the affine space of states after chi. */
	vector<vector<SliceValue> > basis;
	/** The z-coordinate of the active slice for each basis vector. */
	vector<unsigned int> slices;
	/** The number of basis vectors with active bits for each column position. */
	vector<vector<unsigned int> > numPerColumn;
	/** The parity of each basis vector for the affine space of states after chi. */
	vector<vector<RowValue> > basisParity;
	/** The offset of the affine space of states after chi. */
	vector<SliceValue> offset;
	/** The parity of the offset of the affine space of states after chi. */
	vector<RowValue> offsetParity;
	/** The indexes of the basis vectors already added to the current state. */
	vector<int> indexes;
	/** The current state after chi. */
	vector<vector<SliceValue> > stateAtA;
	/** The current state after lambda. */
	vector<vector<SliceValue> > stateAtB;
	/** The current weight of the state after lambda. */
	vector<unsigned int> weight;
	/** The maximum weight to generate states. */
	unsigned int budgetWeight;
	/** The index of the first basis vector in the Kernel. */
	int firstOrbital;
	/** The current parity pattern. */
	vector<vector<RowValue> > C;
	/** The current theta effect. */
	vector<vector<RowValue> > D;
	/** Attribute that indicates whether the iterator has reached the end. */
	bool end;
	/** Attribute that indicates whether the iterator has been initialized. */
	bool initialized;
	/** Attribute that indicates whether the set of states is empty. */
	bool empty;
	/** Number of the current iteration. */
	UINT64 index;
public:
	/** The default constructor.
	* @param DCorLC The propagation context of the parity, as a reference to a KeccakFPropagation object.
	* @param aBasis The basis of the affine space.
	* @param aBudgetWeight The maximum weight to generate states.
	*/
	stateForwardIterator(const KeccakFPropagation& aDCorLC, AffineSpaceOfStates& aBasis, unsigned int aBudgetWeight);

	/** This method initializes the iterator, by calling the method first. */
	void initialize();

	/** This method generates the first valid pattern.
	* @return true if a valid pattern exists, false otherwise.
	*/
	bool first();

	/** This method generats the next valid pattern.
	* @return true if a valid pattern exists, false otherwise.
	*/
	bool next();

	/** This method iterates a given vector index and pushes the corresponding vector to the state.
	* @param workingIndex The index value to iterate.
	* @return true if a valid value exists, false otherwise.
	*/
	bool iterateIndexAndPush(int& workingIndex);

	/** This method computes the next vector index.
	* @param workingIndex The current vector index.
	* @return true if a value exists, false otherwise.
	*/
	bool successorOf(int& workingIndex);
	/** This method check whether the addition of a basis vector is within the budget.
	* @param workingIndex The index of the basis vector to add.
	* @return true if the addition is affordable, false otherwise.
	*/
	bool canAfford(int& workingIndex);
	/** This method pushes a given basis vector to the parity plane.
	* @param workingIndex The index of the basis vector to push.
	*/
	void pushIndex(int& workingIndex);
	/** This method pops the highest basis vector from the state pattern and returns it.
	* @param workingIndex The index of the highest basis vector.
	* @return true if a vector has been popped, false otherwise.
	*/
	bool popIndex(int& workingIndex);
	/** This method checks whether among the remaining basis vectors, there are any that might decrease the cost.
	* These are vectors that are bit-overlapping with other vectors already in the state, or that modify the parity of the state, or that add bits to affected columns of the state.
	* @return true if one of the remaining basis vectors might decrease the cost, false otherwise.
	*/
	bool thereAreSuperEntangledIndexes();
	/** This method checks whether a given basis vector of index @a overlaps any basis vector already in the state.
	* @param a The index of the given basis vector.
	* @return true if the given basis vector overlaps any other basis vector in the state, false otherwise.
	*/
	bool isBitOverlapping(unsigned int& a);
	/** This method checks whether a given basis vector of index @a is adding bits to any affected column of the state.
	* @param a The index of the given basis vector.
	* @return true if the given basis vector adds bits to any affected column, false otherwise.
	*/
	bool isAddingBitToAC(unsigned int& a);

	/** This method checks whether a given basis vector of index @a is turning an unaffected odd column to affected or an affected even column to affected odd.
	* @param a The index of the given basis vector.
	* @return true if the given basis vector turns an unaffected odd column to affected or an affected even column to affected odd, false otherwise.
	*/
	bool isEntangled(unsigned int& a);
	/** This method checks whether a given basis vector of index @a is changing a run of the parity of the current state.
	* @param a The index of the given basis vector.
	* @return true if the given basis vector changes the parity, false otherwise.
	*/
	bool isRunModifying(unsigned int& a);

	/** This method indicates whether the iterator has reached the end of the set.
	*  @return True if there are no more states to consider.
	*/
	bool isEnd();
	/** This method indicates whether the set of states is empty.
	*  @return True if the set is empty.
	*/
	bool isEmpty();
	/** This method indicates whether the iterator knows the cardinality of the set.
	*  @return True if the cardinality is known.
	*/
	bool isBounded();
	/** This method returns the index of the current state pattern encountered,
	* starting from 0 for the first state pattern.
	*  @return The index of the current state pattern.
	*/
	UINT64 getIndex();
	/** If the cardinality of the set is known, this method returns it.
	*  @return The cardinality of the set.
	*/
	UINT64 getCount();
	/** The '++' operator moves the iterator to the next state in the set.
	*/
	void operator++();
	/** The '*' operator gives a constant reference to the current state pattern.
	*  @return A constant reference to the current state pattern.
	*/
	const vector<SliceValue>& operator*();
};





class KeccakFTrailExtensionBasedOnParity : public KeccakFTrailExtension
{
public:

	/** The constructor.
	* See KeccakFTrailExtension::KeccakFTrailExtension(const KeccakFDCLC& aParent, KeccakFPropagation::DCorLC aDCorLC).
	*/
	KeccakFTrailExtensionBasedOnParity(const KeccakFDCLC& aParent, KeccakFPropagation::DCorLC aDCorLC) :
		KeccakFTrailExtension(aParent, aDCorLC) {};

	/** Starting from a given trail (prefix or core), this method
	* appends states to it and systematically looks
	* for all trails with @nrRounds rounds
	* up to total weight @a maxTotalWeight
	* that have the given trail as prefix
	* and whose states after chi are in the kernel.
	* @param  trail   The starting trail core or trail prefix.
	* @param  trailsOut   Where to output the found trails.
	* @param  nrRounds    The target number of rounds.
	* @param  maxTotalWeight  The maximum total weight to consider.
	*/
	void forwardExtendTrailInTheKernel(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight);
	
	/** This function is like forwardExtendTrailInTheKernel(), except that it processes
	* all the trails from @a trailsIn.
	* @param  trailsIn    The starting trail cores or trail prefixes.
	* @param  trailsOut   Where to output the found trails.
	* @param  nrRounds    The target number of rounds.
	* @param  maxTotalWeight  The maximum total weight to consider.
	*/
	void forwardExtendTrailsInTheKernel(TrailIterator& trailsIn, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight);

	/** Starting from a given trail (prefix or core), this method
	* prepends states to it, and systematically looks
	* for all trails with @nrRounds rounds
	* up to total weight @a maxTotalWeight
	* that have the given trail as suffix.
	* and whose states after chi are in the kernel.
	* @param  trail   The starting trail core or trail prefix.
	* @param  trailsOut   Where to output the found trails.
	* @param  nrRounds    The target number of rounds.
	* @param  maxTotalWeight  The maximum total weight to consider.
	*/
	void backwardExtendTrailInTheKernel(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight);
	
	/** This function is like backwardExtendTrails(), except that it processes
	* all the trails from @a trailsIn.
	* @param  trailsIn    The starting trail cores or trail prefixes.
	* @param  trailsOut   Where to output the found trails.
	* @param  nrRounds    The target number of rounds.
	* @param  maxTotalWeight  The maximum total weight to consider.
	*/
	void backwardExtendTrailsInTheKernel(TrailIterator& trailsIn, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight);

	/** Starting from a given trail (prefix or core), this method
	* prepends states to it, and systematically looks
	* for all trails with @a nrRounds rounds
	* up to total weight @a maxTotalWeight
	* that have the given trail as suffix
	* and whose states after chi are outside the kernel.
	* If the given trail is a trail core, then the function looks for all trail cores,
	* unless @a allPrefixes is set to true, in which case all trail prefixes
	* are output.
	* Otherwise if the given trail is a trail prefix, then the function looks for all trail prefixes.
	* @param  trail   The starting trail core or trail prefix.
	* @param  trailsOut   Where to output the found trails.
	* @param  nrRounds    The target number of rounds.
	* @param  maxTotalWeight  The maximum total weight to consider.
	*/
	void backwardExtendTrailOutsideKernel(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight);
	
	/** This function is like backwardExtendTrailOutsideKernel(), except that it processes
	* all the trails from @a trailsIn.
	* @param  trailsIn    The starting trail cores or trail prefixes.
	* @param  trailsOut   Where to output the found trails.
	* @param  nrRounds    The target number of rounds.
	* @param  maxTotalWeight  The maximum total weight to consider.
	*/
	void backwardExtendTrailsOutsideKernel(TrailIterator& trailsIn, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight);

	/** Starting from a given trail (prefix or core), this method
	* appends states to it and systematically looks
	* for all trails with @nrRounds rounds
	* up to total weight @a maxTotalWeight
	* that have the given trail as prefix
	* and whose states after chi are outside the kernel.
	* @param  trail   The starting trail core or trail prefix.
	* @param  trailsOut   Where to output the found trails.
	* @param  nrRounds    The target number of rounds.
	* @param  maxTotalWeight  The maximum total weight to consider.
	*/
	void forwardExtendTrailsOutsideKernel(TrailIterator& trailsIn, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight);
	
	/** This function is like forwardExtendTrailsOutsideKernel(), except that it processes
	* all the trails from @a trailsIn.
	* @param  trailsIn    The starting trail cores or trail prefixes.
	* @param  trailsOut   Where to output the found trails.
	* @param  nrRounds    The target number of rounds.
	* @param  maxTotalWeight  The maximum total weight to consider.
	*/
	void forwardExtendTrailOutsideKernel(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight, vector<AffineSpaceOfRows> basisPerInput);

protected:

	void recurseForwardExtendTrailInTheKernel(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight);
	void recurseBackwardExtendTrailInTheKernel(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight);
	void recurseForwardExtendTrailOutsideTheKernel(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight, vector<AffineSpaceOfRows> basisPerInput);
	void recurseBackwardExtendTrailOutsideKernel(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight, int maxRevWeight, bool allPrefixes);

};

#endif
