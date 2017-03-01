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

/*! \file Keccak-fTree.h

This file defines Orbital-Tree and Run-Tree, used to generate 2-round trail cores in Keccak-f,
as described in the paper.
*/

#ifndef _KECCAKFTREE_
#define _KECCAKFTREE_

#include <stack>
#include <vector>
#include <string>
#include "Keccak-fTrailCoreParity.h"
#include "Keccak-fTree.h"
#include "Tree.h"


class TwoRoundTrailCoreStack;
/**
* \brief OrbitalsSet class : represents the set of orbitals.
* \details This class represents the set of orbitals and defines the order relation among them.
* The orbitals are possibly constrained by the position of other bits already
* present in the state, as captured by the yMin attribute.
*/
class OrbitalsSet {
public:
	/** This attribute indicates whether the orbitals are used to generate states in the kernel or not. */
	bool kernel;
	/** The minimum y position of the lower bit of an orbital for each column of the state. */
	vector<unsigned int> yMin;
	/** The lane size of Keccak-f. */
	unsigned int laneSize;

public:
	/**
	* The default constructor.
	*/
	OrbitalsSet() {};

	/**
	* The constructor, for orbitals in the kernel.
	* @param aLaneSize the lane size.
	*/
	OrbitalsSet(const unsigned int aLaneSize)
		: kernel(true), yMin(5 * aLaneSize, 0), laneSize(aLaneSize) {}

	/**
	* The constructor, with yMin specified.
	* @param ayMin the minimum position of the lower bit of an orbital for each column of the state.
	* @param aLaneSize the lane size.
	*/
	OrbitalsSet(const vector<unsigned int>& ayMin, const unsigned int aLaneSize)
		: kernel(false), yMin(ayMin), laneSize(aLaneSize) {}

	/** This method returns an orbital in the first available position
	* with respect to the order relation [z,x,y0,y1] and restrictions given by yMin.
	* @param unitList the list of units.
	* @return the first available orbital position.
	*/
	OrbitalPosition getFirstChildUnit(const std::vector<OrbitalPosition>& unitList, const TwoRoundTrailCoreStack& cache) const;

	/** This method iterates the current orbital with respect to the order relation [z,x,y0,y1] and restrictions given by yMin.
	* @param unitList the list of units.
	* @param current the current orbital.
	*/
	void iterateUnit(const std::vector<OrbitalPosition>& unitList, OrbitalPosition& current, const TwoRoundTrailCoreStack& cache) const;

	/** This method compares two given orbitals with respect to the order relation [z,x,y0,y1].
	* @param first the first given orbital.
	* @param second the second given orbital.
	* @return 0 if the two orbitals are equal, 1 if the first is smaller, 2 if the second is smaller.
	*/
	unsigned int compare(const OrbitalPosition& first, const OrbitalPosition& second) const;

	/** This method checks if a list of orbitals is z-canonical with respect to the order relation [z,x,y0,y1].
	* @param orbitalList the list of orbitals.
	* @param cache a reference to the cache representing the trail
	* @return true if the given list is z-canonical, false otherwise.
	*/
	bool isCanonical(const std::vector<OrbitalPosition>& orbitalList, TwoRoundTrailCoreStack& cache) const;
};


/**
* \brief Column class : represents column assignments.
* \details This class represents a column assignment for Keccak-f.
*/
class Column {

public:
	/** The (x,z) coordinates of the column. */
	ColumnPosition position;
	/** The value of the column. */
	ColumnValue value;
	/** The parity of the column. */
	bool odd;
	/** The theta-effect on the column. */
	bool affected;
	/** This attribute indicates whether the column is entangled with another one. Namely, if is an odd-0 overlapping an Affected even or vice versa or not.*/
	bool entangled;
	/** Whether the column starts a run or not. */
	bool starting;
	/** Parameter used to optimize the iteration of the column value. */
	unsigned int index;

public:

	/** The default construtor. */
	Column();

	/** The constructor.
	* @param odd the parity of the column.
	* @param affected the theta-effect on the column.
	*/
	Column(bool& odd, bool& affected);

};

/**
* \brief ColumnsSet class : represents the set of column assignments.
* \details This class represents the set of column assignments and defines the order relation among them.
*/
class ColumnsSet {

public:
	/** The lane size. */
	unsigned int laneSize;
	/** The possible values of Affected Even columns, Affected Odd columns and Unaffected Odd columns. */
	static const ColumnValue AEValues[16], AOValues[16], UOValues[5];

public:

	/**
	* The default constructor.
	*/
	ColumnsSet() { };

	/**
	* The constructor.
	* @param aLaneSize the lane size
	*/
	ColumnsSet(const unsigned int aLaneSize)
		: laneSize(aLaneSize) { }

	/** This method returns a column assignment in the first available position
	* with respect to the order relation among runs.
	* @param unitList the list of units.
	* @return the first available column assignment.
	*/
	Column getFirstChildUnit(const std::vector<Column>& unitList, const TwoRoundTrailCoreStack& cache) const;

	/** This method iterates the current unit with respect to the order relation [z,x,value] and restrictions on its type.
	* @param unitList the list of units.
	* @param current the current column assignment.
	*/
	void iterateUnit(const std::vector<Column>& unitList, Column& current, const TwoRoundTrailCoreStack& cache) const;

	/** This method compares two given column assignments with respect to the order relation [z,x,value].
	* @param first the first given column assignment.
	* @param second the second given column assignment.
	* @return 0 if the two column assignments are equal, 1 if the first is smaller, 2 if the second is smaller.
	*/
	unsigned int compare(const Column& first, const Column& second) const;

	/** This method checks if a list of column assignments is z-canonical with respect to the order relation [z,x,value].
	* @param unitList the list of column assignments.
	* @param cache a reference to the cache representing the trail
	* @return true if the given list is z-canonical, false otherwise.
	*/
	bool isCanonical(const std::vector<Column>& unitList, TwoRoundTrailCoreStack& cache) const;

private:
	/** This method checks if a given column assignment overlaps any other column of a given unit list.
	* The given column assignment's entanglement attribute is initialized.
	* @param unitList the list of column assignments.
	* @param current the given column assignment.
	* @return true if the given column overlaps a column in @unitList.
	*/
	bool checkColumnOverlapping(const std::vector<Column>& unitList, Column& current, const TwoRoundTrailCoreStack& cache) const;

};


/**
* \brief TwoRoundTrailCoreStack class : cache representation for 2-round trail cores in Keccak-f.
* \details This class represents a 2-round trail core as a node of a tree.
*/
class TwoRoundTrailCoreStack {
public:
    /** The propagation context of the trail, as a reference to a KeccakFPropagation object. */
	const KeccakFPropagation& DCorLC;
	/** The lane size. */
	unsigned int laneSize;
	/** State at A and state at B. */
	stack<StateAsVectorOfSlices> stack_stateAtA, stack_stateAtB;
	/** The stack for the minimum reverse weight of A and the weight of B. */
	stack<unsigned int> stack_w0, stack_w1;
	/** The stack to indicate if a state is valid or not. */
	stack<bool> stack_complete;
	/** The z-period of the root and the current node. */
	unsigned int rootPeriod, nodePeriod;
	/** The parity pattern of the current state. */
	vector<RowValue> C;
	/** The theta effect of the current state. */
	vector<RowValue> D;
	/** This variable indicates whether the last push was a dummy push or not. */
	bool dummy;

public:

	/**
	* The default constructor.
	*/
	TwoRoundTrailCoreStack();

	/**
	* The constructor.
	* @param aDCorLC The propagation context of the trail, as a reference to a KeccakFPropagation object.
	*/
	TwoRoundTrailCoreStack(const KeccakFPropagation& aDCorLC);

	/**
	* The constructor.
	* @param aDCorLC The propagation context of the trail, as a reference to a KeccakFPropagation object.
	* @param stateA The initial state at A.
	* @param stateB The initial state at B.
	* @param aW0 The initial minimum reverse weight of A.
	* @param aW1 The initial weight of B.
	* @param aComplete Completness of the initial state.
	* @param aRootPeriod The z-period of the root.
	*/
	TwoRoundTrailCoreStack(const KeccakFPropagation& aDCorLC, StateAsVectorOfSlices& stateA, StateAsVectorOfSlices& stateB, unsigned int& aW0, unsigned int& aW1, bool& aComplete, unsigned int& aRootPeriod);

	/**
	* This method pushes an orbital to the cache.
	* @param aOrbital the orbital to be pushed.
	*/
	void push(const OrbitalPosition& aOrbital);

	/**
	* This method pushes a column to the cache.
	* @param aColumn the column to be pushed.
	*/
	void push(const Column& aColumn);

	/**
	* This method performs a dummy push to the cache.
	*/
	void pushDummy();

	/**
	* This method pops the highest unit from the cache when it is a column.
	* @param aColumn the column to pop.
	*/
	void pop(const Column& aColumn);

	/**
	* This method pops the highest unit from the cache when it is an orbital.
	* @param aOrbital the orbital to pop.
	*/
	void pop(const OrbitalPosition& aOrbital);

	/** This method outputs the cache to save it in, e.g., a file.
	* @param fout The stream to save the cache to.
	*/
	void save(ostream& fout);

	/** This methods returns the trail corresponding to the node.
	* @return The current trail.
	*/
	Trail getTrail();

	/** This methods returns the value of a row in position (y,z) in A.
	* @return The row value.
	*/
	RowValue getRowA(unsigned int y, unsigned int z) const;

	/** This methods returns the value of a row in position (y,z) in B.
	* @return The row value.
	*/
	RowValue getRowB(unsigned int y, unsigned int z) const;

private:
	// method that pushes an unaffected odd column
	void pushUnaffectedOddColumn(const Column& aColumn);
	// method that pushes an affected even column
	void pushAffectedEvenColumn(const Column& aColumn);
	// method that sets a bit to one and returns the difference of minimum reverse weight before and after the addition
	int pushBitAndGetDeltaMinReverseWeight(vector<SliceValue>& state, const BitPosition& p);
	// method that sets a bit to one and returns the difference of weight before and after the addition
	int pushBitAndGetDeltaWeight(vector<SliceValue>& state, const BitPosition& p);
};

/**
* \brief TwoRoundTrailCoreCostFunction class : the cost of a 2-round trail core in Keccak-f.
*
* This class represents the cost function for 2-round trail cores in Keccak-f.
* The cost function is defined as: \alpha*w0+\beta*w1.
* Examples are: w0, w1, w0+w1, w0+2w1, 2w0+w1.
*/
class TwoRoundTrailCoreCostFunction {

public:
	/** The weight for the minimum reverse weight. */
	unsigned int alpha;
	/** The weight for the weight. */
	unsigned int beta;

public:

	/** The default constructor. */
	TwoRoundTrailCoreCostFunction()
		: alpha(1), beta(1){}

	/** The constructor.
	* @param aAlpha The weight for the minimum reverse weight.
	* @param aBeta The weight for the weight.
	*/
	TwoRoundTrailCoreCostFunction(unsigned int aAlpha, unsigned int aBeta)
		: alpha(aAlpha), beta(aBeta) {}

	/** This method returns the cost of a given node.
	* @param unitList The list of units representing the given node.
	* @param cache The cache representation of the given node.
	* @return The cost of the given node.
	*/
	unsigned int getCost(const vector<OrbitalPosition>& unitList, const TwoRoundTrailCoreStack& cache) const;

	/** This method checks whether the addition of a new orbital to a given node is within a given budget.
	* @param unitList The list of units representing the given node.
	* @param cache The cache representation of the given node.
	* @newOrbital The orbital to be added.
	* @maxCost The given budget.
	* @return The cost of the given node.
	*/
	bool canAfford(const vector<OrbitalPosition>& unitList, const TwoRoundTrailCoreStack& cache, const OrbitalPosition& newOrbital, unsigned int maxCost, vector<unsigned int> &cost) const;
};


/**
* \brief TwoRoundTrailCoreCostBoundFunction class : the cost bound of a 2-round trail core in Keccak-f.
*
* This class represents a bound on the cost of a 2-round trail core in Keccak-f and of its children.
*/
class TwoRoundTrailCoreCostBoundFunction : public TwoRoundTrailCoreCostFunction {

public:

	/** The default constructor. */
	TwoRoundTrailCoreCostBoundFunction()
		: TwoRoundTrailCoreCostFunction(){};

	/** The constructor.
	* @param aAlpha The weight for the minimum reverse weight.
	* @param aBeta The weight for the weight.
	*/
	TwoRoundTrailCoreCostBoundFunction(unsigned int aAlpha, unsigned int aBeta)
		: TwoRoundTrailCoreCostFunction(aAlpha, aBeta){};

	/** This method returns the cost bound of a given node nd of its children.
	* @param unitList The list of units representing the given node.
	* @param cache The cache representation of the given node.
	* @return The cost bound of the given node.
	*/
	unsigned int getCost(const vector<Column>& unitList, const TwoRoundTrailCoreStack& cache) const;

	/** This method checks whether the addition of a new column to a given node is within a given budget.
	* @param unitList The list of units representing the given node.
	* @param cache The cache representation of the given node.
	* @newColumn The column to be added.
	* @maxCost The given budget.
	* @return The cost bound.
	*/
	bool canAfford(const vector<Column>& unitList, const TwoRoundTrailCoreStack& cache, Column& newColumn, unsigned int maxCost, const vector<unsigned int> &cost) const;
};


/**
* \brief TwoRoundTrailCore class : output representation for 2-round trail cores in Keccak-f.
* \details This class represents a 2-round trail core in Keccak-f.
*/
class TwoRoundTrailCore{

public:
	/** The state at A. */
	StateAsVectorOfSlices stateA;
	/** The state at B. */
	StateAsVectorOfSlices stateB;
	/** The minimum reverse weight of A. */
	unsigned int w0;
	/** The weight of B. */
	unsigned int w1;
	/** The 2-round trail core. */
	Trail trail;
	/** The completness of the trail. */
	bool complete;
	/** The z-period of the trail. */
	unsigned int zPeriod;
	/** The parity pattern of the state before theta. */
	vector<RowValue> C;
	/** The theta-effect of the state before theta. */
	vector<RowValue> D;

public:

	/** the default constructor. */
	TwoRoundTrailCore(){};

	/** This methods outputs the 2-round trail core to save it in, e.g., a file.
	* @param fout The stream to save the trail core to.
	*/
	void save(ostream& fout);

	/** This methods sets the trail.
	* @param unitList The list of orbitals representing the trail.
	* @param cache The cache representation of the trail.
	*/
	void set(const std::vector<OrbitalPosition>& unitList, const TwoRoundTrailCoreStack& cache);

	/** This methods sets the trail.
	* @param unitList The list of columns representing the trail.
	* @param cache The cache representation of the trail.
	*/
	void set(const std::vector<Column>& unitList, const TwoRoundTrailCoreStack& cache);
};


/** The orbital tree. */
typedef GenericTreeIterator<OrbitalPosition, OrbitalsSet, TwoRoundTrailCoreStack, TwoRoundTrailCore, TwoRoundTrailCoreCostFunction> OrbitalTreeIterator;

/** The run tree. */
typedef GenericTreeIterator<Column, ColumnsSet, TwoRoundTrailCoreStack, TwoRoundTrailCore, TwoRoundTrailCoreCostBoundFunction> RunTreeIterator;

#endif
