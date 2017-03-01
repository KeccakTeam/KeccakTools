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

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "Keccak-fTree.h"

OrbitalPosition OrbitalsSet::getFirstChildUnit(const std::vector<OrbitalPosition>& unitList, const TwoRoundTrailCoreStack& cache) const
{
    (void)cache;
	OrbitalPosition newOrbital;
	if (unitList.empty()) {
		if (!newOrbital.first(yMin, laneSize))
			throw EndOfSet();
	}
	else {
		if (!newOrbital.successorOf(unitList.back(), yMin, laneSize))
			throw EndOfSet();
	}
	return newOrbital;
}

void OrbitalsSet::iterateUnit(const std::vector<OrbitalPosition>& unitList, OrbitalPosition& current, const TwoRoundTrailCoreStack& cache) const
{
	(void)unitList;
    (void)cache;
	if (!current.next(yMin, laneSize))
		throw EndOfSet();
}

unsigned int OrbitalsSet::compare(const OrbitalPosition& first, const OrbitalPosition& second) const
{

	if (first.z < second.z)
		return 1;
	else if (first.z == second.z){
		if (first.x < second.x)
			return 1;
		else if (first.x == second.x){
			if (first.y0 < second.y0)
				return 1;
			else if (first.y0 == second.y0){
				if (first.y1 < second.y1)
					return 1;
				else if (first.y1 == second.y1){
					return 0;
				}
			}
		}
	}
	return 2;
}

bool OrbitalsSet::isCanonical(const std::vector<OrbitalPosition>& orbitalList, TwoRoundTrailCoreStack& cache) const
{
	cache.nodePeriod = laneSize;

	if (kernel){
		if (orbitalList[0].z != 0)
			return false;

		unsigned int lastZ = 0;
		for (unsigned int i = 0; i < orbitalList.size(); i++){
			unsigned int z = orbitalList[i].z;
			if (z != 0 && z > lastZ){ // consider translation by z only if it has not been already considered before
				lastZ = z;
				// translate by z
				vector<OrbitalPosition> tauList; // translated list
				for (unsigned int j = i; j < orbitalList.size(); j++){
					OrbitalPosition orbital = orbitalList[j];
					orbital.z -= z;
					tauList.push_back(orbital);
				}
				for (unsigned int j = 0; j < i; j++){
					OrbitalPosition orbital = orbitalList[j];
					orbital.z = orbital.z - z + laneSize;
					tauList.push_back(orbital);
				}
				// compare lists
				unsigned int k = 0;
				for (k = 0; k < orbitalList.size(); k++){
					unsigned int cmp = compare(tauList[k], orbitalList[k]);
					if (cmp == 1)
						return false; // there is a traslated variant smaller than the original
					if (cmp == 2)
						break; // original list is smaller than this translated variant, but still need to check other variants
				}
				// if the two list are identical, then the list is z-periodic.
				// No need to check the other translated variants.
				if (k == orbitalList.size()){
					cache.nodePeriod = z;
					break;
				}
			}
		}
		return true;
	}

	else{
		if (cache.rootPeriod == laneSize)
			return true;
		// in the case outside the kernel, only traslations by the period must be considered
		for (unsigned int z = cache.rootPeriod; z < laneSize; z+=cache.rootPeriod){
			// traslate by z
			vector<OrbitalPosition> tauList; // translated list
			unsigned int i = 0;
			for (i = 0; i < orbitalList.size(); i++){
				if (orbitalList[i].z >= z)
					break;
			}
			for (unsigned int j = i; j < orbitalList.size(); j++){
				OrbitalPosition orbital = orbitalList[j];
				orbital.z -= z;
				tauList.push_back(orbital);
			}
			for (unsigned int j = 0; j < i; j++){
				OrbitalPosition orbital = orbitalList[j];
				orbital.z = orbital.z - z + laneSize;
				tauList.push_back(orbital);
			}
			// compare lists
			unsigned int k = 0;
			for (k = 0; k < orbitalList.size(); k++){
				unsigned int cmp = compare(tauList[k], orbitalList[k]);
				if (cmp == 1)
					return false; // there is a traslated variant smaller than the original
				if (cmp == 2)
					break; // original list is smaller than this translated variant, but still need to check other variants
			}
			// if the two list are identical, then the list is z-periodic
			// no need to check the other translated variants
			if (k == orbitalList.size()){
				cache.nodePeriod = z;
				break;
			}
		}
		return true;
	}
}


Column::Column(){
	position.x = 0;
	position.z = 0;
	value = 0;
	index = 0;
	odd = false;
	affected = false;
	entangled = false;
	starting = false;
}

Column::Column(bool& aOdd, bool& aAffected){

	position.x = 0;
	position.z = 0;
	value = 0;
	index = 0;
	odd = aOdd;
	affected = aAffected;
	entangled = false;
	starting = false;

}

const ColumnValue ColumnsSet::UOValues[5] = {
	0x01, 0x02, 0x04, 0x08, 0x10 };

const ColumnValue ColumnsSet::AEValues[16] = {
	0x00, 0x03, 0x05, 0x06, 0x09, 0x0A, 0x0C, 0x0F,
	0x11, 0x12, 0x14, 0x17, 0x18, 0x1B, 0x1D, 0x1E };

const ColumnValue ColumnsSet::AOValues[16] = {
	0x01, 0x02, 0x04, 0x07, 0x08, 0x0B, 0x0D, 0x0E,
	0x10, 0x13, 0x15, 0x16, 0x19, 0x1A, 0x1C, 0x1F };

bool ColumnsSet::checkColumnOverlapping(const std::vector<Column>& unitList, Column& current, const TwoRoundTrailCoreStack& cache) const
{
    (void)unitList;
	current.entangled = false;
	//unsigned int i = 0;
	
	// AEC cannot overlap another AEC
	// and can overlap only an odd-0 column
			if (current.affected){
		// if already affected
		if ((getBit(cache.D, current.position.x, current.position.z) & 1))
			return true;
		// if odd-0
		if ((getBit(cache.C, current.position.x, current.position.z) & 1)){
			if (getColumn(cache.stack_stateAtA.top(), current.position.x, current.position.z) != 1)
					return true;
				else{
					current.entangled = true;
				return false;
				}
			}
		return false;
	}

	// if UOC cannot overlap another UOC
	//if (current.odd){
    else{    
		if ((getBit(cache.C, current.position.x, current.position.z) & 1))
					return true;
		if ((getBit(cache.D, current.position.x, current.position.z) & 1)){
					current.entangled = true;
			return false;
			}
		return false;
		}
}

Column ColumnsSet::getFirstChildUnit(const std::vector<Column>& unitList, const TwoRoundTrailCoreStack& cache) const
{
	Column newColumn;
	if (unitList.empty()) {
		newColumn.position.x = 0;
		newColumn.position.z = 0;
		newColumn.odd = false;
		newColumn.affected = true;
		newColumn.value = AEValues[0];
		newColumn.entangled = false;
		newColumn.starting = true;
	}
	else{
		// there are different cases:
		// UOC -> successor is the AEC
		// AEC -> if it's the starting AEC then successor is UOC
		//	   -> if it's the ending AEC then successor is another AEC
		// AOC -> never treated, because it's the sum of UOC and AEC
		// each time check if it is overlapping another run

		// case UOC -> new column is ending AEC
		if (!unitList.back().affected && unitList.back().odd){
			newColumn.affected = true;
			newColumn.odd = false;
			newColumn.position.x = (unitList.back().position.x + 4) % 5;
			newColumn.position.z = (unitList.back().position.z + 1) % laneSize;
			newColumn.value = AEValues[0];
			newColumn.starting = false;
			// if it overlaps another AEC then it's the same run
			if (checkColumnOverlapping(unitList, newColumn, cache))
				throw EndOfSet();
		}
		// case AEC
		else if (unitList.back().affected && !unitList.back().odd){
			// case starting AEC -> new column is UOC
			if (unitList.back().starting){
				newColumn.affected = false;
				newColumn.odd = true;
				newColumn.position.x = (unitList.back().position.x + 4) % 5;
				newColumn.position.z = (unitList.back().position.z);
				newColumn.value = UOValues[0];
				newColumn.starting = false;

				if (checkColumnOverlapping(unitList, newColumn, cache))
					throw EndOfSet();
			}
			// case ending AEC -> new column is a starting AEC
			else {
				// the position of the new starting AEC is the position of the new run, that depends on the starting position
				// of the last run. Namely, go back and look for the previous (starting) AEC.
				// for z-canonicity, a new AEC cannot have x-coordinate smaller than x-coordinate of first AEC.
				unsigned int i = 3; // size-1 is ending AEC, size-2 is UOC, check from size-3 on.
				for (i = 3; i <= unitList.size(); i++){
					if (unitList[unitList.size() - i].starting == true)
						break;
				}
				newColumn.affected = true;
				newColumn.odd = false;
				newColumn.position.x = unitList[unitList.size() - i].position.x;
				newColumn.position.z = unitList[unitList.size() - i].position.z;
				newColumn.starting = true;
				do {
					if (newColumn.position.x < 4){
						newColumn.position.x = (newColumn.position.x + 1) % 5;
						newColumn.value = AEValues[0];
						newColumn.index = 0;
					}
					else {
						if (newColumn.position.z < laneSize - 1){
							// optimization for z-canonicity: a starting AEC translated by its z cannot come before the first AEC,
							// otherwise the pattern will not be canonical
							newColumn.position.x = unitList[0].position.x;
							newColumn.position.z = (newColumn.position.z + 1) % laneSize;
							newColumn.index = unitList[0].index;
							newColumn.value = AEValues[newColumn.index];
						}
						else
							throw EndOfSet();
					}
					if (!checkColumnOverlapping(unitList, newColumn, cache))
						break;
				} while (true);
			}
		}
	}
	return newColumn;
}

void ColumnsSet::iterateUnit(const std::vector<Column>& unitList, Column& current, const TwoRoundTrailCoreStack& cache) const
{
	// UOC -> iterate column value. It cannot change position.
	// ending AEC -> iterate column value
	//            -> transform it in a UOC to continue the same run
	// starting AEC -> iterate column value
	//              -> change its position to start the run in a new position

	// UOC
	if (!current.affected && current.odd){
		// iterate column value
		// if it is odd-0 the active bit is only in y=0
		if (current.entangled)
			throw EndOfSet();
		if (current.index < 4){
			current.index++;
			current.value = UOValues[current.index];
		}
		else
			throw EndOfSet();
	}

	// AEC
	else{
		// case starting AEC
		if (current.starting){
			// try to iterate column value
			if (current.index < 15){
				current.index++;
				current.value = AEValues[current.index];
			}
			// else change position to start a new run
			else {
				do {
					if (current.position.x < 4){
						current.position.x = (current.position.x + 1) % 5;
						current.index = 0;
						current.value = AEValues[0];
					}
					else {
						// if it is the first AEC then it is restricted to the first slice
						if (unitList.size()==0)
							throw EndOfSet();
						else if (current.position.z < laneSize - 1){
							// optimization for z-canonicity: a starting AEC translated by its z cannot come before the first AEC,
							// otherwise the pattern will not be canonical
							current.position.x = unitList[0].position.x;
							current.position.z = (current.position.z + 1) % laneSize;
							current.index = unitList[0].index;
							current.value = AEValues[current.index];
						}
						else
							throw EndOfSet();
					}
					if (!checkColumnOverlapping(unitList, current, cache))
						break;
				} while (true);
			}
		}

		// case ending AEC
		else {
			// try to iterate column value
			if (current.index < 15){
				current.index++;
				current.value = AEValues[current.index];
			}
			// transform it in a UOC to continue the same run
			else {
				current.index = 0;
				current.value = UOValues[0];
				current.affected = false;
				current.odd = true;
				current.position.x = (current.position.x + 4) % 5;
				if (checkColumnOverlapping(unitList, current, cache))
					throw EndOfSet();
			}
		}
	}
	if (unitList.empty() && current.position.z>0)
		throw EndOfSet();
}

unsigned int ColumnsSet::compare(const Column& first, const Column& second) const
{
	// shorter run comes before longer run
	if ((!first.starting && first.affected) && second.odd)
		return 1;
	if (first.odd && (!second.starting && second.affected))
		return 2;
	else {
		if (first.position.z < second.position.z)
			return 1;
		else if (first.position.z == second.position.z){
			if (first.position.x < second.position.x)
				return 1;
			else if (first.position.x == second.position.x){
				if (first.value < second.value)
					return 1;
				else if (first.value == second.value){
						return 0;
				}
			}
		}
	}
	return 2;
}

bool ColumnsSet::isCanonical(const std::vector<Column>& unitList, TwoRoundTrailCoreStack& cache) const
{
	cache.nodePeriod = laneSize;

	// optimization for z-canonicity: the first unit is restricted to first slice
	if (unitList[0].position.z != 0)
		return false;

	// z-canonicity is checked only when the highest run is complete
	// namely, when highest unit is ending AEC
	if (unitList.back().odd)
		return true;
	if (unitList.back().starting){
		return true;
	}

	unsigned int lastZ = 0;

	for (unsigned int i = 0; i < unitList.size(); i++){
		if (unitList[i].starting){
			unsigned int z = unitList[i].position.z;
			if (z != 0 && z > lastZ){ // consider translation by z only if it has not been already considered before
				lastZ = z;
				// traslate by z
				vector<Column> tauList; // translated list
				for (unsigned int j = i; j < unitList.size(); j++){
					Column column = unitList[j];
					column.position.z = (column.position.z -z)%laneSize;
					tauList.push_back(column);
				}
				for (unsigned int j = 0; j < i; j++){
					Column column = unitList[j];
					column.position.z = (column.position.z - z) % laneSize;
					tauList.push_back(column);
				}
				// compare lists
				unsigned int j = 0;
				for (j = 0; j < unitList.size(); j++){
					unsigned int cmp = compare(tauList[j], unitList[j]);
					if (cmp == 1)
						return false; // there is a traslated variant smaller than the original
					if (cmp == 2)
						break; // original list is smaller than this translated variant, but still need to check other variants

				}
				// as soon as the original list is equal to a translated variant, then the unit list is z-periodic
				// no need to check the other variants.
				if (j == unitList.size()){
					cache.nodePeriod = z;
					break;
				}
			}
		}
	}
	return true;
}


TwoRoundTrailCoreStack::TwoRoundTrailCoreStack(const KeccakFPropagation& aDCorLC) :
DCorLC(aDCorLC) {

	laneSize = aDCorLC.laneSize;
	StateAsVectorOfSlices emptyState(laneSize, 0);
	stack_stateAtA.push(emptyState);
	stack_stateAtB.push(emptyState);
	stack_w0.push(0);
	stack_w1.push(0);
	stack_complete.push(true);
	rootPeriod = 0;
	nodePeriod = aDCorLC.laneSize;
	vector<RowValue> emptyPlane(laneSize, 0);
	C = emptyPlane;
	D = emptyPlane;
}

TwoRoundTrailCoreStack::TwoRoundTrailCoreStack(const KeccakFPropagation& aDCorLC, StateAsVectorOfSlices& stateA, StateAsVectorOfSlices& stateB, unsigned int& aW0, unsigned int& aW1, bool& aComplete, unsigned int& aRootPeriod) :
DCorLC(aDCorLC) {

	(void)aComplete;

	laneSize = aDCorLC.laneSize;
	stack_stateAtA.push(stateA);
	stack_stateAtB.push(stateB);
	stack_w0.push(aW0);
	stack_w1.push(aW1);
	stack_complete.push(true);
	rootPeriod = aRootPeriod;
	nodePeriod = aDCorLC.laneSize;
	getParity(stateA, C);
	aDCorLC.directThetaEffectFromParities(C, D);
}


void TwoRoundTrailCoreStack::push(const OrbitalPosition& aOrbital){

	unsigned int w0 = stack_w0.top();
	unsigned int w1 = stack_w1.top();

	BitPosition p1(aOrbital.x, aOrbital.y0, aOrbital.z);
	BitPosition p2(aOrbital.x, aOrbital.y1, aOrbital.z);

	w0 -= DCorLC.getMinReverseWeightRow(getRow(stack_stateAtA.top(), p1.y, p1.z)) + DCorLC.getMinReverseWeightRow(getRow(stack_stateAtA.top(), p2.y, p2.z));
	setBitToOne(stack_stateAtA.top(), p1);
	setBitToOne(stack_stateAtA.top(), p2);
	w0 += DCorLC.getMinReverseWeightRow(getRow(stack_stateAtA.top(), p1.y, p1.z)) + DCorLC.getMinReverseWeightRow(getRow(stack_stateAtA.top(), p2.y, p2.z));

	DCorLC.directRhoPi(p1);
	DCorLC.directRhoPi(p2);

	w1 -= DCorLC.getWeightRow(getRow(stack_stateAtB.top(), p1.y, p1.z)) + DCorLC.getWeightRow(getRow(stack_stateAtB.top(), p2.y, p2.z));
	setBitToOne(stack_stateAtB.top(), p1);
	setBitToOne(stack_stateAtB.top(), p2);
	w1 += DCorLC.getWeightRow(getRow(stack_stateAtB.top(), p1.y, p1.z)) + DCorLC.getWeightRow(getRow(stack_stateAtB.top(), p2.y, p2.z));

	stack_w0.push(w0);
	stack_w1.push(w1);

	stack_complete.push(true);
	dummy = false;
}

void TwoRoundTrailCoreStack::push(const Column& aColumn){

	if (aColumn.affected && !aColumn.odd){
		pushAffectedEvenColumn(aColumn);
		setBitToOne(D, aColumn.position.x, aColumn.position.z);
		// if ending AEC then complete is true
		if (aColumn.starting) {
			stack_complete.push(false);
		}
		else
			stack_complete.push(true);
	}
	else if (!aColumn.affected && aColumn.odd){
		pushUnaffectedOddColumn(aColumn);
		setBitToOne(C, aColumn.position.x, aColumn.position.z);
		stack_complete.push(false);
	}
	dummy = false;
}

void TwoRoundTrailCoreStack::pushDummy()
{
	stack_w0.push(stack_w0.top());
	stack_w1.push(stack_w1.top());
	stack_complete.push(true);
	dummy = true;
}

void TwoRoundTrailCoreStack::pushUnaffectedOddColumn(const Column& aColumn)
{
	// the active bit is pushed in A and B
	// the weights are updated

	int delta0 = 0, delta1 = 0;

	for (unsigned int y = 0; y<5; y++) {
		bool activeBit = ((aColumn.value >> y) & 1) != 0;
		if (activeBit) {

			BitPosition p(aColumn.position.x, y, aColumn.position.z);

			DCorLC.reverseRhoPiBeforeTheta(p);
			delta0 += pushBitAndGetDeltaMinReverseWeight(stack_stateAtA.top(), p);

			DCorLC.directRhoPiAfterTheta(p);
			delta1 += pushBitAndGetDeltaWeight(stack_stateAtB.top(), p);

			break;
		}
	}

	int new_w0 = delta0 + stack_w0.top();
	int new_w1 = delta1 + stack_w1.top();

	stack_w0.push(new_w0);
	stack_w1.push(new_w1);
}

void TwoRoundTrailCoreStack::pushAffectedEvenColumn(const Column& aColumn)
{
	// for each bit of the column:
	// - if 0: is pushed in B
	// - if 1: is pushed in A
	// the weights are updated

	int delta0 = 0, delta1 = 0;

	for (unsigned int y = 0; y<5; y++) {
		bool bitBeforeTheta = ((aColumn.value >> y) & 1) != 0;
		BitPosition p(aColumn.position.x, y, aColumn.position.z);
		if (bitBeforeTheta) {
			DCorLC.reverseRhoPiBeforeTheta(p);
			delta0 += pushBitAndGetDeltaMinReverseWeight(stack_stateAtA.top(), p);
		}
		else {
			DCorLC.directRhoPiAfterTheta(p);
			delta1 += pushBitAndGetDeltaWeight(stack_stateAtB.top(), p);
		}
	}
	int new_w0 = delta0 + stack_w0.top();
	int new_w1 = delta1 + stack_w1.top();

	stack_w0.push(new_w0);
	stack_w1.push(new_w1);
}

void TwoRoundTrailCoreStack::pop(const Column& aColumn){

	stack_w0.pop();
	stack_w1.pop();
	stack_complete.pop();

	if (dummy){
		dummy = false;
		return;
	}

	if (aColumn.odd){
		setBitToZero(C, aColumn.position.x, aColumn.position.z);
		for (unsigned int y = 0; y<5; y++) {
			bool activeBit = ((aColumn.value >> y) & 1) != 0;
			if (activeBit) {

				BitPosition p(aColumn.position.x, y, aColumn.position.z);

				DCorLC.reverseRhoPiBeforeTheta(p);
				if ((getBit(stack_stateAtA.top(), p) & 1) == 0)
					setBitToOne(stack_stateAtA.top(), p);
				else
					setBitToZero(stack_stateAtA.top(), p);

				DCorLC.directRhoPiAfterTheta(p);
				if ((getBit(stack_stateAtB.top(), p) & 1) == 0)
					setBitToOne(stack_stateAtB.top(), p);
				else
					setBitToZero(stack_stateAtB.top(), p);
				break;
			}
		}
	}

	else{
		setBitToZero(D, aColumn.position.x, aColumn.position.z);
		for (unsigned int y = 0; y<5; y++) {
			bool bitBeforeTheta = ((aColumn.value >> y) & 1) != 0;
			BitPosition p(aColumn.position.x, y, aColumn.position.z);
			if (bitBeforeTheta) {
				DCorLC.reverseRhoPiBeforeTheta(p);
				if ((getBit(stack_stateAtA.top(), p) & 1) == 0)
					setBitToOne(stack_stateAtA.top(), p);
				else
					setBitToZero(stack_stateAtA.top(), p);
			}
			else {
				DCorLC.directRhoPiAfterTheta(p);
				if ((getBit(stack_stateAtB.top(), p) & 1) == 0)
					setBitToOne(stack_stateAtB.top(), p);
				else
					setBitToZero(stack_stateAtB.top(), p);
			}
		}
	}
}

void TwoRoundTrailCoreStack::pop(const OrbitalPosition& aOrbital){

	stack_w0.pop();
	stack_w1.pop();
	stack_complete.pop();

	if (dummy){
		dummy = false;
		return;
	}

	BitPosition p1(aOrbital.x, aOrbital.y0, aOrbital.z);
	BitPosition p2(aOrbital.x, aOrbital.y1, aOrbital.z);

	setBitToZero(stack_stateAtA.top(), p1);
	setBitToZero(stack_stateAtA.top(), p2);

	DCorLC.directRhoPi(p1);
	DCorLC.directRhoPi(p2);

	setBitToZero(stack_stateAtB.top(), p1);
	setBitToZero(stack_stateAtB.top(), p2);

}


void TwoRoundTrailCoreStack::save(ostream& fout){

	getTrail().save(fout);

}

Trail TwoRoundTrailCoreStack::getTrail(){

	Trail trail;
	trail.setFirstStateReverseMinimumWeight(stack_w0.top());
	trail.append(stack_stateAtB.top(), stack_w1.top());
	return trail;

}

RowValue TwoRoundTrailCoreStack::getRowA(unsigned int y, unsigned int z) const
{
	return getRow(stack_stateAtA.top(), y, z);
}

RowValue TwoRoundTrailCoreStack::getRowB(unsigned int y, unsigned int z) const
{

	return getRow(stack_stateAtB.top(), y, z);
}

int TwoRoundTrailCoreStack::pushBitAndGetDeltaMinReverseWeight(vector<SliceValue>& state, const BitPosition& p)
{
	int weightBefore = DCorLC.getMinReverseWeight(state[p.z]);
	if ((getBit(state, p) & 1) == 0)
		setBitToOne(state, p);
	else
		setBitToZero(state, p);
	return (int)DCorLC.getMinReverseWeight(state[p.z]) - weightBefore;
}

int TwoRoundTrailCoreStack::pushBitAndGetDeltaWeight(vector<SliceValue>& state, const BitPosition& p)
{
	int weightBefore = DCorLC.getWeight(state[p.z]);
	if ((getBit(state, p) & 1) == 0)
		setBitToOne(state, p);
	else
		setBitToZero(state, p);
	return (int)DCorLC.getWeight(state[p.z]) - weightBefore;
}


unsigned int TwoRoundTrailCoreCostFunction::getCost(const vector<OrbitalPosition>& unitList, const TwoRoundTrailCoreStack& cache) const
{
	(void) unitList;
	return alpha*cache.stack_w0.top() + beta*cache.stack_w1.top();
}

bool TwoRoundTrailCoreCostFunction::canAfford(const vector<OrbitalPosition>& unitList, const TwoRoundTrailCoreStack& cache, const OrbitalPosition& newOrbital, unsigned int maxCost, vector<unsigned int> &cost) const
{
	(void) unitList;
	(void) cost;
	unsigned int Gamma = alpha*cache.stack_w0.top() + beta*cache.stack_w1.top();

	if (alpha != 0) {
		// check y_bottom in A
		if (cache.getRowA(newOrbital.y0, newOrbital.z) == 0) {
			Gamma += alpha * 2;
			if (Gamma>maxCost)
				return false;
		}
		// check y_top in A
		if (cache.getRowA(newOrbital.y1, newOrbital.z) == 0) {
			Gamma += alpha * 2;
			if (Gamma>maxCost)
				return false;
	}
	}
	if (beta != 0){
		// check y_bottom in B
		BitPosition p_b(newOrbital.x, newOrbital.y0, newOrbital.z);
		cache.DCorLC.directRhoPi(p_b);
		if (cache.getRowB(p_b.y, p_b.z) == 0) {
			Gamma += beta * 2;
			if (Gamma>maxCost)
				return false;
		}
		// check y_top in B
		BitPosition p_t(newOrbital.x, newOrbital.y1, newOrbital.z);
		cache.DCorLC.directRhoPi(p_t);
		if (cache.getRowB(p_t.y, p_t.z) == 0) {
			Gamma += beta * 2;
			if (Gamma>maxCost)
				return false;
		}
	}
	return true;
}

void TwoRoundTrailCore::set(const std::vector<OrbitalPosition>& unitList, const TwoRoundTrailCoreStack& cache)
{
	(void) unitList;
	stateA = cache.stack_stateAtA.top();
	stateB = cache.stack_stateAtB.top();
	C = cache.C;
	D = cache.D;
	w0 = cache.stack_w0.top();
	w1 = cache.stack_w1.top();
	complete = cache.stack_complete.top();
	zPeriod = cache.nodePeriod;
	trail.clear();
	trail.setFirstStateReverseMinimumWeight(w0);
	trail.append(stateB, w1);
}

void TwoRoundTrailCore::set(const std::vector<Column>& unitList, const TwoRoundTrailCoreStack& cache)
{
	(void) unitList;
	stateA = cache.stack_stateAtA.top();
	stateB = cache.stack_stateAtB.top();
	C = cache.C;
	D = cache.D;
	w0 = cache.stack_w0.top();
	w1 = cache.stack_w1.top();
	complete = cache.stack_complete.top();
	zPeriod = cache.nodePeriod;
	trail.clear();
	trail.setFirstStateReverseMinimumWeight(w0);
	trail.append(stateB, w1);
}

void TwoRoundTrailCore::save(ostream& fout)
{
	trail.save(fout);
}


unsigned int TwoRoundTrailCoreCostBoundFunction::getCost(const vector<Column>& unitList, const TwoRoundTrailCoreStack& cache) const
{

	//Remove from all affected even columns the active bit in y = 0 at a or in the corresponding position at b.
	//For each odd-0 column in slices with only y=0 active row, remove bit at b.
	//Compute the cost of the resulting state pattern.
	//For cost function alpha*wrev(a) + beta*w(b), subtract 2max(alpha,beta)n0 from that cost, with n0
	//the total number of odd-0 columns in slices with row y!=0 active.

	// remove y=0 from all AEC
	StateAsVectorOfSlices stateA = cache.stack_stateAtA.top();
	StateAsVectorOfSlices stateB = cache.stack_stateAtB.top();

	for (unsigned int i = 0; i < unitList.size(); i++){
		if (unitList[i].affected == true){
			if((getBit(cache.C, unitList[i].position.x, unitList[i].position.z) & 1) == 0){
				BitPosition p(unitList[i].position.x, 0, unitList[i].position.z);
				if ( (getBit(stateA, p) & 1)){
					cache.DCorLC.reverseRhoPiBeforeTheta(p);
					setBitToZero(stateA, p);
				}
				else{
					cache.DCorLC.directRhoPiAfterTheta(p);
					setBitToZero(stateB, p);
				}
			}
		}
	}

	// compute the number of odd-0 columns and remove bit after lambda if needed
	unsigned int n0 = 0;
	for (unsigned int i = 0; i < unitList.size(); i++){
		if (unitList[i].odd == true){
			if (!(getBit(cache.D, unitList[i].position.x, unitList[i].position.z) & 1)){
				if (unitList[i].index==0){
					// set this bit to zero and check if the slice value is zero
					SliceValue slice = stateA[unitList[i].position.z];
					slice &= ~((SliceValue)1 << (unitList[i].position.x + 5 * 0));
					if (slice != 0){
						// check if y=0 is the only active row
						if ((getRowFromSlice(slice, 1) == 0) && (getRowFromSlice(slice, 2) == 0) && (getRowFromSlice(slice, 3) == 0) && (getRowFromSlice(slice, 4) == 0)){
							// deactivate its bit after lambda
							BitPosition p(unitList[i].position.x, 0, unitList[i].position.z);
							cache.DCorLC.directRhoPiAfterTheta(p);
							setBitToZero(stateB, p);
						}
						else
					n0++;
					}
					
					else{
						// deactivate its bit after lambda
						BitPosition p(unitList[i].position.x, 0, unitList[i].position.z);
						cache.DCorLC.directRhoPiAfterTheta(p);
						setBitToZero(stateB, p);
					}
				}
			}
		}
	}

	// compute the new cost
	int weight = alpha*cache.DCorLC.getMinReverseWeight(stateA) + beta*cache.DCorLC.getWeight(stateB);

	return max(0,weight-int(2*max(alpha,beta)*n0));
}

bool TwoRoundTrailCoreCostBoundFunction::canAfford(const vector<Column>& unitList, const TwoRoundTrailCoreStack& cache, Column& newColumn, unsigned int maxCost, const vector<unsigned int> &cost) const
{

	if (newColumn.odd){

		unsigned int currentCost = 0;
		if (!cost.empty())
			currentCost = cost.back();

		// the maximum possible contribution is 2*alpha + 2*beta.
		// can avoid the check if the budget is bigger
		if ((maxCost - currentCost) > (2 * alpha + 2 * beta))
			return true;

		// distinguish between two cases: whether the slice is empty or not
		SliceValue slice = cache.stack_stateAtA.top()[newColumn.position.z];

		// case empty slice
		if (slice == 0){
			// if the column is not affected, then the odd column will contribute 2*alpha to w0 wherever its active bit is
			if (!(getBit(cache.D, newColumn.position.x, newColumn.position.z) & 1)){
				unsigned int newCost = currentCost + alpha * 2;
				if (newCost > maxCost){
					// if already above the budget, then this column position is not affordable whatever its value is
					newColumn.index = 4;
					return false;
				}
				// if it is an odd-0 column then the contribution at b is not relevant
				if (newColumn.index == 0)
					return true;
				// else consider the contribution of its bit at b
				// it will contribute at most 2*beta
				if (maxCost - newCost > 2 * beta)
					return true;
				// look for the active bit and add its corresponding at b
				for (unsigned int y = 1; y < 5; y++) {
					bool activeBit = ((newColumn.value >> y) & 1) != 0;
					if (activeBit) {
						BitPosition p(newColumn.position.x, y, newColumn.position.z);
						cache.DCorLC.directRhoPiAfterTheta(p);
						SliceValue sliceB = cache.stack_stateAtB.top()[p.z];
						unsigned int weightBefore = cache.DCorLC.getWeightRow(getRowFromSlice(sliceB, p.y));
						// set the bit to one
						sliceB |= (SliceValue)1 << (p.x + 5 * p.y);
						unsigned int weightAfter = cache.DCorLC.getWeightRow(getRowFromSlice(sliceB, p.y));
						newCost = newCost + beta*(weightAfter - weightBefore);

						if (newCost > maxCost)
							return false;
						else
							return true;
					}
				}
				return true;
			}
			// if the column is AEC, then newColumn can only be odd-0
			// and the corresponding bit at b is canceled
			// the contribution at a is alpha*2
			else{
				BitPosition p(newColumn.position.x, 0, newColumn.position.z);
				cache.DCorLC.directRhoPiAfterTheta(p);
				SliceValue sliceB = cache.stack_stateAtB.top()[p.z];
				int weightBefore = cache.DCorLC.getWeightRow(getRowFromSlice(sliceB, p.y));
				// set the bit to zero
				sliceB &= ~((SliceValue)1 << (p.x + 5 * p.y));
				int weightAfter = cache.DCorLC.getWeightRow(getRowFromSlice(sliceB, p.y));
				int newCost = currentCost + (int)beta*(weightAfter - weightBefore) + alpha*2;
				if (newCost > (int)maxCost)
					return false;
				else
					return true;
			}
		}

		// case non-empty slice
		else{
			// the contribution depends on the position of the active bit
			// if the column is already affected, then new column is odd-0 and may cancel bits and decrease the cost
			// distinguishing the cases is too expensive, return true and do the real push
			if ((getBit(cache.D, newColumn.position.x, newColumn.position.z) & 1))
				return true;
			// case unaffected column
			// if the column is not affected, then it is empty and the column will add new bits
			for (unsigned int y = 0; y < 5; y++) {
				bool activeBit = ((newColumn.value >> y) & 1) != 0;
				if (activeBit) {
					// if it is not odd-0 then both the contribution at a and b are relevant
					// compute the contribution at a
					BitPosition p(newColumn.position.x, y, newColumn.position.z);
					SliceValue sliceA = cache.stack_stateAtA.top()[p.z];
					unsigned int weightBefore = cache.DCorLC.getMinReverseWeightRow(getRowFromSlice(sliceA, p.y));
					// set the bit to one
					sliceA |= (SliceValue)1 << (p.x + 5 * p.y);
					unsigned int weightAfter = cache.DCorLC.getMinReverseWeightRow(getRowFromSlice(sliceA, p.y));
					unsigned int newCost = currentCost + alpha*(weightAfter - weightBefore);
					if (newCost > maxCost)
						return false;
					// if it is odd-0 column, then the contribution at b is not relevant
					if (y == 0)
						return true;
					// else check the contribution at b
					// it will contribute at most 2*beta
					if (maxCost - newCost > 2 * beta)
						return true;
					cache.DCorLC.directRhoPiAfterTheta(p);
					SliceValue sliceB = cache.stack_stateAtB.top()[p.z];
					weightBefore = cache.DCorLC.getWeightRow(getRowFromSlice(sliceB, p.y));
					// set the bit to one
					sliceB |= (SliceValue)1 << (p.x + 5 * p.y);
					weightAfter = cache.DCorLC.getWeightRow(getRowFromSlice(sliceB, p.y));
					newCost = currentCost + beta*(weightAfter - weightBefore);
					if (newCost > maxCost)
						return false;
					else
						return true;
				}
			}
		}
	}

	else if (newColumn.affected) {

		// Compute a lower bound for the cost of a new AEC, whatever its value will be.
		// An AEC will add 5 bits either at A or B, depending on its value.
		// Push all bits in A and B, there will be 5 bits too many.
		// Each of these bits adds at most a factor 2 to weight, that should not be added, so
		// compute new cost - max(alpha,beta)*10 (10 = 2*number of exceeding bits)

		// the contribution will be at most max(alpha,beta)*2*5
		// no need to check if the budget is bigger than that
		unsigned int currentCost = 0;
		if (!cost.empty())
			currentCost = cost.back();
		if ((maxCost - currentCost) > (max(alpha,beta)*10))
			return true;

		// push to A and B
		StateAsVectorOfSlices stateA = cache.stack_stateAtA.top();
		StateAsVectorOfSlices stateB = cache.stack_stateAtB.top();
		for (unsigned int y = 0; y < 5; y++) {
			BitPosition p(newColumn.position.x, y, newColumn.position.z);
			cache.DCorLC.reverseRhoPiBeforeTheta(p);
			setBitToOne(stateA, p.x, p.y, p.z);
			cache.DCorLC.directRhoPiAfterTheta(p);
			setBitToOne(stateB, p.x, p.y, p.z);
		}

		// remove y=0 from all AEC
		for (unsigned int i = 0; i < unitList.size(); i++){
			if (unitList[i].affected == true){
				if ((getBit(cache.C, unitList[i].position.x, unitList[i].position.z) & 1) == 0){
					BitPosition p(unitList[i].position.x, 0, unitList[i].position.z);
					if ((getBit(stateA, p) & 1)){
						cache.DCorLC.reverseRhoPiBeforeTheta(p);
						setBitToZero(stateA, p);
					}
					else{
			cache.DCorLC.directRhoPiAfterTheta(p);
						setBitToZero(stateB, p);
					}
				}
			}
		}

		// count the number of odd-0 
		unsigned int n0 = 0;

		for (unsigned int i = 0; i<unitList.size(); i++){
			if (unitList[i].odd == true){
				if ( !(getBit(cache.D, unitList[i].position.x, unitList[i].position.z) & 1)){
					if (unitList[i].index == 0){
						// set this bit to zero and check if the slice value is zero
						SliceValue slice = stateA[unitList[i].position.z];
						slice &= ~((SliceValue)1 << (unitList[i].position.x + 5 * 0));
						if (slice != 0){
							if ((getRowFromSlice(slice, 1) == 0) && (getRowFromSlice(slice, 2) == 0) && (getRowFromSlice(slice, 3) == 0) && (getRowFromSlice(slice, 4) == 0)){
								// deactivate its bit after lambda
								BitPosition p(unitList[i].position.x, 0, unitList[i].position.z);
								cache.DCorLC.directRhoPiAfterTheta(p);
								setBitToZero(stateB, p);
							}
							else
				n0++;
						}
						else{
							// deactivate its bit after lambda
							BitPosition p(unitList[i].position.x, 0, unitList[i].position.z);
							cache.DCorLC.directRhoPiAfterTheta(p);
							setBitToZero(stateB, p);
						}
		}
				}
			}
		}

		// compute the new cost
		int newWeight = alpha*cache.DCorLC.getMinReverseWeight(stateA) + beta*cache.DCorLC.getWeight(stateB) - max(alpha, beta) * 10;

		if (newWeight -2*max(alpha,beta)*n0 > maxCost){
			newColumn.index = 15;
			return false;
		}
	}

	return true;
}
