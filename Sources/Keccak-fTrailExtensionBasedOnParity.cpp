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

#include <fstream>
#include <iostream>
#include <sstream>
#include "Keccak-fDisplay.h"
#include "Keccak-fParity.h"
#include "Keccak-fParityBounds.h"
#include "Keccak-fTrails.h"
#include "Keccak-fTrailExtensionBasedOnParity.h"

// offset and basis of the affine space
AffineSpaceOfStates buildBasisBeforeChiGivenPatternAfterChi(const vector<SliceValue>& stateAfterChi)
{
	vector<SliceValue> offsetBeforeChi(stateAfterChi.size(), 0);
	vector<vector<SliceValue> > generatorsBeforeChi;   // to store the generators
	vector<vector<RowValue> > generatorsParities;

	for (unsigned int z = 0; z < stateAfterChi.size(); z++){
		if (stateAfterChi[z] != 0){
			for (int y = 0; y < nrRowsAndColumns; y++){
				RowValue row = getRowFromSlice(stateAfterChi[z], y);
				if (row != 0){
					if (row == 0x01 || row == 0x02 || row == 0x04 || row == 0x08 || row == 0x10){
						setRow(offsetBeforeChi, row, RowPosition(y, z));
					}
					for (int x = 0; x < nrRowsAndColumns; x++){
						RowValue value = 1 << x;
						if (value != row){
							vector<SliceValue> generator(stateAfterChi.size(), 0);
							setRow(generator, value, RowPosition(y, z));
							generatorsBeforeChi.push_back(generator);

							vector<RowValue> parity(stateAfterChi.size(), 0);
							getParity(generator, parity);
							generatorsParities.push_back(parity);
						}
					}
				}
			}
		}
	}

	vector<RowValue> OffsetParity(offsetBeforeChi.size(),0);
	getParity(offsetBeforeChi, OffsetParity);

	AffineSpaceOfStates affineSpace(stateAfterChi.size(), generatorsBeforeChi, generatorsParities, offsetBeforeChi, OffsetParity);
	return affineSpace;
}

AffineSpaceOfStates getdBasisAfterThetaGivenPatternBeforeChi(const KeccakFTrailExtension& keccakFTE, const AffineSpaceOfStates& basisBeforeChi)
{
	// build the image of the above space through inverse of pi and rho

	// offset
	vector<SliceValue> offsetAfterTheta(keccakFTE.laneSize, 0);
	keccakFTE.reverseLambdaAfterTheta(basisBeforeChi.offset, offsetAfterTheta);
	vector<RowValue> offsetAfterThetaParity(keccakFTE.laneSize, 0);
	getParity(offsetAfterTheta, offsetAfterThetaParity);

	// basis
	vector<vector<SliceValue> > generatorsAfterTheta;
	vector<vector<RowValue> > generatorsAfterThetaParities;
	for (unsigned int i = 0; i < basisBeforeChi.originalGenerators.size(); i++){
		vector<SliceValue> generator(keccakFTE.laneSize, 0);
		keccakFTE.reverseLambdaAfterTheta(basisBeforeChi.originalGenerators[i], generator);
		generatorsAfterTheta.push_back(generator);
		vector<RowValue> parity(keccakFTE.laneSize, 0);
		getParity(generator, parity);
		generatorsAfterThetaParities.push_back(parity);
	}

	AffineSpaceOfStates affineSpace(keccakFTE.laneSize, generatorsAfterTheta, generatorsAfterThetaParities, offsetAfterTheta, offsetAfterThetaParity);
	return affineSpace;

}


AffineSpaceOfStates buildBasisIntersectionWithKernel(const AffineSpaceOfStates& basis)
{
	vector<SliceValue> offsetIntersection = basis.offset;
	vector<RowValue> offsetParity = basis.offsetParity;
	unsigned int lanesize = offsetIntersection.size();

	vector<vector<SliceValue> > generatorsIntersection;
	vector<vector<RowValue> > generatorsParities;

	vector<ColumnPosition> columns;

	for (unsigned int i = 0; i < basis.originalGenerators.size(); i++){
		vector<SliceValue> generator = basis.originalGenerators[i];
		// 2 possibilities for the active bit:
		// 1) it is in a column that has odd active bits in the offset
		//    => if it is the first, add it to the offset and generate orbitals
		//		 otherwise do nothing, because it will be generated thanks to previous orbitals
		// 2) it is in a column that has even active bits in the offset => generator with orbital => look for other active bits
		for (unsigned int z = 0; z < lanesize; z++){
			if (generator[z] != 0){
				// look for the active column
				unsigned int x = 0;
				for (x = 0; x < nrRowsAndColumns; x++){
					if (getColumn(generator, x, z) != 0)
						break;
				}
				// check if it has already been considered
				bool checkedColumn = false;
				for (unsigned int k = 0; k < columns.size(); k++){
					if (z == columns[k].z && x == columns[k].x){
						checkedColumn = true;
						break;
					}
				}
				if (checkedColumn)
					break;
				else{
					columns.push_back(ColumnPosition(x, z));
					// offset column with an odd number of active bits
					if ((getHammingWeightColumn(getColumn(offsetIntersection, x, z)) % 2) != 0){
						offsetIntersection[z] ^= generator[z];
						offsetParity[z] ^= basis.originalParities[i][z];
					}
					// in any case look for other active bits in same column
					if (i != basis.originalGenerators.size() - 1){
						for (unsigned int j = i + 1; j < basis.originalGenerators.size(); j++){
							if (basis.originalParities[i][z] == basis.originalParities[j][z]){ // generators with active bits in same column
								vector<SliceValue> generatorWithOrbital = generator;
								generatorWithOrbital[z] ^= basis.originalGenerators[j][z];
								generatorsIntersection.push_back(generatorWithOrbital);
								vector<RowValue> parity(lanesize, 0);
								generatorsParities.push_back(parity);
							}
						}
					}
					break;
			}
		}
		}
	}
	AffineSpaceOfStates affineSpace(lanesize, generatorsIntersection, generatorsParities, offsetIntersection, offsetParity);
	return affineSpace;
}

bool intersectionWithKernel(const AffineSpaceOfStates& basis)
{
	unsigned int lanesize = basis.offset.size();

	// check if columns of the offset with an odd number of active bits have corresponding active bits in the basis
	for (unsigned int z = 0; z < lanesize; z++){
		if (basis.offset[z] != 0){
			for (int x = 0; x < nrRowsAndColumns; x++){
				if ((getHammingWeightColumn(getColumn(basis.offset, x, z)) % 2) != 0){
					bool existBasisBit = false;
					for (unsigned int i = 0; i < basis.originalGenerators.size(); i++){
						if (getHammingWeightColumn(getColumn(basis.originalGenerators[i], x, z)) != 0){
							existBasisBit = true;
							break;
						}
					}
					if (existBasisBit == false)
						return false;
				}
			}
		}
	}
	return true;
}


unsigned int index(int x)
{
	x %= 5;
	if (x<0) x += 5;
	return x;
}

void laneRotation(LaneValue& L, int offset, unsigned int aLaneSize)
{
	LaneValue mask = (LaneValue(~0)) >> (64 - aLaneSize);
	offset %= (int)aLaneSize;
	if (offset < 0) offset += aLaneSize;

	if (offset != 0) {
		L &= mask;
		L = (((LaneValue)L) << offset) ^ (((LaneValue)L) >> (aLaneSize - offset));
	}
	L &= mask;
}

void inverseThetaOnParity(vector<LaneValue>& A, unsigned int& aLaneSize)
{
	// parity
	vector<LaneValue> C = A;
	// coeff of Q
	const LaneValue inversePositions64[5] = {
		0xDE26BC4D789AF134ULL,
		0x09AF135E26BC4D78ULL,
		0xEBC4D789AF135E26ULL,
		0x7135E26BC4D789AFULL,
		0xCD789AF135E26BC4ULL };
	// Q for smaller laneSize, reduce mod 1+z^laneSize
	vector<LaneValue> inversePositions(5, 0);
	for (unsigned int z = 0; z<64; z += aLaneSize)
		for (unsigned int x = 0; x<5; x++)
			inversePositions[x] ^= inversePositions64[x] >> z;
	// compute inverse effect and apply
	vector<LaneValue> D(5,0);
	for (unsigned int z = 0; z<aLaneSize; z++) {
		for (unsigned int xOff = 0; xOff<5; xOff++)
			for (int x = 0; x<5; x++)
					if ((inversePositions[xOff] & 1) != 0) //coeff in Q is not zero so add parity
						D[index(x)] ^= C[index(x - xOff)];
		for (unsigned int xOff = 0; xOff<5; xOff++) {
			laneRotation(C[xOff], 1, aLaneSize);
			inversePositions[xOff] >>= 1;
		}
	}
	for (int x = 0; x<5; x++)
		A[index(x)] ^= D[x];
}

void inverseThetaOnParity(vector<RowValue>& A, unsigned int& aLaneSize)
{
	// parity
	vector<LaneValue> C(5, 0);
	fromSlicesToSheetsParity(A, C);
	inverseThetaOnParity(C, aLaneSize);
	fromSheetsToSlicesParity(C, A);
}

unsigned int getStartingSlice(AffineSpaceOfStates& basisAfterTheta){

	unsigned int start = 0;
	unsigned int maxNrEmptySlices = 0;

	unsigned int laneSize = basisAfterTheta.offset.size();
	vector<RowValue> parityOffset = basisAfterTheta.offsetParity;
	vector< vector<RowValue> > parityGenerators = basisAfterTheta.offsetParities;

	// slices containing a x
	vector<bool> generatorSlices(laneSize, false);
	for (unsigned int i = 0; i < parityGenerators.size(); i++){
		for (unsigned int j = 0; j < laneSize; j++){
			if (parityGenerators[i][j] != 0)
				generatorSlices[j] = true;
		}
	}

	// slices containing 1 or x
	vector<bool> activeSlices = generatorSlices;
	for (unsigned int i = 0; i < laneSize; i++){
		if (parityOffset[i] != 0)
			activeSlices[i] = true;
	}

	for (unsigned int i = 0; i < laneSize; i++){
		if (activeSlices[i] == true){
			// look for the number of empty slices before i
			unsigned int j = 1;
			for (j = 1; j < laneSize; j++){
				if (activeSlices[(i - j) % laneSize] == true)
					break;
			}
			unsigned int nrEmptySlices = j - 1;
			if (nrEmptySlices > maxNrEmptySlices){
				maxNrEmptySlices = nrEmptySlices;
				start = i;
			}
		}
	}
	return start;
}

RowValue rotateRow(RowValue r, int offset)
{
	offset %= 5;
	if (offset < 0) offset += 5;

	RowValue a;
	if (offset != 0) {
		r &= 0x1F;
		a = (r << offset) ^ (r >> (5 - offset));
		a &= 0x1F;
	}
	else a = r & 0x1F;
	return a;
}

vector<vector<unsigned int> > getnNrBasisVectorsPerColumn(AffineSpaceOfStates& basis)
{
	unsigned int laneSize = basis.offset.size();
	vector<vector<unsigned int> > basisVectors(laneSize,vector<unsigned int>(5,0));
	for (unsigned int i = 0; i < basis.originalGenerators.size(); i++){
		vector<SliceValue> vector = basis.originalGenerators[i];
		for (unsigned int z = 0; z < laneSize; z++){
			if (vector[z] != 0){
				for (unsigned int x = 0; x < nrRowsAndColumns; x++){
					basisVectors[z][x]+=getHammingWeightColumn(getColumn(vector, x, z));
				}
				break; // this is an optimization when a vector of the basis has a single active bit
			}
		}
	}
	return basisVectors;
}

vector<vector<unsigned int> > getnNrBasisVectorsPerColumnNoOrbitals(AffineSpaceOfStates& basis)
{
	unsigned int laneSize = basis.offset.size();
	vector<vector<unsigned int> > basisVectors(laneSize, vector<unsigned int>(5, 0));
	for (unsigned int i = 0; i < basis.offsetGenerators.size(); i++){
		vector<SliceValue> vector = basis.offsetGenerators[i];
		for (unsigned int z = 0; z < laneSize; z++){
			if (vector[z] != 0){
				for (unsigned int x = 0; x < nrRowsAndColumns; x++){
					basisVectors[z][x] += getHammingWeightColumn(getColumn(vector, x, z));
				}
			}
		}
	}
	return basisVectors;
}

vector< vector<RowValue> > getRowValuesFromBasis(AffineSpaceOfStates& basis)
{
	unsigned int laneSize = basis.offset.size();
	vector< vector<RowValue> > values;
	for (unsigned int z = 0; z < laneSize; z++){
		RowValue offsetForRow = basis.offsetParity[z];
		vector<RowValue> basisForRow;
		for (unsigned int i = 0; i < basis.offsetParities.size(); i++){
			if (basis.offsetParities[i][z] != 0)
				basisForRow.push_back(basis.offsetParities[i][z]);
		}
		values.push_back(getSetOfRowValues(offsetForRow, basisForRow));
	}
	return values;
}



parityBackwardIterator::parityBackwardIterator(const KeccakFPropagation& aDCorLC, vector<SliceValue> &aOffset, vector<RowValue> &aOffsetParity, vector<vector<RowValue> > &aValuesX, vector<vector<unsigned int> > &aNrBasisVectors, unsigned int& aStart, RowValue& aGuess, unsigned int aMaxWeight)
: parityIterator(aDCorLC), start(aStart), offsetParity(aOffsetParity), offset(aOffset), rowsValues(aValuesX), nrBasisVectors(aNrBasisVectors), current(aStart), maxWeight(aMaxWeight), end(false), initialized(false), empty(true) {
	b = offsetParity;
	for (unsigned int i = 0; i < laneSize; i++)
		indexes.push_back(-1);
	a = vector<RowValue>(laneSize, 0);
	a[start] = aGuess;
	lowerBound.push_back(0);
}

void parityBackwardIterator::initialize()
{
	index = 0;
	if (first()) {
		getParity();
		end = false;
		empty = false;
	}
	else {
		end = true;
		empty = true;
	}
	initialized = true;
}

bool parityBackwardIterator::first()
{
	bool firstIteration = true;
	do{
		RowValue workingRow; // Dummy row. iterateRowAndPush will initialize it according to offset and basis
		if(!firstIteration)
			popRow(workingRow);
		do{
			if (iterateRowAndPush(workingRow))
				current=(current - 1 + laneSize) % laneSize;
			else{
				current = (current + 1 ) % laneSize;
				if (!popRow(workingRow))
					return false;
			}
		} while (lowerBound.size() <= laneSize);
		current = (start + 1) % laneSize;
		firstIteration = false;
		if (checkCompatibility())
			return true;
	} while (true);
}


bool parityBackwardIterator::next()
{
	RowValue workingRow;
	do{
		popRow(workingRow);
		do{
			if (iterateRowAndPush(workingRow))
				current = (current - 1 + laneSize) % laneSize;
			else{
				current = (current + 1) % laneSize;
				if (!popRow(workingRow))
					return false;
			}
		} while (lowerBound.size() <= laneSize);
		current = (start + 1) % laneSize;
		if (checkCompatibility())
			return true;
	} while (true);
}



bool parityBackwardIterator::iterateRowAndPush(RowValue& workingRow)
{
	do{
		if (!successorOf(workingRow))
			return false;

		pushRow(workingRow);
		if (lowerBound.back() <= (unsigned)maxWeight)
			return true;
		popRow(workingRow);

	} while (true);
}

bool parityBackwardIterator::successorOf(RowValue& workingRow){

	unsigned int ii = indexes[current]+1;
	if (ii < rowsValues[current].size()){
		workingRow = rowsValues[current][ii];
		indexes[current] = ii;
		return true;
	}
	else
		return false;

}

void parityBackwardIterator::pushRow(RowValue& workingRow){

	b[current] = workingRow;

	for (unsigned int i = 1; i <= (current - start - 1 + laneSize) % laneSize; i++)
		indexes[(current - i + laneSize) % laneSize] = -1;
	lowerBound.push_back(lowerBound.back() + 2 * getNrActiveRowsBeforeTheta(current));

	if (current != ((start + 1) % laneSize))
		a[(current - 1 + laneSize) % laneSize] = rotateRow(b[current] ^ a[current] ^ rotateRow(a[current], 1), 1);

}

unsigned int parityBackwardIterator::getNrActiveRowsBeforeTheta(unsigned int& z){

	unsigned int maxNrRows = 0;

	for (unsigned int x = 0; x < nrRowsAndColumns; x++){
		unsigned int activeBasis = 0;
		unsigned int nrRows = 0;
		bool oddBeforeTheta = getBit(a, x, z) == 1 ? true : false;
		bool oddAfterTheta = getBit(b, x, z) == 1 ? true : false;
		unsigned int nrOffsetBits = getHammingWeightColumn(getColumn(offset,x,z));
		if (oddBeforeTheta){
			if (oddAfterTheta){
				// UOC: at least 1 active row; worst case when orbitals are passive;
				// if there is an odd number of bits of the offset, then it is that number
				// otherwise it means there is at least one vector of the basis active and thus offset+1 active bits
				activeBasis = (nrOffsetBits % 2) == 0 ? 1 : 0;
				nrRows = nrOffsetBits + activeBasis;
			}
			else{
				// AOC: at least 1 active row; worst case when all orbitals are active
				// if there is an even number of bits in the offset then 5-offset-2*maxorbitals
				// otherwise at least one bit of the basis is active, thus 5-(offset+1)-2*maxorbitals
				activeBasis = (nrOffsetBits % 2) == 0 ? 0 : 1;
				unsigned int nrOrbitals = getNrOrbitals(x, z, activeBasis);
				nrRows = max(5 - (nrOffsetBits + activeBasis) - 2 * nrOrbitals, (unsigned int)(1));
			}
		}
		else{
			if (oddAfterTheta){
				// AEC: min active rows is 0; worst case when all orbitals are active at b
				// if odd number bits in the offset, then 5-offset-2*maxorbitals
				// otherwise, at least 1 basis vector is active and thus 5-(offset+1)-2*orbitals
				activeBasis = (nrOffsetBits % 2) == 0 ? 1 : 0;
				unsigned int nrOrbitals = getNrOrbitals(x, z, activeBasis);
				nrRows = max(5 - (nrOffsetBits + activeBasis) - 2 * nrOrbitals, (unsigned int)(0));
			}
			else{
				// UEC: the min number is zero; worst case when orbitals are passive
				// if even number of bits in the offset then it is that number
				// otherwise it means there is at least one vector of the basis active and thus offset+1 at a
				activeBasis = (nrOffsetBits % 2) == 0 ? 0 : 1;
				nrRows = nrOffsetBits + activeBasis;
			}
		}
		if (nrRows > maxNrRows)
			maxNrRows = nrRows;
		if (maxNrRows == 5)
			break;
	}
	return maxNrRows;
}

unsigned int parityBackwardIterator:: getNrOrbitals(unsigned int& x, unsigned int& z, unsigned int& activeBasis){

	return (nrBasisVectors[z][x] - activeBasis) / 2;

}

bool parityBackwardIterator::popRow(RowValue& workingRow){

	if (lowerBound.size() == 1)
		return false;
	workingRow = b[current];
	b[current] = offsetParity[current];
	lowerBound.pop_back();
	return true;
}

bool parityBackwardIterator::checkCompatibility(){

	RowValue temp = rotateRow(b[(start + 1) % laneSize] ^ a[(start + 1) % laneSize] ^ rotateRow(a[(start + 1) % laneSize], 1), 1);
	if (temp == a[start])
		return true;
	else
		return false;
}

void parityBackwardIterator::getParity()
{
	parity.clear();
	parity = b;
}

bool parityBackwardIterator::isEnd()
{
	if (!initialized) initialize();
	return end;
}
bool parityBackwardIterator::isEmpty()
{
	if (!initialized) initialize();
	return empty;
}

bool parityBackwardIterator::isBounded()
{
	if (!initialized) initialize();
	return false;
}

UINT64 parityBackwardIterator::getIndex()
{
	if (!initialized) initialize();
	return index;
}

UINT64 parityBackwardIterator::getCount()
{
	if (!initialized) initialize();
	return (end ? index : 0);
}

void parityBackwardIterator::operator++()
{

	if (!initialized){
		initialize();
	}
	else{
		if (!end) {
			++index;
			if (next())
				getParity();
			else
				end = true;
		}
	}
}

const vector<RowValue>& parityBackwardIterator::operator*()
{
	if (!initialized) initialize();
	return parity;
}

vector<RowValue> getSetOfRowValues(RowValue& offset, vector<RowValue>& basis){

	vector< RowValue > all;

	for (unsigned int n = 0; n < pow(2, basis.size()); n++){
		RowValue value = offset;
		for (unsigned int j = 0; j < basis.size(); j++){
			if (((n >> j) & 1) != 0)
				value ^= basis[j];
		}
		all.push_back(value);
	}
	return all;

}

void setBasisPerInput(const KeccakFTrailExtension& keccakFTE, vector<AffineSpaceOfRows>& aBasisPerInput)
{
	aBasisPerInput = keccakFTE.affinePerInput;

	// change the basis for 0x1F
	AffineSpaceOfRows a;
	a.setOffset(RowValue(0x04));
	a.addGenerator(RowValue(0x06));
	a.addGenerator(RowValue(0x05));
	a.addGenerator(RowValue(0x14));
	a.addGenerator(RowValue(0x0C));
	aBasisPerInput[0x1F] = a;

	// change the offsets for all values
	for (unsigned int i = 1; i < (1 << nrRowsAndColumns) - 1; i++){
		RowValue row = i & 0x1F;
		if (row == 0x4 || row == 0x6 || row == 0x5 || row == 0xD || row == 0x7 || row == 0x17)
			aBasisPerInput[i].setOffset(RowValue(0x04));
		if (row == 0x2 || row == 0x3 || row == 0x12 || row == 0x16 || row == 0x13 || row == 0x1B)
			aBasisPerInput[i].setOffset(RowValue(0x02));
		if (row == 0x1 || row == 0x11 || row == 0x9 || row == 0xB || row == 0x19 || row == 0x1D)
			aBasisPerInput[i].setOffset(RowValue(0x01));
		if (row == 0x10 || row == 0x18 || row == 0x14 || row == 0x15 || row == 0x1C || row == 0x1E)
			aBasisPerInput[i].setOffset(RowValue(0x10));
		if (row == 0x8 || row == 0xC || row == 0x0A || row == 0x1A || row == 0xE || row == 0xF)
			aBasisPerInput[i].setOffset(RowValue(0x08));
	}
}


AffineSpaceOfStates buildBasisAfterChiGivenPatternBeforeChi(vector<AffineSpaceOfRows>& affinePerInput, const vector<SliceValue>& stateBeforeChi)
{
	unsigned int laneSize = stateBeforeChi.size();

	vector<SliceValue> offset(laneSize, 0);
	vector<vector<SliceValue> > generators;   // to store the generators
	vector<vector<RowValue> > parities;

	for (unsigned int z = 0; z < laneSize; z++){
		if (stateBeforeChi[z] != 0){
			for (int y = 0; y < nrRowsAndColumns; y++){
				RowValue row = getRowFromSlice(stateBeforeChi[z], y);
				if (row != 0){
					offset[z] ^= getSliceFromRow(affinePerInput[row].offset,y);
					for (unsigned int i = 0; i < affinePerInput[row].generators.size(); i++){
						vector<SliceValue> generator(laneSize, 0);
						setRow(generator, affinePerInput[row].generators[i], y, z);
						generators.push_back(generator);

						vector<RowValue> genParity(laneSize, 0);
						getParity(generator, genParity);
						parities.push_back(genParity);
					}
				}
			}
		}
	}

	vector<RowValue> OffsetParity(laneSize, 0);
	getParity(offset, OffsetParity);

	AffineSpaceOfStates affineSpace(laneSize, generators, parities, offset, OffsetParity);
	return affineSpace;
}

stateForwardIterator::stateForwardIterator(const KeccakFPropagation& aDCorLC, AffineSpaceOfStates& aBasis, unsigned int aBudgetWeight) :
DCorLC(aDCorLC), laneSize(aDCorLC.laneSize), offset(aBasis.offset), offsetParity(aBasis.offsetParity), budgetWeight(aBudgetWeight), end(false), initialized(false), empty(true)
{
	for (unsigned int i = 0; i < aBasis.offsetGenerators.size(); i++){
		basis.push_back(aBasis.offsetGenerators[i]);
		basisParity.push_back(aBasis.offsetParities[i]);
	}
	firstOrbital = aBasis.offsetGenerators.size();
	for (unsigned int i = 0; i < aBasis.kernelGenerators.size(); i++){
		basis.push_back(aBasis.kernelGenerators[i]);
		vector<RowValue> allZero(laneSize, 0);
		basisParity.push_back(allZero);
	}

	for (unsigned int i = 0; i < basis.size(); i++){
		for (unsigned int z = 0; z < laneSize; z++){
			if (basis[i][z] != 0){
				slices.push_back(z);
				break;
			}
		}
	}
	numPerColumn = getnNrBasisVectorsPerColumnNoOrbitals(aBasis);

}

void stateForwardIterator::initialize()
{
	indexes.push_back(-1);
	stateAtA.push_back(offset);
	vector<SliceValue> B(laneSize,0);
	DCorLC.directLambda(stateAtA.back(), B);
	stateAtB.push_back(B);
	weight.push_back(DCorLC.getWeight(B));
	vector<RowValue> parity(laneSize, 0);
	getParity(stateAtA.back(), parity);
	C.push_back(parity);
	vector<RowValue> effect(laneSize, 0);
	DCorLC.directThetaEffectFromParities(parity, effect);
	D.push_back(effect);
	index = 0;
	if (first()) {
		end = false;
		empty = false;
	}
	else {
		end = true;
		empty = true;
	}
	initialized = true;

}

bool stateForwardIterator::first()
{
	if (weight.back() <= budgetWeight)
		return true;  // return the offset
	else
		return next(); // try to add something
}

bool stateForwardIterator::next()
{
	do{
		int workingIndex = indexes.back();
		do{
			if (iterateIndexAndPush(workingIndex)) // it can be added if weight<budget or if weight>budget but there exists a super-entangled vector or an orbital to be added later
			{
				break;
			}
			else{
				if (!popIndex(workingIndex))
					return false;
			}
		} while (true);
	} while (weight.back() > budgetWeight); // go on adding super-entangled vector or orbitals

	return true; // return only when weight<budget

}

bool stateForwardIterator::iterateIndexAndPush(int& workingIndex)
{
	do{
		if (!successorOf(workingIndex))
			return false;
		pushIndex(workingIndex);
		if (canAfford(workingIndex))
			return true;
		popIndex(workingIndex);
		//
	} while (true);
}

bool stateForwardIterator::successorOf(int& workingIndex)
{
	if (workingIndex == (int)basis.size()-1)
		return false;
	workingIndex += 1;
	return true;
}

bool stateForwardIterator::canAfford(int& workingIndex)
{
	// here check the weight, but also super-entanglement
	if (weight.back() <= budgetWeight)
		return true;

	unsigned int numOrbital = 0;
	// unsigned int numVec = 0;
	unsigned int totalContribution = 0;
	if (workingIndex >= firstOrbital){
		numOrbital = basis.size() - workingIndex - 1;
		totalContribution = 4 * numOrbital;
	}
	else{
		numOrbital = basis.size() - firstOrbital;
		totalContribution = 4 * numOrbital;
		for (unsigned int i = workingIndex + 1; i<(unsigned int)firstOrbital; i++){
			unsigned int z = slices[i];
			unsigned int k = 0;
			unsigned int contribution = 0;
			for (unsigned int x = 0; x < nrRowsAndColumns; x++){

				if (getHammingWeightColumn(getColumn(basis[i], x, z)) != 0){
					k = numPerColumn[z][x]-1;
					if (k == 0){ // parity can not turn from odd to even
						contribution = 2 + 2 * numPerColumn[z][(x + 1) % 5] + 2 * numPerColumn[(z + 1) % laneSize][(x -1) % 5];
					}
					else
						contribution = 22;
					totalContribution += contribution;
				}

			}

		}
	}
	if (weight.back() > (budgetWeight + totalContribution))
		return false;
	if (thereAreSuperEntangledIndexes()) // weight>budget but there exists a super - entangled vector to be added later
		return true;
	return false; // weight exceeds the budget and remaining vectors are not entangled (will only add new bits)
}

void stateForwardIterator::pushIndex(int& workingIndex)
{
	indexes.push_back(workingIndex);
	// push in A
	stateAtA.push_back(stateAtA.back());
	unsigned int z = slices[workingIndex];
	stateAtA.back()[z] ^= basis[workingIndex][z];

	stateAtB.push_back(stateAtB.back()); // this can be optimized
	DCorLC.directLambda(stateAtA.back(),stateAtB.back());
	weight.push_back(DCorLC.getWeight(stateAtB.back())); // this can be optimized

	vector<RowValue> parity(laneSize,0);
	getParity(stateAtA.back(), parity);
	C.push_back(parity);
	vector<RowValue> effect(laneSize, 0);
	DCorLC.directThetaEffectFromParities(C.back(), effect);
	D.push_back(effect);

}

bool stateForwardIterator::popIndex(int& workingIndex)
{
	if (indexes.size() == 1)
		return false;

	workingIndex = indexes.back();

	indexes.pop_back();
	stateAtA.pop_back();
	stateAtB.pop_back();
	weight.pop_back();
	C.pop_back();
	D.pop_back();

	return true;
}

bool stateForwardIterator::thereAreSuperEntangledIndexes()
{
	// for each remaining vector, check if it is entangled or run-modifying

	// orbitals do not change parity,
	// check if they are added to AC or cancel bits
	unsigned int k = max(indexes.back()+1, firstOrbital);
	for (unsigned int i = k; i < basis.size(); i++){
		if (isAddingBitToAC(i))
			return true;
		if (isBitOverlapping(i))
			return true;
	}

	// look if they are bit-overlapping, entangled or runModifying
	for (unsigned int i = indexes.back()+1; i < (unsigned int)firstOrbital; i++)
	{
		if (isBitOverlapping(i))
			return true;
		if (isEntangled(i))
			return true;
		if (isRunModifying(i))
			return true;
	}
	return false;

}


bool stateForwardIterator::isBitOverlapping(unsigned int& a)
{
	unsigned int activeSlice = slices[a];
	if ((basis[a][activeSlice] & stateAtA.back()[activeSlice]) != 0)
		return true;
	return false;
}

bool stateForwardIterator::isAddingBitToAC(unsigned int& a)
{
	unsigned int activeSlice = slices[a];
	// for each column
	for (unsigned int i = 0; i < nrRowsAndColumns; i++){
		// if there are active bits in this column
		if (getHammingWeightColumn(getColumn(basis[a], i, activeSlice)) != 0){
			// if they are in AC
			if ((getBit(D.back(), i, activeSlice)&1) != 0)
				return true;
		}
	}
	return false;
}

bool stateForwardIterator::isEntangled(unsigned int& a)
{
	// the vector adds bits to ACs
	if (isAddingBitToAC(a))
		return true;
	// the vector makes an active UC become AC
	vector<RowValue> effect(laneSize, 0);
	DCorLC.directThetaEffectFromParities(basisParity[a], effect);
	// look in the two slices for AC
	for (unsigned int j = 0; j <= 1; j++){
		unsigned int z = slices[a] + j;
		for (unsigned int i = 0; i < nrRowsAndColumns; i++){
			// look if it is affected
			if ( (getBit(effect, i, z) & 1) != 0) {
				// look if it was UC and active
				if ((getBit(D.back(), i, z) & 1) == 0 && (getHammingWeightColumn(getColumn(stateAtA.back(), i, z)))>0)
					return true;
			}
		}
	}
	return false;
}

bool stateForwardIterator::isRunModifying(unsigned int& a)
{
	unsigned int z = slices[a];
	// extends a run or connects two runs, i.e. a column affected by a is already affected
	vector<RowValue> effect(laneSize, 0);
	DCorLC.directThetaEffectFromParities(basisParity[a], effect);
	// look in the two slices for AC
	for (unsigned int j = 0; j <= 1; j++){
		for (unsigned int i = 0; i < nrRowsAndColumns; i++){
			// look if it is affected and was affected
			if (((getBit(effect, i, z+j) & 1) != 0) && ((getBit(D.back(), i, z+j) & 1) != 0)) {
				return true;
			}
		}
	}

	// cancels or splits a run, i.e. makes an odd column even
	for (unsigned int i = 0; i < nrRowsAndColumns; i++){
		// look if it is odd and was odd
		if (((getBit(basisParity[a], i, z) & 1) != 0) && ((getBit(C.back(), i, z) & 1) != 0)) {
			return true;
		}
	}
	return false;
}


bool stateForwardIterator::isEnd()
{
	if (!initialized) initialize();
	return end;
}
bool stateForwardIterator::isEmpty()
{
	if (!initialized) initialize();
	return empty;
}

bool stateForwardIterator::isBounded()
{
	if (!initialized) initialize();
	return false;
}

UINT64 stateForwardIterator::getIndex()
{
	if (!initialized) initialize();
	return index;
}

UINT64 stateForwardIterator::getCount()
{
	if (!initialized) initialize();
	return (end ? index : 0);
}

void stateForwardIterator::operator++()
{

	if (!initialized){
		initialize();
	}
	else{
		if (!end) {
			++index;
			if (next())
			{
			}
			else
				end = true;
		}
	}

}

const vector<SliceValue>& stateForwardIterator::operator*()
{
	if (!initialized) initialize();
	return stateAtA.back();
}


void KeccakFTrailExtensionBasedOnParity::forwardExtendTrailsInTheKernel(TrailIterator& trailsIn, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight)
{
	progress.stack("File", trailsIn.getCount());
	for (; !trailsIn.isEnd(); ++trailsIn) {
		forwardExtendTrailInTheKernel(*trailsIn, trailsOut, nrRounds, maxTotalWeight);
		++progress;
	}
	progress.unstack();
}

void KeccakFTrailExtensionBasedOnParity::forwardExtendTrailInTheKernel(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight)
{
	if (trail.stateAfterLastChiSpecified)
		throw KeccakException("KeccakFTrailExtension::forwardExtendTrail() can work only with trail cores or trail prefixes.");
	recurseForwardExtendTrailInTheKernel(trail, trailsOut, nrRounds, maxTotalWeight);
}

void KeccakFTrailExtensionBasedOnParity::recurseForwardExtendTrailInTheKernel(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight)
{
	int baseWeight = trail.totalWeight;
	int baseNrRounds = trail.getNumberOfRounds();
	int curNrRounds = baseNrRounds + 1;
	int curWeight = trail.weights.back();
	int maxWeightOut = maxTotalWeight - baseWeight
		- knownBounds.getMinWeight(nrRounds - baseNrRounds - 1);
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
		for (vector<vector<SliceValue> >::const_iterator i = compatibleStates.begin(); i != compatibleStates.end(); ++i) {
			int weightOut = getWeight(*i);
			int curWeight = baseWeight + weightOut;
			if (curNrRounds == (int)nrRounds) {
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
		// build the basis
		AffineSpaceOfStates base = buildStateBase(trail.states.back());
		// generate the iterator
		SlicesAffineSpaceIterator i = base.getIteratorInKernel();
		progress.stack(synopsis + " [affine base]", i.getCount());
		for (; !i.isEnd(); ++i) {
			int weightOut = getWeight(*i);
			int curWeight = baseWeight + weightOut;
			if (curNrRounds == (int)nrRounds) {
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


void KeccakFTrailExtensionBasedOnParity::backwardExtendTrailsInTheKernel(TrailIterator& trailsIn, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight)
{
	progress.stack("File", trailsIn.getCount());
	for (; !trailsIn.isEnd(); ++trailsIn) {
		backwardExtendTrailInTheKernel(*trailsIn, trailsOut, nrRounds, maxTotalWeight);
		++progress;
	}
	progress.unstack();
}

void KeccakFTrailExtensionBasedOnParity::backwardExtendTrailInTheKernel(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight)
{
	bool isPrefix = trail.firstStateSpecified;
	if (isPrefix) {
		recurseBackwardExtendTrailInTheKernel(trail, trailsOut, nrRounds, maxTotalWeight);
	}
	else {
		Trail trimmedTrailPrefix;
		for (unsigned int i = 1; i<trail.states.size(); i++)
			trimmedTrailPrefix.append(trail.states[i], trail.weights[i]);
		recurseBackwardExtendTrailInTheKernel(trimmedTrailPrefix, trailsOut, nrRounds, maxTotalWeight);
	}
}

void KeccakFTrailExtensionBasedOnParity::recurseBackwardExtendTrailInTheKernel(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight)
{
	if (!allPrefixes && (nrRounds == (trail.getNumberOfRounds() + 1))) {
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
	// int baseNrRounds = trail.getNumberOfRounds();
	int maxWeightOut = maxTotalWeight - baseWeight
		- knownBounds.getMinWeight(1);
	if (maxWeightOut < knownBounds.getMinWeight(1))
		return;

	// compute state after chi
	vector<SliceValue> stateAfterChi(laneSize, 0);
	reverseLambda(trail.states[0], stateAfterChi);
	// build basis for the envelope of chi^-1(stateAfterChi)
	AffineSpaceOfStates basisBeforeChi = buildBasisBeforeChiGivenPatternAfterChi(stateAfterChi);
	// apply inverse of pi and rho to obtain basis after theta
	AffineSpaceOfStates basisAfterTheta = getdBasisAfterThetaGivenPatternBeforeChi(*this, basisBeforeChi);
	// check min weight
	int curWeight = baseWeight + getWeight(basisBeforeChi.offset) + getMinReverseWeight(basisAfterTheta.offset);
	if (curWeight > maxTotalWeight)
		return;
	// check whether intersection with the kernel is empty
	if (!intersectionWithKernel(basisAfterTheta))
		return;
	// build basis for the intersection with the kernel
	AffineSpaceOfStates basisIntersection = buildBasisIntersectionWithKernel(basisAfterTheta);
	// build the iterator to all states in the affine space with parity in the kernel
	SlicesAffineSpaceIterator iterator = basisIntersection.getIteratorInKernel();
	//SlicesAffineSpaceIterator iterator = basisAfterTheta.getIteratorInKernel();

	for (; !iterator.isEnd(); ++iterator){
			vector<SliceValue> state = *iterator; // a1 after theta
			vector<SliceValue> stateBeforeChi; // b1 after lambda
		directLambdaAfterTheta(state, stateBeforeChi);
			if (isChiCompatible(stateBeforeChi, stateAfterChi)){ // check whether b1 is compatible with a2 through chi
			Trail newTrail;
			newTrail.setFirstStateReverseMinimumWeight(getMinReverseWeight(state));
			newTrail.append(stateBeforeChi, getWeight(stateBeforeChi));
			for (unsigned int i = 0; i<trail.states.size(); i++)
				newTrail.append(trail.states[i], trail.weights[i]);
			if ((int)newTrail.totalWeight <= maxTotalWeight){
				trailsOut.fetchTrail(newTrail);
			}
		}
	}
}
}

void KeccakFTrailExtensionBasedOnParity::backwardExtendTrailsOutsideKernel(TrailIterator& trailsIn, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight)
{
	progress.stack("File", trailsIn.getCount());
	for (; !trailsIn.isEnd(); ++trailsIn) {
		backwardExtendTrailOutsideKernel(*trailsIn, trailsOut, nrRounds, maxTotalWeight);
		++progress;
	}
	progress.unstack();
}

void KeccakFTrailExtensionBasedOnParity::backwardExtendTrailOutsideKernel(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight)
{
	bool isPrefix = trail.firstStateSpecified;
	if (isPrefix) {
		int maxRevWeight = maxTotalWeight - max(int(trail.totalWeight + 2), knownBounds.getMinWeight(trail.getNumberOfRounds()));
		recurseBackwardExtendTrailOutsideKernel(trail, trailsOut, nrRounds, maxTotalWeight, maxRevWeight, true);
	}
	else {
		Trail trimmedTrailPrefix;
		for (unsigned int i = 1; i<trail.states.size(); i++)
			trimmedTrailPrefix.append(trail.states[i], trail.weights[i]);
		int maxRevWeight = maxTotalWeight - trail.totalWeight;
		recurseBackwardExtendTrailOutsideKernel(trimmedTrailPrefix, trailsOut, nrRounds, maxTotalWeight, maxRevWeight, allPrefixes);
	}
}

void KeccakFTrailExtensionBasedOnParity::recurseBackwardExtendTrailOutsideKernel(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight, int maxRevWeight, bool allPrefixes)
{
	if (!allPrefixes && (nrRounds == (trail.getNumberOfRounds() + 1))) {
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
		int baseNrRounds = trail.getNumberOfRounds();
		int maxWeightOut = maxTotalWeight - baseWeight
			- knownBounds.getMinWeight(nrRounds - baseNrRounds - 1);
		if (maxWeightOut < knownBounds.getMinWeight(1))
			return;
		// compute state after chi by applying inverse of lambda on states[0]
		// this is a2
		vector<SliceValue> stateAfterChi(laneSize, 0);
		reverseLambda(trail.states[0], stateAfterChi);

		// build basis for the envelope of chi^-1(stateAfterChi)
		// this is a basis for b1
		AffineSpaceOfStates basisBeforeChi = buildBasisBeforeChiGivenPatternAfterChi(stateAfterChi);

		// apply inverse of pi and rho to obtain basis after theta
		// this is a basis of theta(a1)
		AffineSpaceOfStates basisAfterTheta = getdBasisAfterThetaGivenPatternBeforeChi(*this, basisBeforeChi);

		// here starts the iteration over parity patterns after theta.
		// for each row of the parity construct the set of possible values based on the basis vectors
		vector< vector<RowValue> > values = getRowValuesFromBasis(basisAfterTheta);
		// for each column get the number of basis vectors with active bits int that column
		vector<vector<unsigned int> > nrVectors = getnNrBasisVectorsPerColumn(basisAfterTheta);

		// the slice from which the iterator starts to construct the parity pattern
		//unsigned int start = 0;
		unsigned int start = getStartingSlice(basisAfterTheta);
		// guess the value of the starting row before theta
		for (unsigned int k = 0; k < 32; k++){
			RowValue guess = k;
			// construct the iterator
			parityBackwardIterator iterator(*this, basisAfterTheta.offset, basisAfterTheta.offsetParity, values, nrVectors, start, guess, maxRevWeight);
			for (; !iterator.isEnd(); ++iterator) {
				vector<RowValue> parity = *iterator;
				// exclude in kernel case
				vector<RowValue> allzero(laneSize, 0);
				if (parity != allzero){
					// if a pattern below the budget is found, then build states with given parity using the basis
					SlicesAffineSpaceIterator stateIterator = basisAfterTheta.getIteratorWithGivenParity(parity);
					for (; !stateIterator.isEnd(); ++stateIterator){
						vector<SliceValue> state = *stateIterator; //theta(a1)
						vector<SliceValue> stateBeforeChi; //b1
						directLambdaAfterTheta(state, stateBeforeChi);
						if (isChiCompatible(stateBeforeChi, stateAfterChi)){// (b1,a2)
							Trail newTrail;
							vector<SliceValue> stateBeforeLambda; //a1
							reverseLambda(stateBeforeChi, stateBeforeLambda);
							newTrail.setFirstStateReverseMinimumWeight(getMinReverseWeight(stateBeforeLambda));
							newTrail.append(stateBeforeChi, getWeight(stateBeforeChi)); // add b1
							for (unsigned int i = 0; i < trail.states.size(); i++)
								newTrail.append(trail.states[i], trail.weights[i]); //add b2...
							if ((int)newTrail.totalWeight <= maxTotalWeight)
								trailsOut.fetchTrail(newTrail);
						}
					}
				}
			}
		}
	}
}

void KeccakFTrailExtensionBasedOnParity::forwardExtendTrailsOutsideKernel(TrailIterator& trailsIn, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight)
{
	progress.stack("File", trailsIn.getCount());
	vector<AffineSpaceOfRows> basisPerInput;
	setBasisPerInput(*this, basisPerInput);
	for (; !trailsIn.isEnd(); ++trailsIn) {
		forwardExtendTrailOutsideKernel(*trailsIn, trailsOut, nrRounds, maxTotalWeight, basisPerInput);
		++progress;
	}
	progress.unstack();
}

void KeccakFTrailExtensionBasedOnParity::forwardExtendTrailOutsideKernel(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight, vector<AffineSpaceOfRows> basisPerInput)
{
	if (trail.stateAfterLastChiSpecified)
		throw KeccakException("KeccakFTrailExtension::forwardExtendTrail() can work only with trail cores or trail prefixes.");
	recurseForwardExtendTrailOutsideTheKernel(trail, trailsOut, nrRounds, maxTotalWeight, basisPerInput);
}

void KeccakFTrailExtensionBasedOnParity::recurseForwardExtendTrailOutsideTheKernel(const Trail& trail, TrailFetcher& trailsOut, unsigned int nrRounds, int maxTotalWeight, vector<AffineSpaceOfRows> basisPerInput)
{
	int baseWeight = trail.totalWeight;
	int baseNrRounds = trail.getNumberOfRounds();
	int curNrRounds = baseNrRounds + 1;
	int curWeight = trail.weights.back();
	int maxWeightOut = maxTotalWeight - baseWeight
		- knownBounds.getMinWeight(nrRounds - baseNrRounds - 1);
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
		for (vector<vector<SliceValue> >::const_iterator i = compatibleStates.begin(); i != compatibleStates.end(); ++i) {
			int weightOut = getWeight(*i);
			int curWeight = baseWeight + weightOut;
			if (curNrRounds == (int)nrRounds) {
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
		vector<unsigned int> activeSlices;
		// build basis for the state after chi given the state before chi
		// this is a basis for b2
		AffineSpaceOfStates base = buildBasisAfterChiGivenPatternBeforeChi(basisPerInput, trail.states.back());
		// generate the iterator on the vaues of state after chi
		stateForwardIterator iterator(*this, base, maxWeightOut);
		progress.stack(synopsis + " [affine base]", iterator.getCount());
		for (; !iterator.isEnd(); ++iterator) {
			// exclude in kernel cases
			vector<RowValue> parity;
			getParity(*iterator, parity);
			vector<RowValue> allzero(laneSize, 0);
			if (parity != allzero){
				// compute state after lambda and construct the extended trail
				vector<SliceValue> stateAfterLambda(laneSize, 0);
				directLambda(*iterator, stateAfterLambda);
				int weightOut = getWeight(stateAfterLambda);
				int curWeight = baseWeight + weightOut;
				if (curNrRounds == (int)nrRounds) {

					if ((curWeight <= maxTotalWeight)) {
						Trail newTrail(trail);
						newTrail.append((stateAfterLambda), weightOut);
						trailsOut.fetchTrail(newTrail);
					}
				}
				else {
					if (weightOut <= maxWeightOut) {
						Trail newTrail(trail);
						newTrail.append((stateAfterLambda), weightOut);
						recurseForwardExtendTrail(newTrail, trailsOut, nrRounds, maxTotalWeight);
					}
				}
				++progress;
			}
		}
		progress.unstack();
	}
}
