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

#ifndef _TRANSLATIONSYMMETRY_H_
#define _TRANSLATIONSYMMETRY_H_

#include <vector>

/** This function defines an order between vectors.
  * @param  a   The first vector to compare.
  * @param  b   The second vector to compare.
  * @return True iff @a a is smaller than @a b.
  */
template<class T>
bool isSmaller(const std::vector<T>& a, const std::vector<T>& b)
{
    for(unsigned int iz=0; iz<a.size(); iz++) {
        unsigned int z = a.size()-1-iz;
        if (a[z] < b[z])
            return true;
        else if (a[z] > b[z])
            return false;
    }
    return false;
}

/** This function returns whether the given vector is smaller
  * after translation than itself.
  * @param  a   The vector to compare.
  * @param  dz  The amount of translation.
  * @return True iff Translate(@a a, @a dz) is smaller than @a a.
  */
template<class T>
bool isSmallerAfterTranslation(const std::vector<T>& a, unsigned int dz)
{
    unsigned int size = a.size();
    unsigned int z = size-1;
    while((z > 0) && (a[z] == a[(z+dz)%size]))
        z--;
    return (a[z] > a[(z+dz)%size]);
}

/** This function returns whether the given vector without translation
  * is the smallest among the translated versions of itself.
  * @param  a   The vector to test.
  * @return True iff @a a is smaller than Translate(@a a, @a dz) for all @a dz ≠ 0.
  */
template<class T>
bool isMinimalSymmetrically(const std::vector<T>& a)
{
    unsigned int laneSize = a.size();
    for(unsigned int dz=1; dz<laneSize; dz++)
        if (isSmallerAfterTranslation(a, dz))
            return false;
    return true;
}

/** This function returns the minimum among the translated
  * versions of the given vector.
  * @param  a   The vector to translate.
  * @param  aMin    The resulting minimum vector, i.e., such that
  *     isMinimalSymmetrically(aMin) is true.
  */
template<class T>
void getSymmetricMinimum(const std::vector<T>& a, std::vector<T>& aMin)
{
    aMin = a;
    unsigned int laneSize = a.size();
    for(unsigned int dz=1; dz<laneSize; dz++) {
        std::vector<T> aDz(a.size());
        for(unsigned int z=0; z<laneSize; z++)
            aDz[(z+dz)%laneSize] = a[z];
        if (isSmaller(aDz, aMin))
            aMin = aDz;
    }
}

#endif
