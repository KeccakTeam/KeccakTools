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

#ifndef _PROGRESS_H_
#define _PROGRESS_H_

#include <time.h>
#include <string>
#include <vector>
#include "types.h"
#include <algorithm>

using namespace std;

class ProgressMeter {
public:
    vector<string> synopsis;
    vector<UINT64> index, count;
    unsigned int height;
    UINT64 topIndex;
    time_t previousDisplay;
    unsigned int lastHeightDisplayed;
    unsigned int nrDisplaysSinceFullDisplay;
public:
    ProgressMeter();
    void stack(UINT64 aCount = 0);
    void stack(const string& aSynopsis, UINT64 aCount = 0);
    void unstack();
    void operator++();
    void clear();
protected:
    void display();
    void displayIfNecessary();
};

#endif
