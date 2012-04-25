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

//#include <cstdlib>
#include <iostream>
//#include <sstream>
//#include <time.h>
#include "progress.h"

ProgressMeter::ProgressMeter()
: height(0), previousDisplay(0), lastHeightDisplayed(0), nrDisplaysSinceFullDisplay(0)
{
}

void ProgressMeter::clear()
{
    height = 0;
    previousDisplay = 0;
    lastHeightDisplayed = 0;
    nrDisplaysSinceFullDisplay = 0;
    index.clear();
    count.clear();
    synopsis.clear();
}

void ProgressMeter::stack(UINT64 aCount)
{
    stack("", aCount);
}

void ProgressMeter::stack(const string& aSynopsis, UINT64 aCount)
{
    if (height > 0)
        index.push_back(topIndex);
    count.push_back(aCount);
    synopsis.push_back(aSynopsis);
    height++;
    topIndex = 0;
}

void ProgressMeter::unstack()
{
    if (height > 0) {
        if (height > 1) {
            topIndex = index[index.size()-1];
            index.pop_back();
        }
        count.pop_back();
        synopsis.pop_back();
        height--;
    }
    if (lastHeightDisplayed > height)
        lastHeightDisplayed = height;
}

void ProgressMeter::operator++()
{
    topIndex++;
    displayIfNecessary();
}

void ProgressMeter::displayIfNecessary()
{
    if (difftime(time(NULL), previousDisplay) >= 10.0)
        display();
}

void ProgressMeter::display()
{
    if (height > 0) {
        unsigned int startHeight = max(int(lastHeightDisplayed)-1, 0);
        if (startHeight >= height)
            startHeight = height-1;
        unsigned int effectiveStartHeight = startHeight;
        if (nrDisplaysSinceFullDisplay >= 100)
            effectiveStartHeight = 0;
        for(unsigned int i=effectiveStartHeight; i<height; i++) {
            for(unsigned int j=0; j<i; j++)
                cout << "  ";
            if (i < startHeight) cout << "(";
            if (synopsis[i].length() > 0)
                cout << synopsis[i] << ": ";
            if (i == (height-1))
                cout << dec << topIndex;
            else
                cout << dec << index[i];
            if (count[i] > 0)
                cout << " / " << dec << count[i];
            if (i < startHeight) cout << ")";
            cout << endl;
        }
        lastHeightDisplayed = height;
        previousDisplay = time(NULL);
        if (effectiveStartHeight > 0)
            nrDisplaysSinceFullDisplay++;
        else
            nrDisplaysSinceFullDisplay = 0;
    }
}
