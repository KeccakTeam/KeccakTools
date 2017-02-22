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

#include <iostream>
#include <sstream>
#include "Keccak-fDisplay.h"

using namespace std;

void getDisplayMap(const vector<SliceValue>& state, vector<unsigned int>& displayMap)
{
    bool optimizeEmptySlices = (state.size() >= 8);
    if (optimizeEmptySlices) {
        unsigned int z=0;
        while(z < state.size()) {
            while((z < state.size()) && (state[z] != 0))
                z++;
            displayMap.push_back(z);
            while((z < state.size()) && (state[z] == 0))
                z++;
            displayMap.push_back(z);
        }
    }
    else {
        displayMap.push_back(state.size());
        displayMap.push_back(state.size());
    }
}

inline unsigned int getDisplayNumberOfSpaces(const unsigned int& delta)
{
    return (delta <= 3) ? delta : ((delta < 10) ? 3 : 4);
}

void displayPlane(ostream& fout, const vector<SliceValue>& state, int offset, unsigned int y, const vector<unsigned int>& displayMap)
{
    unsigned int z=0;
    for(unsigned int i=0; i<displayMap.size(); i+=2) {
        for( ; z<displayMap[i]; z++) {
            RowValue row = getRowFromSlice(state[z], y);
            for(unsigned int sx=0; sx<5; sx++) {
                unsigned int x = KeccakF::index(sx-offset);
                if ((row & (1<<x)) != 0)
                    fout << "X";
                else
                    if ((x == 0) && (y == 0) && (z == 0))
                        fout << "+";
                    else
                        fout << ".";
            }
            if (z < (displayMap[i]-1))
                fout << "   ";
            else if (z < (state.size()-1))
                fout << " ";
        }
        if (displayMap[i] < displayMap[i+1]) {
            int delta = displayMap[i+1] - displayMap[i];
            if (y == 0) {
                ostringstream sout;
                if (delta == 1)
                    sout << "z";
                else if (delta == 2)
                    sout << "zz";
                else
                    sout << "z^" << dec << delta;
                fout << sout.str();
            }
            else {
                unsigned int numberOfSpaces = getDisplayNumberOfSpaces(displayMap[i+1] - displayMap[i]);
                for(unsigned int j=0; j<numberOfSpaces; j++)
                    fout << " ";
            }
            z = displayMap[i+1];
            if (z < state.size())
                fout << " ";
        }
    }
}

void displayNothing(ostream& fout, const vector<SliceValue>& state, const vector<unsigned int>& displayMap)
{
    unsigned int z=0;
    for(unsigned int i=0; i<displayMap.size(); i+=2) {
        for( ; z<displayMap[i]; z++) {
            fout << "     "; // 5 spaces
            if (z < (displayMap[i]-1))
                fout << "   ";
            else if (z < (state.size()-1))
                fout << " ";
        }
        if (displayMap[i] < displayMap[i+1]) {
            unsigned int numberOfSpaces = getDisplayNumberOfSpaces(displayMap[i+1] - displayMap[i]);
            for(unsigned int j=0; j<numberOfSpaces; j++)
                fout << " ";
            z = displayMap[i+1];
            if (z < state.size())
                fout << " ";
        }
    }
}

void displayParity(ostream& fout, const vector<SliceValue>& state, const int& offset, const vector<unsigned int>& displayMap)
{
    unsigned int z=0;
    for(unsigned int i=0; i<displayMap.size(); i+=2) {
        for( ; z<displayMap[i]; z++) {
            RowValue parity = 0;
            for(unsigned int y=0; y<5; y++)
                parity ^= getRowFromSlice(state[z], y);
            for(unsigned int sx=0; sx<5; sx++) {
                unsigned int x = KeccakF::index(sx-offset);
                if ((parity & (1<<x)) != 0)
                    fout << "O";
                else
                    fout << "-";
            }
            if (z < (displayMap[i]-1))
                fout << "   ";
            else if (z < (state.size()-1))
                fout << " ";
        }
        if (displayMap[i] < displayMap[i+1]) {
            unsigned int numberOfSpaces = getDisplayNumberOfSpaces(displayMap[i+1] - displayMap[i]);
            for(unsigned int j=0; j<numberOfSpaces; j++)
                fout << " ";
            z = displayMap[i+1];
            if (z < state.size())
                fout << " ";
        }
    }
}

void displayState(ostream& fout, const vector<SliceValue>& state, bool showParity)
{
    const int offset = 2;
    vector<unsigned int> displayMap;
    getDisplayMap(state, displayMap);
    for(unsigned int sy=0; sy<5; sy++) {
        unsigned int y = KeccakF::index(-1-sy-offset);
        displayPlane(fout, state, offset, y, displayMap);
        fout << endl;
    }
    if (showParity) {
        displayParity(fout, state, offset, displayMap);
        fout << endl;
    }
}

void displaySlice(ostream& fout, SliceValue slice)
{
    vector<SliceValue> slices;
    slices.push_back(slice);
    displayState(fout, slices);
}

void KeccakDisplayInSVG::displaySlice(ostream& fout, SliceValue slice, unsigned int z, double x, double y, double scale) const
{
    fout << "<g transform=\"translate(" << x << ", " << y << ") scale(" << scale << ")\">\n";

    bool activeSlice = (slice != 0);
    for(unsigned int iy=0; iy<5; iy++) {
        unsigned int y = KeccakF::index(7-iy);
        bool activeRow = (getRowFromSlice(slice, y) != 0);
        for(unsigned int ix=0; ix<5; ix++) {
            unsigned int x = KeccakF::index(ix+3);
            bool bit = getRowFromSlice(slice, y) & (1 << x);
            fout << "<rect x=\"" << dec << (ix*bitSize)
                << "\" y=\"" << (iy*bitSize)
                << "\" width=\"" << bitSize
                << "\" height=\"" << bitSize << "\" class=\"bit";
            if (bit)
                fout << " bit_bit_active";
            if (displaySlices && activeSlice)
                fout << " bit_slice_active";
            if (displayRows && activeRow)
                fout << " bit_row_active";
            fout << " " << getAdditionalBitStyles(x, y, z);
            fout << "\"/>\n";
        }
    }
    for(unsigned int iy=0; iy<5; iy++) {
        unsigned int y = KeccakF::index(7-iy);
        bool activeRow = (getRowFromSlice(slice, y) != 0);
        if (displayRows) {
            fout << "<rect x=\"" << dec << 0
                << "\" y=\"" << (iy*bitSize)
                << "\" width=\"" << (5*bitSize)
                << "\" height=\"" << bitSize << "\" class=\"f_row";
            if (activeRow)
                fout << " f_row_row_active";
            if (activeSlice)
                fout << " f_row_slice_active";
            fout << "\"/>\n";
        }
    }
    if (displaySlices) {
        fout << "<rect x=\"" << dec << 0
            << "\" y=\"" << 0
            << "\" width=\"" << (5*bitSize)
            << "\" height=\"" << (5*bitSize) << "\" class=\"f_slice";
        if (activeSlice)
            fout << " f_slice_slice_active";
        fout << "\"/>\n";
    }
    fout << "</g>" << endl;
}


void KeccakDisplayInSVG::displayState(ostream& fout, const vector<SliceValue>& state, double gOffsetX, double gOffsetY, double gScale) const
{
    const double zOffsetX = 0.4*bitSize;
    const double zOffsetY = -0.3*bitSize;
    double offsetY = (zOffsetY < 0.0) ? (-(int)(state.size() - 1) * zOffsetY) : 0.0;

    fout << "<g transform=\"translate(" << gOffsetX << ", " << gOffsetY << ") scale(" << gScale << ")\">\n";

    for(unsigned int iz=0; iz<state.size(); iz++) {
        unsigned int z = state.size()-1-iz;
        displaySlice(fout, state[z], z, z*zOffsetX, offsetY + z*zOffsetY);
    }
    fout << "</g>" << endl;
}

void KeccakDisplayInSVG::displayStateSparsely(ostream& fout, const vector<SliceValue>& state, double gOffsetX, double gOffsetY, double gScale) const
{
    fout << "<g transform=\"translate(" << gOffsetX << ", " << gOffsetY << ") scale(" << gScale << ")\">\n";

    unsigned int activeSlices = 0;
    for(unsigned int z=0; z<state.size(); z++)
        if (state[z] != 0) activeSlices++;
    unsigned int rows = (maxNumberOfHorizontalSlices <= 0) ? 1 : ((activeSlices+maxNumberOfHorizontalSlices-1) / maxNumberOfHorizontalSlices);
    unsigned int slicesPerRow = (activeSlices+rows-1) / rows;
    double px = 0.0;
    double py = 0.0;
    int j = 0;
    for(unsigned int z=0; z<state.size(); z++) {
        if (state[z] != 0) {
            displaySlice(fout, state[z], z, px, py);
            fout << "<text xml:space=\"preserve\" class=\"normal\" x=\"" << dec << (px+2.5*bitSize)
                << "\" y=\"" << (py-0.2*bitSize)
                << "\" text-anchor=\"middle\"><tspan style=\"font-style:italic;\">z</tspan> = " << z
                << "</text>\n";
            px += 7*bitSize;
            j++;
            if (j >= (int)slicesPerRow) {
                px = 0;
                py += 7*bitSize;
                j = 0;
            }
        }
    }
    fout << "</g>" << endl;
}

void KeccakDisplayInSVG::displayTrail(ostream& fout, const KeccakFPropagation& DCorLC, const Trail& trail,
    double x, double y, double scale) const
{
    (void)x;
    (void)y;
    (void)scale;
    vector<vector<SliceValue> > allAfterPreviousChi;
    vector<vector<SliceValue> > allBeforeTheta;
    vector<vector<SliceValue> > allAfterTheta;

    for(unsigned int i=0; i<trail.states.size(); i++) {
        vector<SliceValue> stateAfterChi;
        DCorLC.reverseLambda(trail.states[i], stateAfterChi);
        allAfterPreviousChi.push_back(stateAfterChi);
        vector<SliceValue> stateBeforeTheta;
        DCorLC.directLambdaBeforeTheta(stateAfterChi, stateBeforeTheta);
        allBeforeTheta.push_back(stateBeforeTheta);
        vector<SliceValue> stateAfterTheta;
        DCorLC.directTheta(stateBeforeTheta, stateAfterTheta);
        allAfterTheta.push_back(stateAfterTheta);
    }

    unsigned int py = 0;
    for(unsigned int i=0; i<trail.states.size(); i++) {
        if (i > 0) {
            displayStateSparsely(fout, allAfterPreviousChi[i], 0, py);
            fout << "<path class=\"arrow\" d=\"M " << dec << (-bitSize) << "," << (py+5*bitSize) << " "
                << (-bitSize) << "," << (py+9*bitSize) << "\"/>\n";
            fout << "<text xml:space=\"preserve\" class=\"normal\" x=\"" << dec << (-0.8*bitSize) << "\" y=\"" << (py+7*bitSize)
                << "\"><tspan style=\"font-style:italic;\">\xCE\xB8</tspan>, "
                << "<tspan style=\"font-style:italic;\">\xCF\x81</tspan>, "
                << "<tspan style=\"font-style:italic;\">\xCF\x80</tspan></text>\n";
            py += 9*bitSize;
        }
        {
            displayStateSparsely(fout, trail.states[i], 0, py);
            fout << "<text xml:space=\"preserve\" class=\"normal\" x=\"" << dec << (-0.2*bitSize)
                << "\" y=\"" << (py+3*bitSize)
                << "\" text-anchor=\"end\">weight: " << trail.weights[i]
                << "</text>\n";
            if (i < trail.states.size()-1) {
                fout << "<path class=\"arrow\" d=\"M " << dec << (-bitSize) << "," << (py+5*bitSize) << " "
                    << (-bitSize) << "," << (py+9*bitSize) << "\"/>\n";
            fout << "<text xml:space=\"preserve\" class=\"normal\" x=\"" << dec << (-0.8*bitSize) << "\" y=\"" << (py+7*bitSize)
                    << "\"><tspan style=\"font-style:italic;\">\xCF\x87</tspan></text>\n";
                py += 9*bitSize;
            }
        }
    }
}

void KeccakDisplayInSVG::displayParity(ostream& fout, const KeccakFPropagation& DCorLC, const vector<RowValue>& C, const vector<RowValue>& D, bool displayRuns) const
{
    fout << "<g>\n";

    if (displayRuns) {
        bool previous = false;
        for(unsigned int t=0; t<DCorLC.laneSize*5; t++) {
            unsigned int x, z;
            DCorLC.getXandZfromT(t, x, z);
            unsigned int ix = (x+2)%5;
            unsigned int iz = DCorLC.laneSize-1-z;
            bool odd = ((C[z] & (1 << x)) != 0);
            if (odd) {
                if (!previous)
                    fout << "<path class=\"run\" d=\"M";
                fout << " " << ((ix+0.5)*bitSize) << "," << ((iz+0.5)*bitSize);
            }
            else {
                if (previous)
                    fout << "\"/>\n";
            }
            previous = odd;
        }
    }
    for(unsigned int t=0; t<DCorLC.laneSize*5; t++) {
        unsigned int x, z;
        DCorLC.getXandZfromT(t, x, z);
        unsigned int ix = (x+2)%5;
        unsigned int iz = DCorLC.laneSize-1-z;
        bool odd = ((C[z] & (1 << x)) != 0);
        bool affected = ((D[z] & (1 << x)) != 0);
        fout << "<rect x=\"" << dec << (ix*bitSize)
            << "\" y=\"" << (iz*bitSize)
            << "\" width=\"" << bitSize
            << "\" height=\"" << bitSize << "\" class=\"column";
        if (odd)
            fout << " column_odd";
        if (affected)
            fout << " column_affected";
        fout << "\"/>\n";
        if (odd)
            fout << "<circle cx=\"" << ((ix+0.5)*bitSize) << "\" cy=\"" << ((iz+0.5)*bitSize) << "\" r=\"" << (bitSize*0.35)
                << "\" class=\"odd\"/>\n";
        if (affected)
            fout << "<circle cx=\"" << ((ix+0.5)*bitSize) << "\" cy=\"" << ((iz+0.5)*bitSize) << "\" r=\"" << (bitSize*0.2)
                << "\" class=\"affected\"/>\n";
    }
    fout << "<path class=\"arrow\" d=\"M " << dec << (-bitSize) << "," << (DCorLC.laneSize*bitSize) << " "
        << (-bitSize) << "," << ((DCorLC.laneSize-3.0)*bitSize) << "\"/>\n";
    fout << "<text xml:space=\"preserve\" class=\"normal\" x=\"" << dec << (-bitSize-3) << "\" y=\"" << ((DCorLC.laneSize-3.0)*bitSize-5)
        << "\"><tspan style=\"font-style:italic;\">z</tspan></text>\n";
    fout << "<path class=\"arrow\" d=\"M " << dec << 0 << "," << ((DCorLC.laneSize+1)*bitSize) << " "
        << (3*bitSize) << "," << ((DCorLC.laneSize+1)*bitSize) << "\"/>\n";
    fout << "<text xml:space=\"preserve\" class=\"normal\" x=\"" << dec << (3*bitSize+3) << "\" y=\"" << ((DCorLC.laneSize+1)*bitSize+3)
        << "\"><tspan style=\"font-style:italic;\">x</tspan></text>\n";
    fout << "</g>" << endl;
}

string KeccakDisplayInSVG::getAdditionalBitStyles(unsigned int x, unsigned int y, unsigned int z) const
{
    (void)x;
    (void)y;
    (void)z;
    return "";
}

void KeccakDisplayInSVG::displayRow(ostream& fout)
{
    {
        unsigned int iz = laneSize-1;
        unsigned int z = laneSize-1-iz;
        {
            unsigned int iy = 0;
            unsigned int y = KeccakF::index(7-iy);
            for(unsigned int ix=0; ix<5; ix++) {
                unsigned int x = KeccakF::index(ix+3);
                double px, py;
                getPosition(ix, iy, z, px, py);
                fout << "<rect x=\"" << dec << px
                    << "\" y=\"" << py
                    << "\" width=\"" << bitSize
                    << "\" height=\"" << bitSize << "\" class=\"bit";
                fout << " " << getAdditionalBitStyles(x, y, z);
                fout << "\"/>\n";
            }
        }
    }
}

void KeccakDisplayInSVG::displayColumn(ostream& fout)
{
    {
        unsigned int iz = laneSize-1;
        unsigned int z = laneSize-1-iz;
        for(unsigned int iy=0; iy<5; iy++) {
            unsigned int y = KeccakF::index(7-iy);
            {
                unsigned int ix = 0;
                unsigned int x = KeccakF::index(ix+3);
                double px, py;
                getPosition(ix, iy, z, px, py);
                fout << "<rect x=\"" << dec << px
                    << "\" y=\"" << py
                    << "\" width=\"" << bitSize
                    << "\" height=\"" << bitSize << "\" class=\"bit";
                fout << " " << getAdditionalBitStyles(x, y, z);
                fout << "\"/>\n";
            }
        }
    }
}

void KeccakDisplayInSVG::displayLane(ostream& fout)
{
    for(unsigned int iz=0; iz<(unsigned int)laneSize; iz++) {
        unsigned int z = laneSize-1-iz;
        {
            unsigned int iy = 0;
            unsigned int y = KeccakF::index(7-iy);
            {
                unsigned ix=0;
                unsigned int x = KeccakF::index(ix+3);
                double px, py;
                getPosition(ix, iy, z, px, py);
                fout << "<rect x=\"" << dec << px
                    << "\" y=\"" << py
                    << "\" width=\"" << bitSize
                    << "\" height=\"" << bitSize << "\" class=\"bit";
                fout << " " << getAdditionalBitStyles(x, y, z);
                fout << "\"/>\n";
            }
        }
    }
}

void KeccakDisplayInSVG::displayPlane(ostream& fout)
{
    for(unsigned int iz=0; iz<(unsigned int)laneSize; iz++) {
        unsigned int z = laneSize-1-iz;
        {
            unsigned int iy = 0;
            unsigned int y = KeccakF::index(7-iy);
            for(unsigned int ix=0; ix<5; ix++) {
                unsigned int x = KeccakF::index(ix+3);
                double px, py;
                getPosition(ix, iy, z, px, py);
                fout << "<rect x=\"" << dec << px
                    << "\" y=\"" << py
                    << "\" width=\"" << bitSize
                    << "\" height=\"" << bitSize << "\" class=\"bit";
                fout << " " << getAdditionalBitStyles(x, y, z);
                fout << "\"/>\n";
            }
        }
    }
}

void KeccakDisplayInSVG::displaySlice(ostream& fout)
{
    {
        unsigned int iz = laneSize-1;
        unsigned int z = laneSize-1-iz;
        for(unsigned int iy=0; iy<5; iy++) {
            unsigned int y = KeccakF::index(7-iy);
            for(unsigned int ix=0; ix<5; ix++) {
                unsigned int x = KeccakF::index(ix+3);
                double px, py;
                getPosition(ix, iy, z, px, py);
                fout << "<rect x=\"" << dec << px
                    << "\" y=\"" << py
                    << "\" width=\"" << bitSize
                    << "\" height=\"" << bitSize << "\" class=\"bit";
                fout << " " << getAdditionalBitStyles(x, y, z);
                fout << "\"/>\n";
            }
        }
    }
}

void KeccakDisplayInSVG::displaySheet(ostream& fout)
{
    for(unsigned int iz=0; iz<(unsigned int)laneSize; iz++) {
        unsigned int z = laneSize-1-iz;
        for(unsigned int iy=0; iy<5; iy++) {
            unsigned int y = KeccakF::index(7-iy);
            {
                unsigned ix=0;
                unsigned int x = KeccakF::index(ix+3);
                double px, py;
                getPosition(ix, iy, z, px, py);
                fout << "<rect x=\"" << dec << px
                    << "\" y=\"" << py
                    << "\" width=\"" << bitSize
                    << "\" height=\"" << bitSize << "\" class=\"bit";
                fout << " " << getAdditionalBitStyles(x, y, z);
                fout << "\"/>\n";
            }
        }
    }
}

void KeccakDisplayInSVG::displayState(ostream& fout)
{
    for(unsigned int iz=0; iz<(unsigned int)laneSize; iz++) {
        unsigned int z = laneSize-1-iz;
        for(unsigned int iy=0; iy<5; iy++) {
            unsigned int y = KeccakF::index(7-iy);
            for(unsigned int ix=0; ix<5; ix++) {
                unsigned int x = KeccakF::index(ix+3);
                double px, py;
                getPosition(ix, iy, z, px, py);
                fout << "<rect x=\"" << dec << px
                    << "\" y=\"" << py
                    << "\" width=\"" << bitSize
                    << "\" height=\"" << bitSize << "\" class=\"bit";
                fout << " " << getAdditionalBitStyles(x, y, z);
                fout << "\"/>\n";
            }
        }
    }
}

void KeccakDisplayInSVG::getPosition(unsigned int x, unsigned int y, unsigned int z, double &px, double &py) const
{
    const double zOffsetX = 0.4*bitSize;
    const double zOffsetY = -0.3*bitSize;

    px = (x*bitSize) + z*zOffsetX;
    py = (y*bitSize) + z*zOffsetY;
}

void displayStates(ostream& fout,
                   const vector<SliceValue>& state1, bool showParity1,
                   const vector<SliceValue>& state2, bool showParity2)
{
    const int offset = 2;
    vector<unsigned int> displayMap1, displayMap2;
    getDisplayMap(state1, displayMap1);
    getDisplayMap(state2, displayMap2);
    for(unsigned int sy=0; sy<5; sy++) {
        unsigned int y = KeccakF::index(-1-sy-offset);
        displayPlane(fout, state1, offset, y, displayMap1);
        fout << "  |  ";
        displayPlane(fout, state2, offset, y, displayMap2);
        fout << endl;
    }
    if (showParity1 || showParity2) {
        if (showParity1)
            displayParity(fout, state1, offset, displayMap1);
        else
            displayNothing(fout, state1, displayMap1);
        if (showParity2) {
            fout << "     ";
            displayParity(fout, state2, offset, displayMap2);
        }
        fout << endl;
    }
}

void displayStates(ostream& fout,
                   const vector<SliceValue>& state1, bool showParity1,
                   const vector<SliceValue>& state2, bool showParity2,
                   const vector<SliceValue>& state3, bool showParity3)
{
    const int offset = 2;
    vector<unsigned int> displayMap1, displayMap2, displayMap3;
    getDisplayMap(state1, displayMap1);
    getDisplayMap(state2, displayMap2);
    getDisplayMap(state3, displayMap3);
    for(unsigned int sy=0; sy<5; sy++) {
        unsigned int y = KeccakF::index(-1-sy-offset);
        displayPlane(fout, state1, offset, y, displayMap1);
        fout << "  |  ";
        displayPlane(fout, state2, offset, y, displayMap2);
        fout << "  |  ";
        displayPlane(fout, state3, offset, y, displayMap3);
        fout << endl;
    }
    if (showParity1 || showParity2 || showParity3) {
        if (showParity1)
            displayParity(fout, state1, offset, displayMap1);
        else
            displayNothing(fout, state1, displayMap1);
        fout << "     ";
        if (showParity2)
            displayParity(fout, state2, offset, displayMap2);
        else
            displayNothing(fout, state2, displayMap2);
        if (showParity3) {
            fout << "     ";
            displayParity(fout, state3, offset, displayMap3);
        }
        fout << endl;
    }
}
