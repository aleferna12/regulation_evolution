/*

Copyright 1996-2006 Roeland Merks

This file is part of Tissue Simulation Toolkit.

Tissue Simulation Toolkit is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

Tissue Simulation Toolkit is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tissue Simulation Toolkit; if not, write to the Free
Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
02110-1301 USA

*/
#include <math.h>
#include <cstdlib>
#include <functional> // for bind, see InitIncreaseVal()
#include "crash.h"
#include "parameter.h"
#include "ca.h"
#include "intplane.h"
#include "misc.h"

using std::placeholders::_1;

extern Parameter par;

/** PRIVATE **/

IntPlane::IntPlane(const int sx, const int sy) {
  sigma = nullptr;
  sizex = sx;
  sizey = sy;
  diagonal = sqrt(2 * sizex * sizey);

  sigma = new int[sizex * sizey];
}


IntPlane::IntPlane() {
  sigma = nullptr;
  sizex = 0;
  sizey = 0;
  diagonal = 0;
}

// destructor (virtual)
IntPlane::~IntPlane() {
  if (sigma) {
    delete[] sigma;
    sigma = nullptr;
  }
}

void IntPlane::Plot(Graphics *g2) {
  // l=layer: default layer is 0
  for (int x = 1; x < sizex - 1; x++)
    for (int y = 1; y < sizey - 1; y++) {
      // Make the pixel four times as large
      // to fit with the CPM plane
      g2->Point(Sigma(x, y), 2 * x, 2 * y);
      g2->Point(Sigma(x, y), 2 * x + 1, 2 * y);
      g2->Point(Sigma(x, y), 2 * x, 2 * y + 1);
      g2->Point(Sigma(x, y), 2 * x + 1, 2 * y + 1);
    }

}

// Plot the value of the intplane only in the medium of the CPM
void IntPlane::Plot(Graphics *g2, CellularPotts *cpm) {

  // this has to take into account stuff from cpm plane (maybe x,y info should give a tau of the cpm plane)

  // cpm->Sigma(x, y) returns sigma, which I can use to indicise the vector of cells... can I?
  // not really, food doesn't know about cell

  // maybe this function should really be in dish...

  // suspend=true suspends calling of DrawScene
  for (int x = 1; x < sizex - 1; x++)
    for (int y = 1; y < sizey - 1; y++)
      //if (cpm->Sigma(x,y)==0) {
      if (Sigma(x, y) != 0) {
        int colorindex;
        if (Sigma(x, y) > 0) colorindex = 10 + Sigma(x, y);
        else colorindex = 5;
        // Make the pixel four times as large
        // to fit with the CPM plane
        g2->Point(colorindex, 2 * x, 2 * y);
        g2->Point(colorindex, 2 * x + 1, 2 * y);
        g2->Point(colorindex, 2 * x, 2 * y + 1);
        g2->Point(colorindex, 2 * x + 1, 2 * y + 1);
      }
}

//copy of this function in ca.cpp
int IntPlane::SetNextVal(int pos) {
  //the plane has a 1 px boundary on all size, therefore we place the pixels
  //within that boundary
  static int xcount = 1, ycount = 1;

  if (xcount >= sizex - 1 || ycount >= sizey - 1) {
    return 1;
  }

  setSigma(xcount, ycount, pos);
  ycount++;
  if (ycount == sizey - 1) {
    ycount = 1;
    xcount++;
  }
  return 0;
}


