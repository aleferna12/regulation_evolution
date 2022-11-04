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

#ifndef _INTPL_HH_
#define _INTPL_HH_

#include <stdio.h>
#include <float.h>
#include <functional>
#include <climits>
#include "graph.h"
#include "random.h"
#include "foodpatch.h"

class IntPlane; //forward declaration

class CellularPotts;

class IntPlane {

    friend class Info;

public:

    /*!
     * \brief Constructor for PDE object containing arbitrary number of planes.
     * \param layers: Number of PDE planes
     * \param sx: horizontal size of PDE planes
     * \param sy: vertical size of PDE planes
    */

    IntPlane(const int sx, const int sy);


    // destructor must also be virtual
    virtual ~IntPlane();

    /*!
     * \brief Plots one layer of the PDE plane to a Graphics window.
     * \param g: Graphics window.
     * \param layer: The PDE plane to be plotted. Default layer 0.
    */
    void Plot(Graphics *g);

    /*! \brief Plots one layer of the PDE to a Graphics window, but not over the cells.
      \param g: Graphics window.
      \param cpm: CellularPotts object containing the cells.
      \param layer: The PDE plane to be plotted. Default layer 0.
    */

    void Plot(Graphics *g, CellularPotts *cpm);

    /*! \brief Plots the PDE field using contour lines.

    \param g: Graphics window.
    \param layer: The PDE plane to be plotted. Default layer 0.
    \param colour: Color to use for the contour lines, as defined in the "default.ctb" color map file, which should be in the same directory as the executable. Default color 1 (black in the default color map).
    */

    //! \brief Returns the horizontal size of the PDE planes.
    inline int SizeX() const {
      return sizex;
    }

    //! \brief Returns the vertical size of the PDE planes.
    inline int SizeY() const {
      return sizey;
    }

    /*! \brief Returns the value of grid point x,y of PDE plane "layer".

    Warning, no range checking done.

    \param layer: the PDE plane to probe.
    \param x, y: grid point to probe.
    */
    inline int Sigma(const int x, const int y) const {
      return sigma[x * sizex + y];
    }

    /*! \brief Sets grid point x,y of PDE plane "layer" to value "value".

    \param layer: PDE plane.
    \param x, y: grid point
    \param value: new contents

    */
    inline void setSigma(const int x, const int y, const int value) {
      sigma[x * sizex + y] = value;
    }

    int SetNextVal(int sig);

    pair<int, int> getMinMax() const {
      int maxval = 0;
      int minval = INT_MAX;
      for (int i = 1; i < sizex; ++i) for (int j = 1; j < sizex; ++j) {
        minval = min(Sigma(i, j), minval);
        maxval = max(Sigma(i, j), maxval);
      }
      return {minval, maxval};
    }

protected:
    // Protected member functions

    /*! \brief Used in Plot. Takes a color and turns it into a grey value.

    \param val: Value from PDE plane.

    Implement this function in you main simulation code. See e.g. vessel.cpp.
    */
    //virtual int MapColour(double val);

    //! empty constructor (necessary for derivation)
    IntPlane();

private:
    int *sigma;
    int sizex;
    int sizey;
    // Diagonal length of the lattice (calculated)
    double diagonal;
};

#endif
