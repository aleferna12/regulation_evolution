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
#include "graph.h"
#include "random.h"

class IntPlane; //forward declaration

class CellularPotts;
class IntPlane {

  friend class Info;

 public:

  /*!
   * \brief Constructor for PDE object containing arbitrary number of planes.
   * \param layers: Number of PDE planes
   * \param sizex: horizontal size of PDE planes
   * \param sizey: vertical size of PDE planes
  */

  IntPlane(const int sizex, const int sizey);


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
    return sigma[x][y];
  }

  /*! \brief Sets grid point x,y of PDE plane "layer" to value "value".

  \param layer: PDE plane.
  \param x, y: grid point
  \param value: new contents

  */
  inline void setValue(const int x, const int y, const int value) {
    sigma[x][y]=value;
  }

  /*! \brief Adds a number to a PDE grid point.

  \param layer: PDE plane.
  \param x, y: grid point
  \param value: value to add
  */
  inline void addtoValue(const int x, const int y, const int value) {
    sigma[x][y]+=value;
  }


  //IncreaseVal should be a function pointer to either
  // IncreaseValEverywhere <- which takes no arguments, or
  // IncreaseValIfEmpty, which takes CellularPotts *cpm as argument
  // a way to achieve this, is to set IncreaseVal as a function from bind
  // but what is the right way to declare it?

//   std::function<void(int)> IncreaseVal;
  // f_display = print_num; probably I can do this in InitIncreaseVal()

  typedef std::function<void(IntPlane&)> my_fun_t; //apparently this is legal declaration
  my_fun_t IncreaseVal;

  // Add values to random positions in plane (e.g. food influx)
  int SetNextVal(int sig);
  void IncreaseValEverywhere(void);
  // Add values to random positions in plane if not occupied (e.g. food influx)
  void IncreaseValIfEmpty(CellularPotts *cpm); //set if notonprey parameter
  void NotIncreaseVal(void);                   // set if nowhere
  void IncreaseValPatchy(CellularPotts *cpm);         // patchy food increase
  void IncreaseValSomewhere(CellularPotts *cpm);         // food increases only somewhere, always there
  void IncreaseValSomewhereIfEmpty(CellularPotts *cpm);
  void IncreaseValPatchyRandom(CellularPotts *cpm); // food increases in random patches that are depleted
  void IncreaseValPatchyRandomPersistence(CellularPotts *cpm);
  void IncreaseValSelfGrowth(CellularPotts *cpm);
  void IncreaseValBoundaryGrad(CellularPotts *cpm);
  void IncreaseValSpecifiedExp(CellularPotts *cpm);
  void InitIncreaseVal(CellularPotts *cpm);// set if everywhere parameter
                                           //Set function pointer for food update,
                                           // pointer to cpm is not used if IncreaseVal
                                           // points to IncreaseValEverywhere

  /*! \brief Gets the maximum value of PDE layer l.

  \param l: layer
  \return Maximum value in layer l.
  */
  inline int Max(void ) {
    int max=sigma[0][0];
    int loop=sizex*sizey;
    for (int i=1;i<loop;i++)
      if (sigma[0][i]>max) {
	max=sigma[0][i];
      }
    return max;
  }
  /*! \brief Returns the minimum value of PDE layer l.

  \param l: layer
  \return Minimum value in layer l.
  */
  inline int Min(void) {
    int min=sigma[0][0];
    int loop=sizex*sizey;
    for (int i=1;i<loop;i++)
      if (sigma[0][i]<min) {
	min=sigma[0][i];
      }
    return min;
  }

  /*! \brief Carry out $n$ diffusion steps for all PDE planes.

  We use a forward Euler method here. Can be replaced for better algorithm.

  \param repeat: Number of steps.

  Time step dt, space step dx, diffusion coefficient diff_coeff and
  boundary conditions (bool periodic_boundary) are set as global
  parameters in a parameter file using class Parameter.

  */
  void Diffuse(int repeat);

  void DiffuseParticles(void);

  /*! \brief Implementation of no-flux boundaries.

  Called internally (optionally) by Diffuse(). */
  void NoFluxBoundaries(void);

  /*! \brief Implementation of absorbing boundaries.

  Called internally (optionally) by Diffuse(). */
  void AbsorbingBoundaries(void);

  /*! \brief Implementation of periodic boundaries.

  Called internally (optionally) by Diffuse(). */
  void PeriodicBoundaries(void);

  /*! \brief Reaction and interaction of CPM plane with PDE planes.

  \param cpm: CellularPotts plane the PDE plane interacts with

  You should implement this member function as part of your main
  simulation code. See for an example vessel.cpp.

  */
  void Secrete(CellularPotts *cpm);

  /*! \brief Returns cumulative "simulated" time,
    i.e. number of time steps * dt. */
  inline double TheTime(void) const {
    return thetime;
  }

  inline int GetPeakx(void){
    return peakx;
  }
  inline int GetPeaky(void){
    return peaky;
  }
  inline void SetPeakx(int px){
    peakx=px;
  }
  inline void SetPeaky(int py){
    peaky=py;
  }

 protected:

  int **sigma;

  int sizex;
  int sizey;
  int peakx,peaky; // location of peak of gradient,
                   // for gradient experiments

  // Protected member functions

  /*! \brief Used in Plot. Takes a color and turns it into a grey value.

  \param val: Value from PDE plane.

  Implement this function in you main simulation code. See e.g. vessel.cpp.
  */
  //virtual int MapColour(double val);

  //! empty constructor (necessary for derivation)
  IntPlane(void);

  /*! \brief Allocates a PDE plane (internal use).

  For internal use, can be reimplemented in derived class to change
  method of memory allocation.
  */
  virtual int **AllocateSigma(const int sx, const int sy);

 private:
  static const int nx[9], ny[9];
  double thetime;

};

#endif
