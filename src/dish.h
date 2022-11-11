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

/*! \class Dish
  \brief The virtual Petri dish.

  Hosts the cells with states and the CA-plane.
*/

#ifndef CRITTER_H_
#define CRITTER_H_

#include <vector>
#include <set>
#include "graph.h"
#include "random.h"
#include "pde.h"
#include "intplane.h"
#include "cell.h"
#include "ca.h"
#include "foodpatch.h"

#define PREY 1
#define PREDATOR 2

namespace ColourMode {
    enum {
        State, CellType, Sigma, Auxilliary
    };
}

class Dish {

    friend class Info;

public:
    Dish();

    /*! \brief Init defines the initial state of the virtual
      cell culture.

      Define Init() in your main file describing the simulation set up,
      within the block INIT { }. See for examples vessel.cpp and
      sorting.cpp.
    */
    void Init();

    virtual ~Dish();

    /*! \brief Plot the Dish to graphics window g.

    Simply calls CPM->Plot.
    And also colors depnd on food.
    */
    void ChemPlot(Graphics *g) const;

    void Plot(Graphics *g, int colour);

    void InitKeyLock();

    static int CalculateJwithMedium(vector<int> key);

    static int CalculateJfromKeyLock(vector<int> key1, vector<int> lock1, vector<int> key2, vector<int> lock2);

    void InitVectorJ(); //Initialise vector of J values for each cell
    void UpdateVectorJ(const vector<int> &sigma_to_update);

    void MutateCells(const vector<int> &sigma_to_update);

    void InitContactLength();

    void UpdateNeighDuration();

    //! Returns the number of completed Monte Carlo Steps.
    int Time() const;

    //! Returns the number of cells in the dish, excluding apoptosed cells.
    int CountCells() const;

//count how many cells of each group there are (when running competition), return 1 when one is extinct
    int CountCellGroups() const;

    void CellsEat(int time); // Based on the old CellsEat2

    void InitCellMigration();

    void CellMigration();

    void UpdateCellParameters(int Time);

    void GradientBasedCellKill(int currentsize);

    //! \brief. Returns the summed area of all cells in the dish
    int Area() const;

    //! \brief Returns the summed of all cells target area in the dish
    int TargetArea() const;

    int SaveData(int Time);

    void SaveNetworks(int Time);

    void SaveAdheringNeighbours(int Time);

    void SaveAncestry(int Time);

    void MakeBackup(int Time);

    // TODO: Test if it is working with multiple fpatches
    int ReadBackup(char *filename);

    int ReadCompetitionFile(char *filename);

    //! \brief Returns the horizontal size of the dish.
    int SizeX() const;

    //! \brief Returns the horizontal size of the dish.
    int SizeY() const;

    //! \brief Returns a reference to cell number "c"
    inline Cell &getCell(int c) {
      return cell[c];
    }

    const FoodPatch &getFPatch(int id) {
      return fpatches[id];
    }

    inline int GetGradSources() const {
      return grad_sources;
    }

    // Get distance of coordinate to the closest source of resources
    pair<int, double> closestFPatch(int x, int y);

    // ChemPlane in relation to the distance from a single peak F(x)
    double FoodEquation(double dist_from_peak) const;

    // Determines how much food is in a specific position
    double FoodAtPosition(int x, int y);

    // Iterates lattice and finds max food
    // Writes a few rows of the sigma lattice to par.peaksdatafile
    int WritePeaksData() const;

    static double DetermineMinDist(int n) ;

    double distMostIsolatedPoint();

    void updateChemPlane();

    int addFPatch(int x, int y);

    int addRandomFPatch();

    void removeFPatch(int id);

    int getFoodLeft() {
      int food = 0;
      for (auto &fp : fpatches)
        food += fp.getFoodLeft();
      return food;
    }

    PDE *PDEfield;
    IntPlane *ChemPlane;
    IntPlane *FoodPlane;
    CellularPotts *CPM;

protected:
    //! The cells in the Petri dish; accessible to derived classes
    std::vector<Cell> cell;

    int sizex, sizey;
    // Number of resource sources
    int grad_sources;
    vector<FoodPatch> fpatches;
    // Coefficient for how the resources decrease according to their distance from the sources
    // TODO: Should this be a parameter?
    double dist_coef;

    void FoodPlot(Graphics *g, int colori) const;

};

#define INIT void Dish::Init(void)

#endif
