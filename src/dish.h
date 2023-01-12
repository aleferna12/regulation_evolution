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

#define MIGRATE 1
#define DIVIDE 2


class Dish {
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

    //! Master function that should be modified whenever we add new subplots or change the colortable format.
    void makePlots(int Time, Graphics *g);

    void plotChemPlane(Graphics *g, int start_index, int n_colors) const;

    void plotFoodPLane(Graphics *g, int color_index) const;

    //! makePlots the food for each cell (but differently for each cell type).
    //! \param start_index: First color of the two color gradients in the colortable (they must be in tandem).
    //! \param n_colors: Size of both gradients summed.
    void plotCellFood(Graphics *g, int start_index, int n_colors);

    //! makePlots migrating and dividing cells.
    void plotCellTau(Graphics *g, int div_index, int mig_index);

    void plotCellVectors(Graphics *g);

    //! Used to decide whether a cell border should be drawn at this position.
    void drawCellBorderIfNeeded(Graphics *g, int i, int j) const;

    static int CalculateJwithMedium(vector<int> key);

    static int CalculateJfromKeyLock(vector<int> key1, vector<int> lock1, vector<int> key2, vector<int> lock2);

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

    //! \brief. Returns the summed area of all cells in the dish
    int Area() const;

    //! \brief Returns the summed of all cells target area in the dish
    int TargetArea() const;

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

    static double DetermineMinDist(int n);

    double distMostIsolatedPoint();

    void updateChemPlane();

    int addFPatch(int x, int y);

    int addRandomFPatch();

    void removeFPatch(int id);

    int getFoodLeft() {
        int food = 0;
        for (auto &fp: fpatches)
            food += fp.getFoodLeft();
        return food;
    }

    PDE *PDEfield;
    IntPlane *ChemPlane;
    IntPlane *FoodPlane;
    CellularPotts *CPM;

    void saveLattice(int Time) const;
    void readLattice();

    //! Saves information about the cells as a CSV file in the directory specified by par.datadir
    int saveCellData(int Time);
    int readCellData();

    void saveFoodData(int Time);
    int readFoodData();

protected:
    //! The cells in the Petri dish; accessible to derived classes
    std::vector<Cell> cell;

    int sizex, sizey;
    // Number of resource sources
    int grad_sources;
    vector<FoodPatch> fpatches;
};

#define INIT void Dish::Init(void)

#endif
