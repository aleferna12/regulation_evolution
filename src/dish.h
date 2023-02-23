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


struct CellGravestone {
    int sigma;
    int tau;
    int time_since_birth;
    int time_death;
    double self_gamma;
    string reason;
};


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

    //! Plots the chemotactic gradient
    void plotChemPlane(Graphics *g, int start_index, int n_colors) const;

    //! Plots circles irradiating from the food patches
    //! \param bg_index: Color of the background, only used if not drawing the chemotactic gradient behind the circles
    void plotChemPlaneCircles(Graphics *g, int fg_index, int bg_index) const;

    void plotFoodPLane(Graphics *g, int color_index) const;

    //! Plots the food for each cell (but differently for each cell type).
    //! \param start_index: First color of the two color gradients in the colortable (they must be in tandem).
    //! \param n_colors: Size of both gradients summed.
    void plotCellFood(Graphics *g, int start_index, int n_colors);

    //! Plots migrating and dividing cells.
    void plotCellTau(Graphics *g, int mig_index, int div_index);

    //! Plots cell color according to their 'group' attribute. Only supports two groups right now.
    void plotCellGroup(Graphics *g, int group1_tau1, int group1_tau2, int group2_tau1, int group2_tau2);

    void plotCellVectors(Graphics *g);

    void plotCellBorders(Graphics *g);

    void MutateCells(const vector<int> &sigma_to_update);

    void InitContactLength();

    void UpdateNeighDuration();

    //! Returns the number of completed Monte Carlo Steps.
    int Time() const;

    //! Returns the number of cells in the dish, excluding apoptosed cells.
    int CountCells() const;

    //! Count how many cells of each group there are (when running competition), return 1 when one is extinct.
    //! Assumes two groups!
    bool groupExtinction() const;

    void CellsEat(int time); // Based on the old CellsEat2

    void InitCellMigration();

    void CellMigration();

    void UpdateCellParameters(int Time);

    //! \brief. Returns the summed area of all cells in the dish
    int Area() const;

    //! \brief Returns the summed of all cells target area in the dish
    int TargetArea() const;

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
    pair<int, double> closestFPatch(int x, int y) const;

    // ChemPlane in relation to the distance from a single peak F(x)
    double FoodEquation(double dist_from_peak) const;

    // Determines how much food is in a specific position
    double FoodAtPosition(int x, int y) const;

    static double DetermineMinDist(int n);

    double distMostIsolatedPoint() const;

    void updateChemPlane() const;

    int addFPatch(int x, int y);

    int addRandomFPatch();

    void removeFPatch(int id);

    int getFoodPatches() const {
        int fps = 0;
        for (auto &fp : fpatches)
            if (not fp.empty)
                ++fps;
        return fps;
    }

    int getFoodLeft() const {
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

    //! Saves information about the cells as a CSV file in the directory specified by par.celldatadir
    int saveCellData(int Time);

    int readCellData();

    void saveCellGraveData(int Time);

    void saveFoodData(int Time);

    int readFoodData();

    // Info regarding cells that have died since last time we saved data
    vector<CellGravestone> cell_graves;

private:
    //! The cells in the Petri dish; accessible to derived classes
    std::vector<Cell> cell;
    int sizex, sizey;
    // Number of resource sources
    int grad_sources;
    vector<FoodPatch> fpatches;
    static const string cell_headers;
    static const string cellgrave_headers;
};

#define INIT void Dish::Init(void)

#endif
