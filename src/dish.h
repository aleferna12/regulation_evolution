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

#define PREY 1
#define PREDATOR 2

namespace ColourMode {
  enum { State,CellType,Sigma,Auxilliary };
}

class Dish {

  friend class Info;

public:
  Dish(void);

  /*! \brief Init defines the initial state of the virtual
    cell culture.

    Define Init() in your main file describing the simulation set up,
    within the block INIT { }. See for examples vessel.cpp and
    sorting.cpp.
  */
  void Init(void);

  void ConstructorBody(void);

  virtual ~Dish();

  /*! \brief Plot the Dish to graphics window g.

  Simply calls CPM->Plot.
  And also colors depnd on food.
  */
  void FoodPlot(Graphics *g);
  void Plot(Graphics *g, int colour);

  void InitKeyLock(void);
  void InitTargetArea(void);
  int CalculateJwithMedium( vector<int> key );
  int CalculateJfromKeyLock( vector<int> key1, vector<int> lock1, vector<int> key2, vector<int> lock2 );
  void InitVectorJ(void); //Initialise vector of J values for each cell
  void UpdateVectorJ(vector<int> sigma_to_update);
  void MutateCells(vector<int> sigma_to_update);

  void InitContactLength(void);
  void UpdateNeighDuration(void);
  void PrintReality(int which);
  void PrintContactList(int which=-1);
  void PrintCellParticles(void);

  int ZygoteArea(void) const;

  //! Returns the number of completed Monte Carlo Steps.
  int Time(void) const;

  //! Returns the number of cells in the dish, excluding apoptosed cells.
  int CountCells(void) const;

  /*! \brief Stretched induced cell growth and division.

  See Hogeweg (2000), Journal of Theoretical Biology.

  Find stretched cells, and increase their target area.
  Find enlarged cells, and divide them.*/
  void CellGrowthAndDivision(void);
  void CellsEat(void);
  void CellsEat2(void); //chenges cells direction vector based on where more food is

  void UpdateCellExpression(void);

  void Predate(void);

  void InitCellMigration(void);
  void CellMigration(void);
  //find cells that ate enough and let them grow; divide big cells and kill small cells
  void CellGrowthAndDivision2(void);
  void UpdateCellParameters(int Time);
  int CheckWhoMadeitLinear(void);
  int CheckWhoMadeitRadial(void);

  double FitnessFunction(int particles, double meanx, double meany);
  void ReproduceEndOfSeason(void);
  void RemoveCellsUntilPopIs(int popsize);

  void RemoveWhoDidNotMakeIt(void);
  void ReproduceWhoMadeIt(void);
  void ReproduceWhoMadeIt2(void); //with particles dependent reproduction
  void ReproduceWhoMadeIt3(void); //trying to save cells that reproduce a lot
  inline void ClearWhoMadeItSet(void){
    who_made_it.clear();
    for(auto &c:cell) c.particles=0;
  }

  //void CalculateJTable(void);

  //! \brief. Returns the summed area of all cells in the dish
  int Area(void) const;

  //! \brief Returns the summed of all cells target area in the dish
  int TargetArea(void) const;

  //return number of preys and predators
  // basically just how many are in Tau = 1 and 2
  //there's no error handling, so be careful
  int CountPreys(void);
  int CountPredators(void);

  int SaveData(int Time);
  void MakeBackup(int Time);
  int ReadBackup(char *filename);
  //! \brief Returns the horizontal size of the dish.
  int SizeX(void);

  //! \brief Returns the horizontal size of the dish.
  int SizeY(void);

  //! \brief Returns a reference to cell number "c"
  inline Cell &getCell(int c) {
    return cell[c];
  }

  PDE *PDEfield;
  IntPlane *Food;
  CellularPotts *CPM;
  std::set<int> who_made_it;
  //if the_line is 100 -> semicircle of area ~ 15000
  // in which 300 cells of area 50 can fit
  //int the_line = 80;

  //unsigned int howmany_makeit_for_nextgen = 30;
  //unsigned int popsize = 100;

  // Was used for gradient measurements, not functional now.
  void ClearGrads(void);

  void MeasureChemConcentrations(void);
protected:
  //! Assign a the cell to the current Dish
  void SetCellOwner(Cell &which_cell);

private:
  bool CellLonelyP(const Cell &c, int **neighbours) const;

protected:
  //! The cells in the Petri dish; accessible to derived classes
  std::vector<Cell> cell;
};

#define INIT void Dish::Init(void)

#endif
