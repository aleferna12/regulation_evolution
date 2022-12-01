/*

Copyright 1996-2006 Roeland Merks, Paulien Hogeweg

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
// #include <iostream>
// #include <fstream>
#include <map>
#include <vector>
#include <list>
#include <set>
#include <algorithm>
#include <fstream>
#include <string.h>
#include <sstream>
#include <iterator>
#include <iostream>
#include <errno.h>
#include <math.h>
#include "dish.h"
#include "sticky.h"
#include "parameter.h"
#include "info.h"
#include "crash.h"
#include "pde.h"
#include "intplane.h"
#include "misc.h"
#include <chrono>

using namespace std::chrono;

#define EXTERNAL_OFF

extern Parameter par;

using namespace std;

Dish::Dish() {
  sizex = par.sizex;
  sizey = par.sizey;
  // This should be always 1 if we keep the gradients linear
  dist_coef = 1;

  grad_sources = par.foodpatches;

  CPM = new CellularPotts(&cell, par.sizex, par.sizey);

  FoodPlane = new IntPlane(par.sizex, par.sizey, -1);
  fpatches = vector<FoodPatch>{};

  ChemPlane = new IntPlane(par.sizex, par.sizey);

  Cell::maxsigma = 0;
  // Allocate the first "cell": this is the medium (tau=0)
  cell.push_back(*(new Cell(*this, 0)));

  // indicate that the first cell is the medium
  cell.front().sigma = 0;
  cell.front().tau = 0;

  if (par.n_chem)
    PDEfield = new PDE(par.n_chem, par.sizex, par.sizey);

  cout << "Starting the dish. Initialising..." << endl;
  // Initial cell distribution is defined by user in INIT {} block
  Init();
}


Dish::~Dish() {
  cell.clear();

  delete CPM;
  delete ChemPlane;
  delete FoodPlane;
}

// Adds a new FoodPatch at a semi-random position (still takes into account mindist)
int Dish::addRandomFPatch() {
  int x = 0, y = 0;
  if (fpatches.empty()) {
    // Ok positions = [2, side - 3] (inclusive both sides)
    x = (int) RandomNumber(sizex - 4) + 1;
    y = (int) RandomNumber(sizey - 4) + 1;
    return addFPatch(x, y);
  }
  double dist = 0;
  double mindist = DetermineMinDist(int(fpatches.size()) + 1);
  while (dist < mindist) {
    x = (int) RandomNumber(sizex - 4) + 1;
    y = (int) RandomNumber(sizey - 4) + 1;
    dist = closestFPatch(x, y).second;
  }
  return addFPatch(x, y);
}


pair<int, double> Dish::closestFPatch(int x, int y) {
  double mindist_sq = sizex * sizex + sizey * sizey;
  int res_id = -1;

  for (auto &fp : fpatches) {
    if (not fp.empty) {
      double dx = fp.getCenterX() - x;
      double dy = fp.getCenterY() - y;
      double dist_sq = dx * dx + dy * dy;
      if (dist_sq < mindist_sq) {
        mindist_sq = dist_sq;
        res_id = fp.getId();
      }
    }
  }

  return {res_id, sqrt(mindist_sq)};
}

// Alternatively, we could use the most isolated point to know where to put next peak at each iteration
// I think that doing this would be worse, as it is more computationally intensive and probably will tend to accumulate
// fpatches in the corners (?), which may be problematic for small grad_sources numbers
double Dish::DetermineMinDist(int n) {
  double ratio = (par.sizey - 2) / (double) (par.sizex - 2);
  // ratio * sepx = sepy
  // sepx = sepy / ratio
  // (sepx + 1) * (sepy + 1) = grad_sources - 1
  // ratio * pow(sepx, 2) + sepx * (1 + ratio) + 2 - grad_sources = 0
  // Do the same for sepy and solve quadradic equations
  double sepx = SolveQuadradic(ratio, 1 + ratio, 2 - n);
  double sepy = SolveQuadradic(1 / ratio, 1 + 1 / ratio, 2 - n);
  double mindistx = (par.sizex - 2) / (sepx * 2 + 2);
  double mindisty = (par.sizey - 2) / (sepy * 2 + 2);
  return sqrt(mindistx * mindistx + mindisty * mindisty);
}

void Dish::InitKeyLock() {
  auto key_lock_length = (size_t) par.key_lock_length; //take size of key and lock
  for (const auto &kl_pair: par.keylock_list) {
    if (kl_pair.tau == 0) continue;
    if (kl_pair.key.size() != key_lock_length || kl_pair.lock.size() != key_lock_length) {
      cerr << "Dish::InitKeyLock(): error. Initial key and lock vectors are not of size par.key_lock_length = "
           << par.key_lock_length << endl;
      exit(1);
    }
  }

  vector<Cell>::iterator c;
  for (c = cell.begin(), ++c; c != cell.end(); ++c) {
    int current_tau = c->getTau();
    if (current_tau == PREY) {
      //cerr<<"Hello, got prey"<<endl;
      //cerr<<"Now: "<<PREY<<endl;

      c->setJkey(par.keylock_list[PREY].key);
      c->setJlock(par.keylock_list[PREY].lock);

      //cerr<<"Check "<< par.keylock_list[ PREY ].key[0] <<endl;


    } else if (current_tau == PREDATOR) {
      //cerr<<"Hello, got predator"<<endl;
      //cerr<<"Now: "<<PREDATOR<<endl;

      c->setJkey(par.keylock_list[PREDATOR].key);
      c->setJlock(par.keylock_list[PREDATOR].lock);
      //cerr<<"Check "<< par.keylock_list[ PREDATOR ].key[0] <<endl;

    } else {
      cerr << "Dish::InitKeyLock(): error. Are there more than two types? got tau = " << c->getTau() << endl;
      exit(1);
    }

  }

  //exit(1);
}

// In this version we can let user define the lookup table for J with medium
//notice we use Paulien's method of summing powers
// Also, is this half J value? think so... - YES
int Dish::CalculateJwithMedium(vector<int> key) {
  int keypos_formedium = par.Jmedr.keypos_formedium;
  vector<int> lookup_table = par.Jmedr.lookup_table;
  int offset = par.Jmedr.offset;
  int Jval = 0;


  // Per renske suggestion, I will try a different range, more contained
  // if I just add 10 to the range (0->15)
  //I get range (10->25)
  for (int i = 0; i < keypos_formedium; i++)
    //Jval += key[i]*pow(2.0,keypos_formedium-i-1); //so that zeroth bit is most significant
    Jval += key[i] * lookup_table[i];
  Jval += offset;//so that interaction with medium can't be 0

  return (int) Jval;
}

int Dish::CalculateJfromKeyLock(vector<int> key1, vector<int> lock1, vector<int> key2, vector<int> lock2) {
  int score = 0;

  for (int i = 0; i < par.key_lock_length; i++) {
    score += (key1[i] != lock2[i]) ? 1 : 0;
    score += (key2[i] != lock1[i]) ? 1 : 0;
  }
  //now perfect score is 20, make... a sigmoid response?
  // with 20 you should get very low J val (high adhesion)
  // with 0 you should get high J val (low adh)
  //This is a arbitrary function 3+40*exp(-0.01*x^2)
  //add 0.5 before truncation to int makes it apprx to closest integer

  int Jfromkeylock = (int) (52. - 48. * ((double) score) / (2. * par.key_lock_length));
  // int Jfromkeylock = 3 + (int)(0.5+ 40.*exp( -pow( (score/double(par.key_lock_length)) , 2.) ));


  /*
  cout<<"CalculateJfromKeyLock: I got:"<<endl<< "key1: ";
  for (auto i: key1)
    cout << i << " ";
  cout<<"lok1: ";
  for (auto i: lock1)
    cout << i << " ";
  cout<<endl<<"lok2: ";
  for (auto i: lock2)
    cout << i << " ";
  cout<<"key2: ";
  for (auto i: key2)
    cout << i << " ";
  cout<<endl<<"score = "<<score<<" Jval = "<<Jfromkeylock<<endl;
  */

  return Jfromkeylock;
}

void Dish::InitVectorJ() //Initialise vector of J values for each cell
{
  std::vector<Cell>::iterator ci;
  std::vector<Cell>::iterator cj;


  int thisval;
  for (ci = cell.begin(); ci != cell.end(); ++ci) {
    for (cj = ci, ++cj; cj != cell.end(); ++cj) {

      if (ci->Sigma() == MEDIUM) {
        if (cj->Sigma() == MEDIUM) thisval = 0;
        else thisval = CalculateJwithMedium(cj->getJkey());  //J with medium is calculated from part of cell's
        ci->setVJ_singleval(cj->Sigma(), thisval);
        cj->setVJ_singleval(MEDIUM, thisval); //this cj.sigma=0, we update its J values as well

      } else {
        thisval = CalculateJfromKeyLock(ci->getJkey(), ci->getJlock(), cj->getJkey(), cj->getJlock());
        ci->setVJ_singleval(cj->Sigma(), thisval);
        cj->setVJ_singleval(ci->Sigma(), thisval);

//         if(thisval<=1){
//           cerr<<"ci="<<ci->Sigma()<<", cj->sigma="<<cj->Sigma()<<", J=" <<thisval<<endl;
//           cerr<<"keylock"<<endl;
//           for (auto i = ci->getJkey().begin(); i != ci->getJkey().end(); ++i)
//             std::cout << *i ;
//           std::cout << ' ';
//           for (auto i = ci->getJlock().begin(); i != ci->getJlock().end(); ++i)
//             std::cout << *i ;
//           std::cout<<endl;
//           for (auto i = cj->getJkey().begin(); i != cj->getJkey().end(); ++i)
//             std::cout << *i ;
//           std::cout << ' ';
//           for (auto i = cj->getJlock().begin(); i != cj->getJlock().end(); ++i)
//             std::cout << *i ;
//           std::cout<<endl;
//
//           exit(1);
//         }

      }
    }
  }

  cerr << "InitVectorJ done" << endl;
//   for(auto c: cell){
//     cerr<<"Cell "<<c.Sigma()<<", tau: "<<c.getTau()<<": " ;
//     for(auto jval: c.getVJ())
//       cerr<<jval<<" ";
//     cerr<<endl;
//   }

}

//This function and the one above could be merged
// sigma_to_update is a vector of int, there are a lot of zeros,
// but the numbers are the new sigmas (see DivideCells in ca.cpp)
// Note: some optimisation is possible here because we are updating twice the new cells (all the upd_sigma)
void Dish::UpdateVectorJ(const vector<int> &sigma_to_update) {
  //cerr<<"UpdateVectorJ begin, cell vector size: "<<cell.size()<<endl;

  vector<Cell>::iterator c;

  for (auto upd_sigma: sigma_to_update) {
    if (upd_sigma != 0) {
      for (c = cell.begin(); c != cell.end(); ++c) {
        if (0 == c->Area() && 0 == c->AliveP())
          continue; //used to be if !c->AliveP(), but some cells are dead but not disappeared yet.
        if (c->Sigma() != MEDIUM) {
          if (c->Sigma() == upd_sigma) continue;
          //update for each cell, their interactions with cell at position sigma to update,
          int jval = CalculateJfromKeyLock(c->getJkey(), c->getJlock(), cell[upd_sigma].getJkey(),
                                           cell[upd_sigma].getJlock());  //k1,l1,k2,l2
//           if(jval<=1){
//             cerr<<"upd_sigma="<<upd_sigma<<", c->sigma="<<c->Sigma()<<", J=" <<jval<<endl;
//             cerr<<"keylock"<<endl;
//             for (auto i = cell[upd_sigma].getJkey().begin(); i != cell[upd_sigma].getJkey().end(); ++i)
//               std::cout << *i ;
//             std::cout << ' ';
//             for (auto i = cell[upd_sigma].getJlock().begin(); i != cell[upd_sigma].getJlock().end(); ++i)
//               std::cout << *i ;
//             std::cout<<endl;
//             for (auto i = c->getJkey().begin(); i != c->getJkey().end(); ++i)
//               std::cout << *i ;
//             std::cout << ' ';
//             for (auto i = c->getJlock().begin(); i != c->getJlock().end(); ++i)
//               std::cout << *i ;
//             std::cout<<endl;
//
//             exit(1);
//           }
          c->setVJ_singleval(upd_sigma, jval);
          //the same number can be used by cell at which sigma is being updated
          cell[upd_sigma].setVJ_singleval(c->Sigma(), jval);
        } else {
          int jval = CalculateJwithMedium(cell[upd_sigma].getJkey()); // needs only key

          c->setVJ_singleval(upd_sigma, jval);
          cell[upd_sigma].setVJ_singleval(c->Sigma(), jval); //update nrg of medium with cell c

        }
      }
    }
  }

}


// sigma_newcells is an int vector as lognas there are cells,
// it is zero everywhere, except at the positions of a mother cell's sigma,
// where it contains as value the sigma of the daughter cell
void Dish::MutateCells(const vector<int> &sigma_to_update) {
  for (auto upd_sigma: sigma_to_update) {
    if (upd_sigma != 0) {
      cell[upd_sigma].MutateKeyAndLock();
      //assign new gextiming
      cell[upd_sigma].setGTiming((int) (RANDOM() * par.scaling_cell_to_ca_time));
      if (par.evolreg) {
        cell[upd_sigma].MutateGenome(par.mu, par.mustd);
      }
    }
  }
}

// Notice that at this stage cells are completely isolated,
// thus completely surrounded by medium
// Also, they should be far enough from borders,
//nevertheless, we are checking that
void Dish::InitContactLength() {
  int k;
  int celltype, celltypeneigh;
  int sigma, sigmaneigh;
  std::vector<Cell>::iterator c;

  for (c = cell.begin(); c != cell.end(); ++c) {
    c->clearNeighbours();
  }

  for (int x = 1; x < par.sizex - 1; x++) {
    for (int y = 1; y < par.sizey - 1; y++) {
      sigma = CPM->Sigma(x, y); //focus is on a cell - provided it is not medium
      if (sigma) {
        for (k = 1; k <= CPM->n_nb; k++)//go through neighbourhood of the pixel
        {

          int neix = x + CellularPotts::nx[k];
          int neiy = y + CellularPotts::ny[k];
          if (neix <= 0 || neix >= par.sizex - 1 || neiy <= 0 || neiy >= par.sizey - 1) {
            cerr << "InitContactLength(): warning. Neighbour is beyond borders" << endl;
            if (par.periodic_boundaries) {
              cerr << "Wrapped boundary condition applies" << endl;
              if (neix <= 0) neix = par.sizex - 2 + neix;
              if (neix >= par.sizex - 1) neix = neix - par.sizex + 2;
              if (neiy <= 0) neiy = par.sizey - 2 + neiy;
              if (neiy >= par.sizey - 1) neiy = neiy - par.sizey + 2;
            } else {
              cerr << "Fixed boundary condition applies: neighbour contact discarded." << endl;
              continue;
            }
          }

          sigmaneigh = CPM->Sigma(neix, neiy);
          if (sigmaneigh != sigma)//medium can also be a neighbour!
          {
            cell[sigma].updateNeighbourBoundary(sigmaneigh, 1);
            //cout<<"We updated cell "<<sigma<<" to have boundary with "<<sigmaneigh<<endl;
          }
        }
      }
    }
  }
  //PrintContactList();
  // PRINTING INITIALISED CONTACTS - THIS SHOULD PRINT ONLY MEDIUM -- TRUE
  //   cout<<"cell 1 has "<<cell[1].neighbours[0].first<<" contacts with cell 0"<<endl;
}

//call this function after every MCS to update the contact duration between neighbours and to erase neighbours
//with whom contact has been lost
void Dish::UpdateNeighDuration() {
  std::vector<Cell>::iterator c;
  std::map<int, pair<int, int> >::iterator n, prev;

  //cerr<<"Hello UpdateNeighDuration begin"<<endl;
  for (c = cell.begin(), c++; c != cell.end(); ++c) {
    if (!c->AliveP()) continue;
    n = c->neighbours.begin();
    //cerr<<"Hello UpdateNeighDuration 0"<<endl;
    while (n != c->neighbours.end()) {
      if (n->first == -1) {
        cerr << "We got a cell that has boundary as neighbour, n=-1:" << endl;
        exit(1);
      }
      if (n->second.first == 0) {
        prev = n;
        n++;
        //cerr<<"Hello UpdateNeighDuration 1"<<endl;
        c->neighbours.erase(prev);
        //cerr<<"Hello UpdateNeighDuration 2"<<endl;
        //cout<<"erasing contact of cell "<<c->Sigma()<<" with cell "<<n->first<<endl;
      } else {
        n->second.second++;
        n++;
      }

    }

  }
  //cerr<<"Hello UpdateNeighDuration end"<<endl;
}

// Colors for food are indicised from 10 to 60, with some simple calculations it should be easy
// to make them pretty
void Dish::ChemPlot(Graphics *g) const {
  // cpm->sigma[x][y] returns sigma, which I can use to indicise the vector of cells... can I? yes_
  int startcolorindex = 16;
  int ncolors = 29;
  auto minmaxfood = ChemPlane->getMinMax();

  // suspend=true suspends calling of DrawScene
  for (int x = 1; x < par.sizex - 1; x++)
    for (int y = 1; y < par.sizey - 1; y++)

      if (ChemPlane->Sigma(x, y) != 0) {
        if (CPM->Sigma(x, y) == 0) {
          if (ChemPlane->Sigma(x, y) < 0) {
            cerr << "chemplane below zero!!" << endl;
          }
          int colori;
          if (minmaxfood.first == minmaxfood.second) {
            colori = startcolorindex;
          } else {
            colori = startcolorindex +
                     ncolors * (ChemPlane->Sigma(x, y) - minmaxfood.first) / (minmaxfood.second - minmaxfood.first);
          }
          // Make the pixel four times as large
          // to fit with the CPM plane
          g->Point(colori, 2 * x, 2 * y);
          g->Point(colori, 2 * x + 1, 2 * y);
          g->Point(colori, 2 * x, 2 * y + 1);
          g->Point(colori, 2 * x + 1, 2 * y + 1);
        } else { ;
          // it's getting a bit cumbersome to look at this, for now I'll do without
          // g->Point(60+ChemPlane->Sigma(x,y),2*x,2*y);
          // g->Point(60+ChemPlane->Sigma(x,y),2*x+1,2*y);
          // g->Point(60+ChemPlane->Sigma(x,y),2*x,2*y+1);
          // g->Point(60+ChemPlane->Sigma(x,y),2*x+1,2*y+1);
        }
      }
}

void Dish::FoodPlot(Graphics *g, int colori) const {
  for (int x = 1; x < par.sizex - 1; x++)
    for (int y = 1; y < par.sizey - 1; y++)
      if (FoodPlane->Sigma(x, y) != -1 and CPM->Sigma(x, y) == 0) {
        g->Point(colori, 2 * x, 2 * y);
        g->Point(colori, 2 * x + 1, 2 * y);
        g->Point(colori, 2 * x, 2 * y + 1);
        g->Point(colori, 2 * x + 1, 2 * y + 1);
      }
}

void Dish::Plot(Graphics *g, int colour) {
  if (CPM) {
    CPM->Plot(g, colour);
  }

  //here chem grad plotting, with info from cpm and cell
  ChemPlot(g);
  FoodPlot(g, 100);

  //Plot direction arrows, with line function from X11?
  if (par.startmu > 0) {
    for (auto &c: cell) {
      if (c.sigma == 0 or !c.alive) continue;
      int x1 = 2 * c.meanx;
      int y1 = 2 * c.meany;
      int x2 = 2 * (c.meanx + 5 * c.tvecx);
      if (x2 >= 2 * par.sizex) x2 = 2 * par.sizex; //if too large or too small, truncate it
      else if (x2 < 0) x2 = 0;
      int y2 = 2 * (c.meany + 5 * c.tvecy);
      if (y2 >= 2 * par.sizey) y2 = 2 * par.sizey;
      else if (y2 < 0) y2 = 0;
      //now we have to wrap this
      // do we really? we could just truncate vectors up to the max size..
      g->Line(x1, y1, x2, y2, 1); //notice that Line just calls Point for drawing,
      // so it does not intrinsically suffer from x,y inversion
    }
  }
}


[[deprecated]]
int Dish::WritePeaksData() const {
  int num_rows = 5;
  ofstream file;
  file.open(par.peaksdatafile);
  for (int i = 0; i < num_rows; i++) {
    int row = i * sizex / num_rows + 1;
    for (int col = 1; col < sizey - 1; col++) {
      file << ChemPlane->Sigma(row, col) << ",";
    }
    file << endl;
  }
  file.close();
  return 0;
}


void Dish::CellsEat(int time) {
  for (auto &c: cell) {
    if (c.AliveP()) {
      if (time % par.metabperiod == 0)
        --c.food;

      if (c.getTau() == PREY) {
        int chemsumx = 0, chemsumy = 0, chemtotal = 0;
        BoundingBox bb = c.getBoundingBox();
        int pixel_count = 0;
        for (int x = bb.getMinX(); x < bb.getMaxX(); ++x) {
          for (int y = bb.getMinY(); y < bb.getMaxY(); ++y) {
            if (CPM->Sigma(x, y) == c.Sigma()) {
              ++pixel_count;
              if (time - c.last_meal > par.eatperiod) {
                int fp_id = FoodPlane->Sigma(x, y);
                if (fp_id != -1) {
                  c.food += fpatches[fp_id].consumeFood(x, y);
                  c.last_meal = time;
                }
              }
              int chem_xy = ChemPlane->Sigma(x, y);
              chemsumx += x * chem_xy;
              chemsumy += y * chem_xy;
              chemtotal += chem_xy;
            }
          }
        }
        if (pixel_count != c.Area()) {
          cerr << "Cell area is " << c.Area() << " but only " << pixel_count
               << " pixels were found inside bounding box";
          cerr << "Terminating the program";
          exit(1);
        }

        if (chemtotal) {
          double xvector = chemsumx / (double) chemtotal - c.meanx;
          double yvector = chemsumy / (double) chemtotal - c.meany;
          c.grad_conc = chemtotal / c.Area();
          double hyphyp = hypot(xvector, yvector);

          // in a homogeneous medium, gradient is zero
          // we then pick a random direction
          if (hyphyp > 0.0001) {
            xvector /= hyphyp;
            yvector /= hyphyp;
            c.setChemVec(xvector, yvector);
          } else {
            double theta = 2. * M_PI * RANDOM();
            c.setChemVec(cos(theta), sin(theta));
          }
        } else {
          double theta = 2. * M_PI * RANDOM();
          c.setChemVec(cos(theta), sin(theta));
        }

        if (c.chemvecx > 1 || c.chemvecy > 1) {
          std::cerr << ", vector: " << c.chemvecx << " " << c.chemvecy << '\n';
          exit(1);
        }
      }
    }
  }

  // TODO: This needs big changes to only update the right sections (call to removeFPatch)
  bool update_chem = false;
  for (auto &fp : fpatches) {
    if (fp.empty and not fp.removed) {
      fp.removed = true;
      update_chem = true;
    }
  }
  if (update_chem) {
    updateChemPlane();
  }
}


double Dish::distMostIsolatedPoint() {
  double dist = 0;
  for (int i = 1; i < sizex - 1; i++)
    for (int j = 1; j < sizey - 1; j++) {
      double closest_dist = closestFPatch(i, j).second;
      if (closest_dist > dist) {
        dist = closest_dist;
      }
    }
  return dist;
}


void Dish::updateChemPlane() {
  for (int i = 1; i < sizex - 1; i++)
    for (int j = 1; j < sizey - 1; j++) {
      double dfood = FoodAtPosition(i, j);
      int local_maxfood = (int) dfood;
      ChemPlane->setSigma(i, j, local_maxfood);
      if (RANDOM() < dfood - local_maxfood) local_maxfood++;
      if (RANDOM() < par.gradnoise)
        ChemPlane->setSigma(i, j, local_maxfood);
    }
  // If we want to see transversal profiles of the plane
  // WritePeaksData();
}

//to initialise cells' mu, perstime and persdur
void Dish::InitCellMigration() {
  auto icell = std::begin(cell);
  ++icell;  //discard first element of vector cell, which is medium

  //when the initialisation period has passed: start up the vectors and migration
  for (auto end = std::end(cell); icell != end; ++icell) {
    if (icell->getTau() == 1) {
      icell->setMu(par.startmu);
      icell->startTarVec();
      if (par.persduration < par.mcs) {
        icell->setPersTime(
                int(par.persduration * RANDOM())); //so that cells don't all start turning at the same time...
      } else {
        icell->setPersTime(0); //special type of experiment
        icell->tvecx = 1.;
        icell->tvecy = 0.;
      }
      icell->setPersDur(par.persduration);

      icell->setChemMu(par.init_chemmu);
      icell->startChemVec();
    } else {
      icell->setMu(0.);
      icell->setPersTime(0);
      icell->setPersDur(0);
      icell->setChemMu(0.);
    }

    //cerr<<"Cell "<<icell->sigma<<" vecx="<<icell->tvecx<<" vecy="<<icell->tvecy<<endl;
    //cerr<<"Cell persdur "<<icell->persdur<<" perstime "<<icell->perstime<<endl;
  }
  cerr << "init chemmu is" << par.init_chemmu << endl;
}

//function to have cells update their persistence time (perstime);
//In the future perhaps also their persistence duration (persdur), or how long they remember their preferred direction;
//and force of migration (mu)
void Dish::CellMigration() {
  auto icell = std::begin(cell);
  ++icell;  //discard first element of vector cell, which is medium
  for (auto end = std::end(cell); icell != end; ++icell) {
    if (!icell->AliveP()) continue; //if cell is not alive, continue

    icell->updatePersTime();

  }
}


void Dish::UpdateCellParameters(int Time) {
  vector<Cell>::iterator c; //iterator to go over all Cells
  vector<int> to_divide;
  vector<int> to_kill;
  array<double, 2> inputs = {0., 0.}; //was inputs(2,0.);
  array<int, 2> output = {0, 0};
  int interval;

  //cout<<"Update Cell parameters "<<Time<<endl;
  //update networks asynchronously f
  for (c = cell.begin(), ++c; c != cell.end(); ++c) {
    if (c->AliveP()) {
      // Mark cell to die
      if (c->food <= 0 or RANDOM() < 0.00001) {
        to_kill.push_back(c->Sigma());
        continue;
      }

      c->time_since_birth++;
      interval = Time + c->Gextiming();
      //update the network withing each cell, if it is the right time
      if (!(interval % par.scaling_cell_to_ca_time)) {
        //calculate inputs
        inputs[0] = (double) c->grad_conc;
        double division_cost = par.scaling_cell_to_ca_time * (par.divtime + par.divdur) / (double) par.metabperiod;
        inputs[1] = c->food / division_cost;
        c->UpdateGenes(inputs, true);
        c->FinishGeneUpdate();
        //what is the state of the output node of the cell?
        c->GetGeneOutput(output);

        //cell decides to divide
        if (output[0] == 1) {
          //cout<<"cell "<<c->Sigma()<<" wants to divide"<<endl;
          c->dividecounter++;

          if (c->dividecounter >= par.divtime + par.divdur) {
            //divide
            if (c->Area() > 30) {
              //cout<<"cell "<<c->Sigma()<<" will divide"<<endl;
              if (!par.nodivisions) {
                // Divide cells latter. Updating params while dividing did not work (Segmentation faults)
                to_divide.push_back(c->Sigma());
              } else {
                c->AddTimesDivided();
              }
            }
            //we already set the target area back to normal. We won't run any AmoebaeMove in between this and division
            //like this both daughter cells will inherit the normal size
            //and if the cell was too small, it needs to start all over anyway. (Hopefully a rare case)
            c->SetTargetArea(par.target_area);
            c->dividecounter = 0;
            c->ClearGenomeState(); //reset the GRN!
          }
          //not time to divide yet, but do stop migrating and start growing
          else if (c->dividecounter > par.divtime) {
            //cout<<"cell "<<c->Sigma()<<" starting to divide"<<endl;
            if (c->TargetArea() < par.target_area * 2)
              c->SetTargetArea(c->TargetArea() + 1);
            c->setMu(0.);
            c->setChemMu(0.0);
            c->setTau(2); //basically only for color right now...
          }
        }
          //this is a migratory cell
        else {
          //if (c->dividecounter) cout<<"cell "<<c->Sigma()<<" stopped division program"<<endl;
          c->dividecounter = 0;
          c->setMu(par.startmu);
          c->setChemMu(par.init_chemmu);
          c->SetTargetArea(par.target_area);
          c->setTau(1);
          //cout<<"cell "<<c->Sigma()<<" is a migratory cell"<<endl;
        }
      }

      //check area:if cell is too small (whether alive or not) we remove its sigma
      // notice that this keeps the cell in the cell array, it only removes its sigma from the field
      if (c->Area() < par.min_area_for_life) {
        to_kill.push_back(c->Sigma());
      }
    }
  }
  for (auto c_sigma : to_kill) {
    CPM->killCell(c_sigma);
  }

  vector<int> sigma_newcells;
  for (int c_sigma: to_divide) {
    int new_sigma = CPM->DivideCell(c_sigma, cell[c_sigma].getBoundingBox());
    sigma_newcells.push_back(new_sigma);
  }
  MutateCells(sigma_newcells);
  UpdateVectorJ(sigma_newcells);
}

void Dish::removeFPatch(int id) {
  // TODO
}

int Dish::addFPatch(int x, int y) {
  for (auto &fp : fpatches) {
    if (fp.removed) {
      int id = fp.getId();
      fpatches[id] = FoodPatch(this, id, x, y, fp.getLength(), fp.getFoodPerSpot());
      return id;
    }
  }
  fpatches.emplace_back(this, fpatches.size(), x, y, par.foodpatchlength, par.foodperspot);
  return fpatches.back().getId();
}

//death based on distance from the gradient until popsize is restored
void Dish::GradientBasedCellKill(int popsize) {
  std::vector<pair<int, double> > sig_dist;
  double distance, rn;
  double deathprob;
  //double fitness, totalfit, thisfit;
  //totalfit=0.;
  //int removei,removesig;


  //determine final distance of those cells that are alive
  for (auto &c: cell) {
    if (c.Sigma() > 0 && c.AliveP()) {
      distance = closestFPatch((int) round(c.getXpos()), (int) round(c.getYpos())).second;
      //fitness=1./(1.+pow(distance,3.)/par.fitscale); //for relative version
      //totalfit+=fitness; //for relative version
      //sig_dist.push_back(make_pair(c.Sigma(),totalfit)); //for relative version
      deathprob = par.mindeathprob + (par.maxdeathprob - par.mindeathprob) * pow(distance, 3.) /
                                     (pow(par.fitscale, 3.) + pow(distance, 3.));
      //cerr<<"distance is "<<distance<<", deathprob is "<<deathprob<<endl;
      sig_dist.emplace_back(c.Sigma(), deathprob);
    }
  }

  //this is a nonrelative version, based on an individually-determined death rate
  for (auto &n: sig_dist) {
    rn = RANDOM();//should be between 0 and 1

    if (rn < n.second && cell[n.first].AliveP()) {
      CPM->killCell(n.first);
    } else if (cell[n.first].AliveP()) {
      cell[n.first].ResetTimesDivided();
      cell[n.first].ClearGenomeState();
      cell[n.first].setTau(1);
    }
  }
}

int Dish::CountCells() const {

  int amount = 0;
  vector<Cell>::const_iterator i;
  for ((i = cell.begin(), ++i); i != cell.end(); ++i) {
    if (i->AliveP()) {
      // cerr<<"Time since birth: "<<i->time_since_birth<<endl;
      amount++;
    } //else {
    //cerr << "Dead cell\n";
    //}
  }
  return amount;
}

int Dish::CountCellGroups() const {

  int amount2 = 0, amount3 = 0;
  vector<Cell>::const_iterator i;
  for ((i = cell.begin(), ++i); i != cell.end(); ++i) {
    if (i->AliveP() && i->Colour() == 2) {
      amount2++;
    }
    if (i->AliveP() && i->Colour() == 3) {
      amount3++;
    }

  }
  if (amount2 == 0 || amount3 == 0) {
    return 1;
  }
  return 0;
}

int Dish::Area() const {

  int total_area = 0;

  vector<Cell>::const_iterator i;
  for ((i = cell.begin(), i++);
       i != cell.end();
       ++i) {

    total_area += i->Area();

  }
  return total_area;
}

int Dish::TargetArea() const {

  int total_area = 0;

  vector<Cell>::const_iterator i;
  for ((i = cell.begin(), i++); i != cell.end(); ++i) {

    if (i->AliveP())
      total_area += i->TargetArea();

  }
  return total_area;
}


double Dish::FoodEquation(double dist_from_peak) const {
  return par.gradscale * FoodPlane->getDiagonal() / 100 * (1 - dist_from_peak / FoodPlane->getDiagonal());
}


double Dish::FoodAtPosition(int x, int y) {
  // double pfood_j = 0.125;

  // makes gradient
  // int maxval = 3;
  // int maxval = 1+5.* (1. - dist_from_peak/(double)sizey);

  //This is how it was before, worked for field size of 500
  // double dfood = 1+5.* (1. - dist_from_peak/(double)sizey); //this the usable line
  // so maybe - to standardize gradients across field sizes, I could do:
  // dfood = 1 + sizey/100 * (1. - dist_from_peak/(double)sizey)
  // so that the local slope of the gradient stays the same?
  // also- the 1+ part of the equation could go...
  // or even better counter balanced by a lesser gradient in the variable part
  double dfood = 0;
  // TODO: Solve the interference (can it be safely implemented as parameter? Do we even want that?)
  double dist_from_peak = closestFPatch(x, y).second;
  dfood = FoodEquation(dist_from_peak);
  dfood++;
  return dfood;
}

int Dish::SaveData(int Time) {
  std::ofstream ofs;
  ofs.open(par.datafile, std::ofstream::out | std::ofstream::app);
  int pred = 0, prey = 0;

  //make file where evey so often you dump everybody - EASY!
  // save for each indiv
  // Time pred/prey key lock J with medium, neighbors
  auto icell = std::begin(cell);
  ++icell;  //discard first element of vector cell, which is medium
  for (auto end = std::end(cell); icell != end; ++icell) {
    if (!icell->AliveP()) continue; //if cell is not alive, continue

    int itau = icell->getTau();
    int isigma = icell->Sigma();
    if (itau == PREDATOR) pred++;
    else if (itau == PREY) prey++;
    else {
      cerr << "SaveData(): Error. Got cell that is neither prey, nor predator" << endl;
      exit(1);
    }
    ofs << Time << " " << isigma << " " << itau << " "; // Time now, tau of me
    ofs << icell->getXpos() << " " << icell->getYpos() << " ";
    ofs << icell->getXvec() << " " << icell->getYvec() << " ";
    ofs << icell->getChemXvec() << " " << icell->getChemYvec() << " ";
    ofs << icell->mu << " ";
    ofs << icell->chemmu << " ";
    ofs << icell->TimesDivided() << " ";
    ofs << icell->grad_conc << " ";
    // ... and all this seems to work fine. except here...
    ofs << icell->GetTimeSinceBirth() << " "; // used to be date of birth, now it's time since birth
    ofs << icell->Colour() << " "; //used during competition experiments, to delineate groups
    for (auto x: icell->getJkey()) ofs << x; //key
    ofs << " ";
    for (auto x: icell->getJlock()) ofs << x; //lock
    ofs << " ";

    //ofs << icell->getXvec()<<" "<< icell->getYvec()<<" ";
    ofs << icell->food << " ";

    // YOU SHOULD KEEP THIS AT THE LAST, because it's not constant
    for (auto i: icell->neighbours) {
      int thisj = icell->getVJ()[i.first];
      ofs << cell[i.first].getTau() << " " << thisj << " "; // get J val with neighbors
    }

    ofs << endl;
  }

  //prey=CountPreys();
  //pred=CountPredators();
  //save
  //ofs << Time << " "<< prey << " " << pred << endl;
  ofs.flush();
  ofs.close();
  return (prey + pred);
}

void Dish::SaveAncestry(int Time) {
  char fname[300];
  sprintf(fname, "%s/anc_t%010d.txt", par.networkdir, Time);

  ofstream ofs;
  ofs.open(fname, ofstream::out);
  for (auto &c: cell) {
    if (c.AliveP()) {
      ofs << c.Sigma() << " " << c.getAncestor() << endl;
      c.resetAncestor();
    }
  }
}

void Dish::SaveNetworks(int Time) {
  char fname[300];

  for (auto &c: cell) {
    if (c.AliveP()) {
      sprintf(fname, "%s/t%010d_c%04d.txt", par.networkdir, Time, c.Sigma());
      c.WriteGenomeToFile(fname);
    }
  }
}

void Dish::SaveAdheringNeighbours(int Time) {
  char fname[300];
  sprintf(fname, "%s/neigh_t%010d.txt", par.networkdir, Time);

  ofstream ofs;
  ofs.open(fname, ofstream::out);
  for (auto &c: cell) {
    if (c.AliveP()) {
      ofs << c.Sigma();
      double e_medium = c.EnergyDifference(cell[0]);
      for (auto &n: c.neighbours) {
        // Skip cell 0 (medium)
        if (n.first == 0)
          continue;
        Cell &nc = cell[n.first];
        double e_cell = c.EnergyDifference(nc) / 2.0;
        if (nc.AliveP() and e_cell <= e_medium)
          ofs << " " << nc.Sigma();
      }
      ofs << endl;
    }
  }
  ofs.close();
}

void Dish::MakeBackup(int Time) {
  std::ofstream ofs;

  //filename, c++11 strings are concatenated by summing them
  // but because we live in the middle age we are going to use sprintf;
  char filename[300];
  sprintf(filename, "%s/backup_t%09d.txt", par.backupdir, Time);

  ofs.open(filename, std::ofstream::out | std::ofstream::app);
  ofs << Time << endl;
  // Can't check using id because first few patches could be empty
  bool first = true;
  for (auto &fp : fpatches) {
    if (not fp.empty) {
      if (not first) {
        ofs << ",";
      } else {
        first = false;
      }
      ofs << fp.getX() << " " << fp.getY();
    }
  }
  ofs << endl;
  //each cell's variables are written on a single line.
  //don't store area or neighbours: can be inferred from the stored plane
  for (auto &c: cell) {
    if (c.sigma == 0) continue;
    ofs << c.sigma << " " << c.tau << " " << c.alive << " " << c.time_since_birth << " " << c.tvecx << " " << c.tvecy
        << " " << c.prevx << " " << c.prevy << " "
        << c.persdur << " " << c.perstime << " " << c.mu << " " << c.chemmu << " " << c.chemvecx << " " << c.chemvecy
        << " " << c.target_area << " "
        << c.half_div_area << " " << c.length << " " << c.last_meal << " " << c.food << " " << c.growth << " "
        << c.gextiming << " " << c.dividecounter << " " << c.grad_conc << " ";
    for (auto x: c.jkey) ofs << x; //key
    ofs << " ";
    for (auto x: c.jlock) ofs << x; //lock
    ofs << endl;
  }

  // particle plane
  // ca plane
  // posix files can be at most 2048 characters long on this laptop (LINE_MAX)
  // so a big field cannot be represented in );an obvious way,
  // however we can save some lines if we print a few sigmas in the same line
  //(anyway, field size is already specified above)... but howmany?
  // well, conceivably I will never run anything larger than in the thousands by thousands pixels
  // and so population size should be less than 10^5 <- so 5 characters+1 space,
  // so I can print 100 sigmas, which should be max 6*100=600 << 2048
  // ok!
  int counter = 0;
  int maxperline = 100;
  ofs << endl;
  for (int x = 1; x < par.sizex - 1; x++) {
    for (int y = 1; y < par.sizey - 1; y++) {
      int isigma = CPM->Sigma(x, y);
      ofs << isigma << " ";
      counter++;
      if (counter >= maxperline) {
        ofs << endl;
        counter = 0;
      }
    }
  }
  ofs << endl;
}

int Dish::ReadCompetitionFile(char *filename) {
  std::ifstream ifs;
  string line, key1, lock1, key2, lock2;
  Cell *rc;

  char grnfile1[500], grnfile2[500];

  // copying the contents of the
  // string to char array

  int placement, groupsize;

  ifs.open(filename, std::ifstream::in);

  if (ifs.is_open()) {
    //first read the placement and group sizes
    getline(ifs, line);
    stringstream strstr(line);
    strstr >> placement >> groupsize;

    strstr.clear();
    strstr.str(std::string());
    getline(ifs, line);
    strstr << line;
    strstr >> grnfile1 >> grnfile2;
    //strcpy(char_array1,grnfile1.c_str());
    //strcpy(char_array2,grnfile2.c_str());

    strstr.clear();
    strstr.str(std::string());
    getline(ifs, line);
    strstr << line;
    strstr >> key1 >> key2;

    strstr.clear();
    strstr.str(std::string());
    getline(ifs, line);
    strstr << line;
    strstr >> lock1 >> lock2;
  } else {
    return 1;
  }

  //place cells
  CPM->Place2Groups(placement, par.size_init_cells, groupsize);
  CPM->ConstructInitCells(*this);

  for (auto &c: cell) {
    if (c.Sigma()) {
      //cerr<<"Setting competition for cell "<<c.Sigma()<<endl;
      c.setGTiming((int) (RANDOM() * par.scaling_cell_to_ca_time));
      c.dividecounter = 0;
      c.SetTargetArea(par.target_area); //sets target area
      if (c.getXpos() <= par.sizex / 2) {
        c.ReadGenomeFromFile(grnfile1);
        //read the key and lock into the cell
        for (char &k: key1) {
          c.jkey.push_back(k - '0');
          //cout<<"key1 "<<k<<endl;
        }
        for (char &l: lock1) {
          c.jlock.push_back(l - '0');
          //cout<<"lock1 "<<l<<endl;
        }
        c.SetColour(2);
      } else {
        c.ReadGenomeFromFile(grnfile2);
        for (char &k: key2) {
          c.jkey.push_back(k - '0');
          //cout<<"key2 "<<k<<endl;
        }
        for (char &l: lock2) {
          c.jlock.push_back(l - '0');
          //cout<<"lock2 "<<l<<endl;
        }
        c.SetColour(3);
      }
      c.ClearGenomeState();
    }
  }
  return 0;
}


int Dish::ReadBackup(char *filename) {
  //for file reading
  std::ifstream ifs;
  string line;
  Cell *rc;
  //return value: time at which the simulation restarts
  int starttime;

  string jkey, jlock; //read them first as a string, then put into vector?
  int pos;

  ifs.open(filename, std::ifstream::in);

  if (ifs.is_open()) {
    //first read the time stamp
    getline(ifs, line);
    stringstream strstr(line);
    strstr >> starttime;

    //second read peakx and peaky
    getline(ifs, line);
    stringstream strstr2(line);

    vector<string> peakss;
    while (strstr2.good()) {
      string substr;
      getline(strstr2, substr, ',');
      peakss.push_back(substr);
    }
    for (auto &peak : peakss) {
      int x, y;
      stringstream(peak) >> x >> y;
      addFPatch(x, y);
    }

    //now read all the cell variables
    getline(ifs, line);
    while (line.length()) {
      rc = new Cell(*this); //temporary pointer to cell object
      strstr.clear();
      strstr.str(std::string());
      strstr << line;
      //read the straightforward cell variables from the line
      strstr >> rc->sigma >> rc->tau >> rc->alive >> rc->time_since_birth >> rc->tvecx >> rc->tvecy >> rc->prevx
             >> rc->prevy >> rc->persdur
             >> rc->perstime >> rc->mu >> rc->chemmu >> rc->chemvecx >> rc->chemvecy >> rc->target_area
             >> rc->half_div_area >> rc->length >>
             rc->last_meal >> rc->food >> rc->growth >> rc->gextiming >> rc->dividecounter >> rc->grad_conc
             >> jkey >> jlock;
      //read the key and lock into the cell
      for (char &c: jkey) {
        rc->jkey.push_back(c - '0');
      }
      for (char &c: jlock) {
        rc->jlock.push_back(c - '0');
      }
      rc->meanx = 0.5 * par.sizex;
      rc->meany = 0.5 * par.sizey;
      cell.push_back(*rc);

      //read the next line
      getline(ifs, line);
    }

    //now read the planes
    //hopefully, plane allocation already happened during dish creation
    //this does assume size parameters match, so be careful!
    getline(ifs, line);
    while (line.length()) {
      strstr.clear();
      strstr.str(std::string());
      strstr << line;
      while (strstr >> pos) {
        int wrong = CPM->SetNextSigma(pos);
        if (wrong) {
          cerr << "ReadBackup error in reading CA plane. more values than fit in plane." << endl;
          exit(1);
        }
      }
      getline(ifs, line);
    }

    //read the food intplane
    getline(ifs, line);
    while (line.length()) {
      strstr.clear();
      strstr.str(std::string());
      strstr << line;
      while (strstr >> pos) {
        int wrong = ChemPlane->SetNextVal(pos);
        if (wrong) {
          cerr << "ReadBackup error in reading ChemPlane plane. more values than fit in plane." << endl;
          exit(1);
        }
      }
      getline(ifs, line);
    }
  } else {
    cerr << "ReadBackup error: could not open file. exiting..." << endl;
    exit(1);
  }


  return starttime;

}

int Dish::SizeX() const { return CPM->SizeX(); }

int Dish::SizeY() const { return CPM->SizeY(); }

int Dish::Time() const {
  return CPM->Time();
}
