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
#include <fstream>
#include <iostream>
#include <cmath>
#include "dish.h"
#include "sticky.h"
#include "crash.h"
#include "pde.h"
#include "intplane.h"
#include "misc.h"

#define EXTERNAL_OFF

using namespace std;


Dish::Dish() {
    sizex = par.sizex;
    sizey = par.sizey;

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

    for (auto &fp: fpatches) {
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


// sigma_newcells is an int vector as lognas there are cells,
// it is zero everywhere, except at the positions of a mother cell's sigma,
// where it contains as value the sigma of the daughter cell
void Dish::MutateCells(const vector<int> &sigma_to_update) {
    for (auto upd_sigma: sigma_to_update) {
        if (upd_sigma != 0) {
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


void Dish::makePlots(int Time, Graphics *g) {
    g->BeginScene();
    g->ClearImage();
    plotChemPlane(g, 8, 64);
    plotFoodPLane(g, 2);

    for (auto &plot : stringToVector<string>(par.plots, ' ')) {
        // The idea is to change these indexes if we ever need to update the colortable format
        if (plot == string("tau"))
            plotCellTau(g, 72, 72 + 8);
        else if (plot == string("food"))
            plotCellFood(g, 72, 16);
        else
            throw runtime_error("Unrecognized plot name in parameter 'plots'");
        plotCellBorders(g);
        plotCellVectors(g);

        char fname[300];
        sprintf(fname, "%s/%s%09d.png", par.moviedir, plot.c_str(), Time);
        g->Write(fname);
    }
    g->EndScene();
}


void Dish::plotChemPlane(Graphics *g, int start_index, int n_colors) const {
    auto minmaxfood = ChemPlane->getMinMax();

    // suspend=true suspends calling of DrawScene
    for (int x = 1; x < par.sizex - 1; x++)
        for (int y = 1; y < par.sizey - 1; y++)

            if (ChemPlane->Sigma(x, y) != 0) {
                if (ChemPlane->Sigma(x, y) < 0) {
                    cerr << "chemplane below zero!!" << endl;
                }
                int colori;
                if (minmaxfood.first == minmaxfood.second) {
                    colori = start_index;
                } else {
                    double perc = (ChemPlane->Sigma(x, y) - minmaxfood.first)
                                  / double(minmaxfood.second - minmaxfood.first);
                    colori = start_index + int(round((n_colors - 1) * perc));
                }
                // Make the pixel four times as large
                // to fit with the CPM plane
                g->Point(colori, x, y);
            }
}


void Dish::plotFoodPLane(Graphics *g, int color_index) const {
    for (int x = 1; x < par.sizex - 1; x++)
        for (int y = 1; y < par.sizey - 1; y++)
            if (FoodPlane->Sigma(x, y) != -1) {
                g->Point(color_index, x, y);
            }
}


void Dish::plotCellTau(Graphics *g, int div_index, int mig_index) {
    for (auto &c : cell) {
        int color_i;
        if (c.getTau() == DIVIDE)
            color_i = div_index;
        else
            color_i = mig_index;
        auto bb = c.getBoundingBox();
        for (int i = bb.getMinX(); i < bb.getMaxX(); ++i) for (int j = bb.getMinY(); j < bb.getMaxY(); ++j) {
            if (CPM->Sigma(i, j) == c.Sigma()) {
                g->Point(color_i, i, j);
            }
        }
    }
}


void Dish::plotCellFood(Graphics *g, int start_index, int n_colors) {
    for (auto &c : cell) {
        int tau_index = start_index;
        if (c.getTau() == MIGRATE)
            tau_index += n_colors / 2;
        double perc = min(1., c.food / par.foodstart);
        int color_i = tau_index + int(round((n_colors/2. - 1) * (1 - perc)));

        auto bb = c.getBoundingBox();
        for (int i = bb.getMinX(); i < bb.getMaxX(); ++i) for (int j = bb.getMinY(); j < bb.getMaxY(); ++j) {
            if (CPM->Sigma(i, j) == c.sigma) {
                g->Point(color_i, i, j);
            }
        }
    }
}


void Dish::plotCellVectors(Graphics *g) {
    //Plot direction arrows, with line function from X11?
    if (par.startmu > 0) {
        for (auto &c: cell) {
            if (c.sigma == 0 or !c.alive) continue;
            int x1 = int(c.meanx);
            int y1 = int(c.meany);
            int x2 = int(c.meanx + 5 * c.tvecx);
            if (x2 >= par.sizex) x2 = par.sizex; //if too large or too small, truncate it
            else if (x2 < 0) x2 = 0;
            int y2 = int(c.meany + 5 * c.tvecy);
            if (y2 >= par.sizey) y2 = par.sizey;
            else if (y2 < 0) y2 = 0;
            //now we have to wrap this
            // do we really? we could just truncate vectors up to the max size..
            g->Line(x1, y1, x2, y2, 1); //notice that Line just calls Point for drawing,
            // so it does not intrinsically suffer from x,y inversion
        }
    }
}


void Dish::plotCellBorders(Graphics *g) {
    for (auto &c : cell) {
        auto bb = c.getBoundingBox();
        for (int i = bb.getMinX(); i < bb.getMaxX(); ++i) for (int j = bb.getMinY(); j < bb.getMaxY(); ++j) {
            int sigma = c.Sigma();
            if (sigma != CPM->Sigma(i, j))
                continue;
            if (i < SizeX() and sigma != CPM->Sigma(i + 1, j) and not (CPM->Sigma(i, j - 1) and sigma != CPM->Sigma(i, j - 1)))
                g->Point(1, i + 1, j);
            if (j < SizeY() and sigma != CPM->Sigma(i, j + 1))
                g->Point(1, i, j + 1);
        }
    }
}


void Dish::CellsEat(int time) {
    for (auto &c: cell) {
        if (not c.AliveP())
            continue;

        if (time % par.metabperiod == 0)
            --c.food;

        int chemtotal = 0, chemsumx = 0, chemsumy = 0;
        auto bb = c.getBoundingBox();
        int pixel_count = 0;
        for (int x = bb.getMinX(); x < bb.getMaxX(); ++x) {
            for (int y = bb.getMinY(); y < bb.getMaxY(); ++y) {
                if (c.Sigma() != CPM->Sigma(x, y))
                    continue;

                if (time - c.last_meal > par.eatperiod) {
                    int fp_id = FoodPlane->Sigma(x, y);
                    if (fp_id != -1) {
                        c.food += fpatches[fp_id].consumeFood(x, y);
                        c.last_meal = time;
                    }
                }

                if (c.getTau() == MIGRATE) {
                    int chem_xy = ChemPlane->Sigma(x, y);
                    chemsumx += x * chem_xy;
                    chemsumy += y * chem_xy;
                    chemtotal += chem_xy;
                }

                ++pixel_count;
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

    // TODO: This needs big changes to only update the right sections (call to removeFPatch)
    bool update_chem = false;
    for (auto &fp: fpatches) {
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
    int interval;

    //cout<<"Update Cell parameters "<<Time<<endl;
    for (c = cell.begin(), ++c; c != cell.end(); ++c) {
        if (c->AliveP()) {
            // Death checks
            int death_period = 25;
            string death_reason;
            // notice that this keeps the cell in the cell array, it only removes its sigma from the field
            if (c->Area() < par.min_area_for_life) {
                death_reason = "squeezed";
            } else if (c->food <= 0) {
                death_reason = "starved";
            }
            // Mark cell to die
            // Only calculate prob every 25 mcs
            // TODO: Test if its still working after making asynchronous (maybe plot average lifespan or something)
            else if (c->time_since_birth % death_period == 0 and RANDOM() < par.gompertz_alpha * pow(M_E, par.gompertz_beta * c->time_since_birth / death_period)) {
                death_reason = "old";
            }
            if (not death_reason.empty()) {
                to_kill.push_back(c->Sigma());
                CellGravestone cg = {
                    c->Sigma(),
                    c->getTau(),
                    c->time_since_birth,
                    Time,
                    CPM->calculateGamma(c->Sigma(), c->Sigma()),
                    death_reason
                };
                cell_graves.push_back(std::move(cg));
                continue;
            }

            c->time_since_birth++;
            interval = Time + c->Gextiming();
            //update the network withing each cell, if it is the right time
            if (!(interval % par.scaling_cell_to_ca_time)) {
                //calculate inputs
                inputs[0] = (double) c->grad_conc;
                double division_cost =
                    par.scaling_cell_to_ca_time * (par.divtime + par.divdur) / (double) par.metabperiod;
                inputs[1] = c->food / division_cost;
                c->UpdateGenes(inputs, true);
                c->FinishGeneAndJDecsUpdate();

                //cell decides to divide
                if (c->genome.outputnodes[0].Boolstate == 1) {
                    //cout<<"cell "<<c->Sigma()<<" wants to divide"<<endl;
                    c->dividecounter++;

                    if (c->dividecounter >= par.divtime + par.divdur) {
                        //divide
                        if (c->Area() > 30) {
                            //cout<<"cell "<<c->Sigma()<<" will divide"<<endl;
                            if (!par.nodivisions) {
                                // Divide cells later. Updating params while dividing did not work (Segmentation faults)
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
                        // TODO: Shouldn't this be out of this condition?
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
        }
    }
    for (auto c_sigma: to_kill) {
        CPM->killCell(c_sigma);
    }

    vector<int> sigma_newcells;
    for (int c_sigma: to_divide) {
        int new_sigma = CPM->DivideCell(c_sigma, cell[c_sigma].getBoundingBox());
        sigma_newcells.push_back(new_sigma);
    }
    MutateCells(sigma_newcells);
}

void Dish::removeFPatch(int id) {
    // TODO
}

int Dish::addFPatch(int x, int y) {
    for (auto &fp: fpatches) {
        if (fp.removed) {
            int id = fp.getId();
            fpatches[id] = FoodPatch(this, id, x, y, fp.getLength(), fp.getFoodPerSpot(), nullptr);
            return id;
        }
    }
    fpatches.emplace_back(this, fpatches.size(), x, y, par.foodpatchlength, par.foodperspot);
    return fpatches.back().getId();
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


int Dish::readFoodData() {
    int cur_time;
    ifstream file(par.fooddatafile);
    if (not file)
        throw runtime_error("Failed to open file");

    string line;
    // Skip headers
    getline(file, line);
    int id = 0;
    while (getline(file, line)) {
        auto attrs = stringToVector<string>(line, ',');
        auto sigmas = stringToVector<int>(attrs[3], ' ');
        fpatches.emplace_back(
            this,
            id,
            stoi(attrs[0]),
            stoi(attrs[1]),
            stoi(attrs[2]),
            par.foodperspot,
            &sigmas[0]
        );
        fpatches.back().updateFoodLeft();
        cur_time = stoi(attrs[4]);
    }
    return cur_time;
}


void Dish::saveFoodData(int Time) {
    char filename[300];
    // TODO add parameter
    sprintf(filename, "%s/t%09d.csv", par.fooddatadir, Time);
    ofstream file(filename);
    if (not file)
        throw runtime_error("Failed to open file");

    vector<string> col_names {"x", "y", "length", "sigma_list", "time"};
    file << vectorToString(col_names, ',') << endl;

    for (auto &fp : fpatches) {
        if (fp.empty)
            continue;

        file << fp.getX() << ',' << fp.getY() << ',' << fp.getLength() << ',';
        file << vectorToString(fp.getSigmasAsVector(), ' ') << ',';
        file << Time << endl;
    }
}


void Dish::readLattice() {
    if (cell.size() <= 1) {
        cerr << "readLattice must be called after readData";
        exit(1);
    }

    ifstream file(par.latticefile);
    if (not file)
        throw runtime_error("Failed to open file");
    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        while (ss.good()) {
            string val;
            getline(ss, val, ',');
            CPM->SetNextSigma(stoi(val));
        }
    }
    updateChemPlane();
}

void Dish::saveLattice(int Time) const {
    char filename[300];
    sprintf(filename, "%s/t%09d.csv", par.latticedir, Time);

    ofstream file(filename);
    if (not file)
        throw runtime_error("Failed to open file");

    for (int x = 1; x < par.sizex - 1; x++) {
        for (int y = 1; y < par.sizey - 1; y++) {
            int isigma = CPM->Sigma(x, y);
            // Skips dead cells (necessary because we are not saving info about them so cant interpret these sigmas)
            if (not cell[isigma].AliveP()) {
                isigma = 0;
            }
            file << isigma;
            if (y < par.sizey - 2)
                file << ",";
        }
        file << endl;
    }
}


int Dish::readCellData() {
    int cur_time;
    ifstream file(par.celldatafile);
    if (not file)
        throw runtime_error("Failed to open file");

    string line;
    // Skip headers
    getline(file, line);
    int last_sigma = 0;
    while (getline(file, line)) {
        auto attrs = stringToVector<string>(line, ',');
        auto it = attrs.begin();

        int sigma = stoi(*it); ++it;
        for (int j = 1; j < sigma - last_sigma; ++j) {
            Cell *rc = new Cell(*this);
            rc->alive = false;
            rc->sigma = last_sigma + j;
            cell.push_back(*rc);
        }
        last_sigma = sigma;

        Cell *rc = new Cell(*this);
        rc->alive = true;
        rc->sigma = sigma;
        rc->tau = stoi(*it); ++it;
        rc->time_since_birth = stoi(*it); ++it;
        rc->tvecx = stod(*it); ++it;
        rc->tvecy = stod(*it); ++it;
        rc->prevx = stod(*it); ++it;
        rc->prevy = stod(*it); ++it;
        rc->persdur = stoi(*it); ++it;
        rc->perstime = stoi(*it); ++it;
        rc->mu = stod(*it); ++it;
        rc->half_div_area = stoi(*it); ++it;
        rc->length = stod(*it); ++it;
        rc->last_meal = stoi(*it); ++it;
        rc->food = stod(*it); ++it;
        rc->growth = stod(*it); ++it;
        rc->gextiming = stoi(*it); ++it;
        rc->dividecounter = stoi(*it); ++it;
        rc->grad_conc = stoi(*it); ++it;
        rc->meanx = stod(*it); ++it;
        rc->meany = stod(*it); ++it;
        rc->chemvecx = stod(*it); ++it;
        rc->chemvecy = stod(*it); ++it;
        rc->target_area = stoi(*it); ++it;
        rc->chemmu = stod(*it); ++it;
        rc->times_divided = stoi(*it); ++it;
        rc->colour = stoi(*it); ++it;
        rc->ancestor = stoi(*it); ++it;
        // Skip Jmed
        ++it;
        rc->jkey_dec = stoi(*it); ++it;
        rc->jlock_dec = stoi(*it); ++it;
        // Skip neighbour info
        it += 2;

        int innr = stoi(*it); ++it;
        int regnr = stoi(*it); ++it;
        int outnr = stoi(*it); ++it;
        rc->genome.innr = innr;
        rc->genome.regnr = regnr;
        rc->genome.outnr = outnr;
        rc->genome.inputscale = stringToVector<double>(*it, ' '); ++it;

        vector<double> reg_thres = stringToVector<double>(*it, ' '); ++it;
        vector<double> reg_w_in = stringToVector<double>(*it, ' '); ++it;
        vector<double> reg_w_reg = stringToVector<double>(*it, ' '); ++it;
        for (int i = 0; i < regnr; ++i) {
            Gene *gene = new Gene(1, i + innr, innr, regnr);
            gene->threshold = reg_thres[i];
            for (int j = 0; j < innr; ++j)
                gene->w_innode[j] = reg_w_in[i * innr + j];
            for (int j = 0; j < regnr; ++j)
                gene->w_regnode[j] = reg_w_reg[i * regnr + j];
            rc->genome.regnodes.push_back(*gene);
        }

        vector<double> out_thres = stringToVector<double>(*it, ' '); ++it;
        vector<double> out_w_reg = stringToVector<double>(*it, ' '); ++it;
        for (int i = 0; i < outnr; ++i) {
            Gene *gene = new Gene(2, i + innr + regnr, innr, regnr);
            gene->threshold = out_thres[i];
            for (int j = 0; j < regnr; ++j)
                gene->w_regnode[j] = out_w_reg[i * regnr + j];
            rc->genome.outputnodes.push_back(*gene);
        }

        cur_time = stoi(*it);
        cell.push_back(*rc);
    }
    return cur_time;
}


int Dish::saveCellData(int Time) {
    int n_cells = 0;

    char filename[300];
    sprintf(filename, "%s/t%09d.csv", par.celldatadir, Time);
    ofstream file(filename);
    if (not file)
        throw runtime_error("Failed to open file");

    vector<string> col_names{"sigma", "tau", "time_since_birth", "tvecx", "tvecy", "prevx", "prevy", "persdur",
                             "perstime", "mu", "half_div_area", "length", "last_meal", "food", "growth", "gextiming",
                             "dividecounter", "grad_conc", "meanx", "meany", "chemvecx", "chemvecy", "target_area",
                             "chemmu", "times_divided", "colour", "ancestor", "Jmed", "jkey_dec", "jlock_dec",
                             "neighbour_list", "Jneighbour_list", "innr", "regnr", "outnr", "in_scale_list",
                             "reg_threshold_list", "reg_w_innode_list", "reg_w_regnode_list", "out_threshold_list",
                             "out_w_regnode_list", "time"};
    file << vectorToString(col_names, ',') << endl;

    for (auto &c: cell) {
        if (not c.AliveP() or c.sigma == 0)
            continue;
        ++n_cells;
        // Cant use the helper function because the values are of different types
        file << c.sigma << ',' << c.tau << ',' << c.time_since_birth << ',' << c.tvecx << ',' << c.tvecy << ','
             << c.prevx << ',' << c.prevy << ',' << c.persdur << ',' << c.perstime << ',' << c.mu << ','
             << c.half_div_area << ',' << c.length << ',' << c.last_meal << ',' << c.food << ',' << c.growth << ','
             << c.gextiming << ',' << c.dividecounter << ',' << c.grad_conc << ',' << c.meanx << ',' << c.meany << ','
             << c.chemvecx << ',' << c.chemvecy << ',' << c.target_area << ',' << c.chemmu << ','
             << c.times_divided << ',' << c.colour << ',' << c.ancestor << ',' << par.Jmed << ',' << c.jkey_dec << ','
             << c.jlock_dec <<  ',';
        // Need to reset everytime we save data
        c.resetAncestor();

        vector<int> neighs {};
        vector<double> Jneighs {};
        for (auto &n : c.neighbours) {
            if (cell[n.first].AliveP()) {
                neighs.push_back(n.first);
                Jneighs.push_back(CPM->energyDifference(c.sigma, n.first));
            }
        }
        file << vectorToString(neighs, ' ') << ',';
        file << vectorToString(Jneighs, ' ') << ',';

        file << c.genome.innr << ',' << c.genome.regnr << ',' << c.genome.outnr << ',';
        file << vectorToString(c.genome.inputscale, ' ') << ',';

        vector<double> reg_thres {};
        vector<double> reg_w_in {};
        vector<double> reg_w_reg {};
        for (auto &g: c.genome.regnodes) {
            reg_thres.push_back(g.threshold);
            // Unpack the weight vectors
            reg_w_in.insert(reg_w_in.end(), g.w_innode.begin(), g.w_innode.end());
            reg_w_reg.insert(reg_w_reg.end(), g.w_regnode.begin(), g.w_regnode.end());
        }
        file << vectorToString(reg_thres, ' ') << ',';
        file << vectorToString(reg_w_in, ' ') << ',';
        file << vectorToString(reg_w_reg, ' ') << ',';

        vector<double> out_thres {};
        vector<double> out_w_reg {};
        for (auto &g: c.genome.outputnodes) {
            out_thres.push_back(g.threshold);
            // Unpack the weight vectors
            out_w_reg.insert(out_w_reg.end(), g.w_regnode.begin(), g.w_regnode.end());
        }
        file << vectorToString(out_thres, ' ') << ',';
        file << vectorToString(out_w_reg, ' ') << ',';

        file << Time << endl;
    }
    return n_cells;
}


void Dish::saveCellGraveData(int Time) {
    char filename[300];
    sprintf(filename, "%s/t%09d.csv", par.cellgravesdatadir, Time);
    ofstream file(filename);
    if (not file)
        throw runtime_error("Failed to open file");

    vector<string> col_names{"sigma", "tau", "age", "time_death", "reason", "self_gamma"};
    file << vectorToString(col_names, ',') << endl;

    for (auto &cg : cell_graves) {
        file << cg.sigma << ',' << cg.tau << ',' << cg.age << ',' << cg.time_death << ',' << cg.reason << ','
        << cg.self_gamma << endl;
    }
    // Cant forget to clear the vector!
    cell_graves.clear();
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
            if (c.getXpos() <= par.sizex / 2.) {
                c.ReadGenomeFromFile(grnfile1);
                //read the key and lock into the cell
                c.SetColour(2);
            } else {
                c.ReadGenomeFromFile(grnfile2);
                c.SetColour(3);
            }
            c.ClearGenomeState();
        }
    }
    return 0;
}

int Dish::SizeX() const { return CPM->SizeX(); }

int Dish::SizeY() const { return CPM->SizeY(); }

int Dish::Time() const {
    return CPM->Time();
}
