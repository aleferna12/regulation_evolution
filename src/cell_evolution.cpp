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

#ifndef __APPLE__

#include <malloc.h>

#endif

#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include "dish.h"
#include "random.h"
#include "cell.h"
#include "output.h"
#include "misc.h"

#ifdef QTGRAPHICS
#include "qtgraph.h"
#else

#include "x11graph.h"

#endif

//NOTE: The bookkeeping for cell contacts is very extensive:
// When cells are initially placed (call dish->InitContactLength afterwards)
// When cells divide (in cpm->dividecells)
// When cells are killed and removed (cpm->removecells)
// During CPM updates (cpm->convertspin)
// I added pieces of code to take care of this in the various applicable functions
// We may want to add a parameter to make these parts optional in case we don't need it --it's a bit more costly

using namespace std;

INIT {
    try {

        // Define initial distribution of cells
        //CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.subfield);

        // THIS IS JUST FOR EXPERIMENTS
        //CPM->PlaceOneCellsAtXY(par.sizex/2,par.sizey/2., par.size_init_cells, 1);
        //CPM->PlaceOneCellsAtXY(par.sizex/4,par.sizey/4, par.size_init_cells, 2);
        if (!strlen(par.latticefile) && !strlen(par.competitionfile)) {
            //THIS IS TO USE FOR NORMAL INITIALISATION
            if (par.scatter_start) {
                CPM->PlaceCellsRandomly(par.n_init_cells, par.size_init_cells);
            } else {
                CPM->PlaceCellsOrderly(par.n_init_cells, par.size_init_cells);
            }
            CPM->InitializeEdgeList(false);
            cout << "done initialising edge list" << endl;

            CPM->ConstructInitCells(*this); //within an object, 'this' is the object itself

            // Assign a random type to each of the cells, i.e. PREYS and PREDATORS
            CPM->SetRandomTypes();
            cout << "done setting types" << endl;
            //Initialise the contactlength bookkeeping now that the cells are placed
            // at this stage, cells are only surrounded by medium
            InitContactLength();  // see dish.cpp - you don't need dish->InitContactLength because this part IS in dish
            cout << "done setting contact length" << endl;

            cout << "Going to initialise genome" << endl;
            for (auto &c: cell) {
                if (c.Sigma()) {
                    c.setGTiming((int) (RANDOM() * par.scaling_cell_to_ca_time));
                    //c.setGTiming(0);
                    c.dividecounter = 0;
                    c.SetTargetArea(
                        par.target_area); //sets target area because in dividecells the new target area = area
                    //initialise a cell's timing for gex Updating
                    //c.setGTiming((int)(RANDOM()*par.scaling_cell_to_ca_time));
                    //creates a cell's genome, either randomly or from file
                    if (strlen(par.genomefile)) {
                        c.ReadGenomeFromFile(par.genomefile);
                    } else {
                        c.CreateRandomGenome(2, par.nr_regnodes, 1 + par.key_lock_len * 2);
                    }
                    c.ClearGenomeState();
                }
            }

            cout << "Done initialising genome" << endl;
            for (int i = 0; i < par.foodpatches; ++i)
                addRandomFPatch();
            cout << "done with food" << endl;
            updateChemPlane();
            //run CPM for some time without persistent motion
            for (int init_time = 0; init_time < 10; init_time++) {
                CPM->AmoebaeMove2(PDEfield);  //this changes neighs
            }
            cout << "done with update" << endl;
            InitCellMigration();
            UpdateCellParameters(0);//update cell status //UpdateCellParameters2();
            par.starttime = 0;
        } else if (strlen(par.competitionfile)) {
            ReadCompetitionFile(par.competitionfile);
            CPM->InitializeEdgeList(false);
            cout << "done initialising edge list for competition" << endl;
            InitContactLength();
            cout << "done initialising contacts for competition" << endl;
            // Initialises food plane (now the gradient plane)
            for (int init_time = 0; init_time < 10; init_time++) {
                CPM->AmoebaeMove2(PDEfield);  //this changes neighs
            }
            cout << "done with update" << endl;
            InitCellMigration();
            UpdateCellParameters(0);//update cell status //UpdateCellParameters2();
            par.starttime = 0;
        } else {
            cout << "Reading backfile" << endl;
            cout << "backup file is " << par.latticefile << endl;
            par.starttime = readCellData();
            if (!par.fooddatafile) {
                for (int i = 0; i < par.foodpatches; i++) {
                    addRandomFPatch();
                }
            } else if (readFoodData() != par.starttime)
                cerr << "Food data and cell data date from different times!" << endl;
            readLattice();
            CPM->InitializeEdgeList(false);
            InitContactLength();
        }
    } catch (const char *error) {
        cerr << "Caught exception\n";
        std::cerr << error << "\n";
        exit(1);
    }
}

static int last_added_fp = 0;

TIMESTEP {
    try {
        static Dish *dish = new Dish(); //here ca planes and cells are constructed
        static int i = par.starttime; //starttime is set in Dish. Not the prettiest solution, but let's hope it works.

        if (!(i % 100000)) cerr << "TIME: " << i << endl;

        //auto start = high_resolution_clock::now();
        dish->CellsEat(i);
        //auto stop = high_resolution_clock::now();
        //auto duration = duration_cast<microseconds>(stop - start);
        //sum+=duration.count();
        //nr++;
        //cout << duration.count() << endl;
        //This function updates the network and deals with the consequences of the output (motility vs division)

        dish->UpdateCellParameters(i); // for continuous GRN updating and reproduction

        dish->CellMigration();//updates persistence time and targetvectors

        dish->CPM->AmoebaeMove2(dish->PDEfield);  //this changes neighs

        dish->UpdateNeighDuration();

        if (i % 25 == 0) {
            double emptiness = 1. - dish->getFoodLeft() / (double) par.maxfood;
            if (emptiness >= 0) {
                int timer = int(par.foodpatchperiod / emptiness);
                if (i - last_added_fp > timer) {
                    last_added_fp = i;
                    dish->addRandomFPatch();
                    // TODO: Change to only update around gradient (actually do this inside addFPatch)
                    dish->updateChemPlane();
                }
            }

            if (par.evolsim) {
                if (i > par.starttime && i % par.season_duration == 0) {
                    std::cerr << "Time = " << i << '\n';
                    std::cerr << "There are " << dish->CountCells() << " cells" << '\n';

                    if (strlen(par.competitionfile)) {
                        //check if one of the groups is extinct. if yes, end simulation
                        bool extinct = dish->CountCellGroups();

                        if (extinct) {
                            std::cout << "Group extinct after " << i << " time steps. ending simulation..." << endl;
                            exit(0);
                        }
                    }
                }
            } else {
                // //not evolutionary simulation: before, used to check when enough cells passed an arbitrary boundary. don't want that now
                // if( ((strcmp(par.food_influx_location,"boundarygradient") == 0) && dish->CheckWhoMadeitLinear() ) ||
                //     ((strcmp(par.food_influx_location,"specified_experiment") == 0) && dish->CheckWhoMadeitRadial() )){
                //   //for printing switching times
                //   //write switching time to file
                //   static char timename[300];
                //   sprintf(timename,"%s/finaltime.txt",par.datadir);
                //   static ofstream myfile(timename, ios::out | ios::app);
                //   myfile << i << endl;
                //   myfile.close();
                if (i > par.starttime && i % par.season_duration == 0) {
                    exit(0);
                }
            }
        }

        // TO FILE FOR MOVIE
        if (par.store && !(i % par.storage_stride)) {
            dish->makePlots(i, this);
        }
        if (!(i % par.save_data_period)) {
            dish->saveFoodData(i);
            if (not dish->cell_graves.empty())
                dish->saveCellGraveData(i);
            int popsize = dish->saveCellData(i);
            if (not popsize) {
                cerr << "Global extinction after " << i << " time steps, simulation terminates now" << endl;
                exit(0);
            }
        }
        // TO FILE FOR BACKUP
        if (!(i % par.save_lattice_period)) {
            dish->saveLattice(i);
        }

        i++;
    } catch (const char *error) {
        cerr << "Caught exception\n";
        std::cerr << error << "\n";
        exit(1);
    }

}

int PDE::MapColour(double val) {

    return (((int) ((val / ((val) + 1.)) * 100)) % 100) + 155;
}

//////////////////////////////
// ------------------------ //
// ---       MAIN       --- //
// ------------------------ //
//////////////////////////////
int main(int argc, char *argv[]) {


    try {

#ifdef QTGRAPHICS
        //QCoreApplication a(argc, argv);
        QApplication a(argc, argv);
        QTimer g;
        //QApplication a2(argc, argv);
#endif

        par.Read(argv[1]); // Read parameters from file

        //command line arguments overwrite whatever is in the parameter file
        if (argc > 2) {
            int exit_valarg = par.ReadArguments(argc, argv);
            if (0 != exit_valarg) {
                par.PrintWelcomeStatement(); //see parameter.h
                exit(1);
            }
        }

        cerr << endl << "Warning, this version is ***NOT*** suitable for pde field!!!" << endl;
        //Depends on this: AddSiteToMoments (and Remove), FindCellDirections2, etc...
        cerr << endl << "WARNING, use wrapped boundaries if cells are A LOT smaller than sizex and sizey" << endl
             << endl;
        cerr << endl
             << "WARNING: DO NOT EVOLVE CHEMMU, or if you do, change the replication function (where it is always reset to init_chemmu)"
             << endl << endl;

        //check if directory for movies exists, create it if not, exit otherwise
        DoesDirExistsIfNotMakeit(par.moviedir);  //see output.cpp
        DoesDirExistsIfNotMakeit(par.latticedir);  //see output.cpp
        DoesDirExistsIfNotMakeit(par.celldatadir);  //see output.cpp
        DoesDirExistsIfNotMakeit(par.cellgravesdatadir);  //see output.cpp
        DoesDirExistsIfNotMakeit(par.fooddatadir);  //see output.cpp

        if (par.periodic_boundaries && par.lambda2 > 0.) {
            cerr << "main(): Error. Cannot have wrapped periodic boundaries and lambda2>0" << endl;
            cerr << "(because I cannot calculate second moment for cells crossing boundaries)" << endl;
            exit(1);
        }

        Seed(par.rseed);

        //QMainWindow mainwindow w;
#ifdef QTGRAPHICS
        cerr<<"wat is deze? "<<par.readcolortable<<endl;
        //exit(1);
        //QtGraphics g(par.sizex*2,par.sizey*2);

        QImage image(par.sizex*2,par.sizey*2, QImage::Format_ARGB32);
        QPainter painter(&image);
        QPaintDevice *device = painter.device();

        //QtGraphics g2(par.sizex*2,par.sizey*2);
        cerr<<"Hello 1"<<endl;
        //a->setMainWidget( &g );
        //a->connect(&g, SIGNAL(SimulationDone(void)), SLOT(quit(void)) );
        g.connect (&g, SIGNAL(timeout()), &a, SLOT(quit()));
        cerr<<"Hello 1.1"<<endl;

        //a2.setMainWidget( &g2 );
        //a2.connect(&g2, SIGNAL(SimulationDone(void)), SLOT(quit(void)) );

        if (par.graphics)
        {
          //g.show();
          //   // g2.show();
          cerr<<"Hello 2"<<endl;
        }
        a.exec();
        //a2.exec();
        cerr<<"Hello 3"<<endl;
#else
        cerr << "Using X11 graphics (batch mode). sizex and y are " << par.sizex << " " << par.sizey << endl;
        X11Graphics g(par.sizex, par.sizey);
        int t;

        for (t = 0; t <= par.mcs; t++) {
            //cerr<<"Time: "<<t<<endl;
            g.TimeStep();

        }
#endif

    } catch (const char *error) {
        std::cerr << error << "\n";
        exit(1);
    }
    return 0;
}
