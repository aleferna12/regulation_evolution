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


#include "parameter.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cerrno>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include "output.h"
#include "parse.h"


//parameter constructor - initialise default
// these are the values things takes when par file or command line does not specify them
Parameter::Parameter() {

    T = 50.;
    target_area = 50;
    half_div_area = 50;
    half_div_area_2 = -1;
    target_length = 60;
    lambda = 50;
    lambda2 = 5.0;
    //Jtable = strdup("J.dat");
    conn_diss = 2000;
    vecadherinknockout = false;
    extensiononly = false;
    chemotaxis = 1000;
    border_energy = 100;
    neighbours = 2;
    min_area_for_life = 5;
    Jmed = 14;
    Jalpha = 7;
    key_lock_len = 6;
    n_chem = 1;
    diff_coeff = new double[1];
    diff_coeff[0] = 1e-13;
    decay_rate = new double[1];
    decay_rate[0] = 1.8e-4;
    secr_rate = new double[1];
    secr_rate[0] = 1.8e-4;
    saturation = 0;
    dt = 2.0;
    dx = 2.0e-6;
    pde_its = 15;
    n_init_cells = 100;
    size_init_cells = 25;
    sizex = 200;
    sizey = 200;
    mcs = 10000;
    rseed = -1;
    subfield = 1.0;
    relaxation = 0;
    storage_stride = 10;
    graphics = true;
    store = false;
    divisioncolour = false;
    genomefile = strdup("");
    nr_regnodes = 1;
    mu = 0.05;
    mustd = 0.1;
    divtime = 50;
    divdur = 10;
    gompertz_alpha = 0.0000075;
    gompertz_beta = 0.0004;
    scatter_start = true;
    chemgrad = true;
    groupextinction = false;
    evolvable_adh = true;
    chemcircles = true;
    circle_dist = 40;
    circle_segments = 15;
    circle_thickness = 1;
    motiledeath = 1.0;
    dividingdeath = 0.;
    moviedir = strdup("data_film");
    celldatafile = strdup("");
    fooddatafile = strdup("");
    save_data_period = 100;
    food_influx_location = strdup("nowhere");
    metabperiod = 20;
    growth = 0;
    ardecay = 0.;
    gradnoise = 0.1;
    gradscale = 1.0;
    foodpatches = 1;
    foodpatchperiod = 1000;
    foodpatchlength = 1;
    foodperspot = 1;
    maxfoodpatches = 1;
    foodstart = 1000;
    eatperiod = 1000;
    nodivisions = false;
    min_contact_duration_for_preying = 10;
    frac_contlen_eaten = 1.;
    metabolic_conversion = 0.5;
    colortable_filename = "default.ctb";
    plots = "tau";
    key_lock_weights = "1 2 3 4 5 6";
    mut_rate = 0.01;
    circle_coverage = 0.5;
    startmu = 0.0;
    persduration = 0;
    scaling_cell_to_ca_time = 1;
    latticedir = strdup("lattice");
    celldatadir = strdup("celldata");
    cellgravesdatadir = strdup("cellgravedata");
    fooddatadir = strdup("fooddata");
    save_lattice_period = 0;
    init_chemmu = 0.;
    latticefile = strdup("");
    starttime = 0;

    howmany_makeit_for_nextgen = 30;
    popsize = 100;
    the_line = 50;
    evolsim = 0;
    is_there_food = false;
    zero_persistence_past_theline = false;
    season_experiment = true;
    season_duration = 100000;
    init_cell_config = 0;
    cell_placement = 0;
}

Parameter::~Parameter() {

    // destruct parameter object

    // free string parameter

    CleanUp();

}

void Parameter::CleanUp() const {
    if (diff_coeff)
        free(diff_coeff);
    if (decay_rate)
        free(decay_rate);
    if (secr_rate)
        free(secr_rate);
    if (moviedir)
        free(moviedir);
    if (genomefile)
        free(genomefile);
    if (latticedir)
        free(latticedir);
    if (celldatadir)
        free(celldatadir);
    if (cellgravesdatadir)
        free(cellgravesdatadir);
    if (fooddatadir)
        free(fooddatadir);
    if (latticefile)
        free(latticefile);

}

void Parameter::PrintWelcomeStatement() {
    cerr << "CellEvol: v0.something (very much a prototype)" << endl;
    cerr << "Usage is: " << endl;
    cerr << "./cell_evolution path/to/data [optional arguments]" << endl;
    cerr << "Arguments: " << endl;
    cerr
            << " -name path/to/name_for_all_output # gives a name to all output, alternative to -celldatafile -fooddatafile -moviedir -latticedir -celldatadir -cellgravedatadir -fooddatadir"
            << endl;
    cerr << " -celldatafile path/to/celldatafile # output file" << endl;
    cerr << " -fooddatafile path/to/fooddatafile # output file" << endl;
    cerr << " -moviedir path/to/moviedir # output movie dir" << endl;
    cerr << " -latticedir path/to/latticedir # output backup lattice dir" << endl;
    cerr << " -celldatadir path/to/celldatadir # output cell data dir" << endl;
    cerr << " -cellgravedatadir path/to/cellgravedatadir # output cell data dir" << endl;
    cerr << " -fooddatadir path/to/fooddatadir # output food data dir" << endl;
    cerr << " -networkdir path/to/networkdir # output network dir" << endl;
    cerr << " -store # store pictures" << endl;
    cerr << " -keylockfilename path/to/keylockfilename" << endl;
    cerr << " -seed INT_NUMBER # for random number generator" << endl;
    cerr << " -maxtime INT_NUMBER" << endl;
    // cerr<<" -halfdiv_area_predator INT_NUMBER"<<endl;
    cerr << " -persmu FLOAT_NUMBER [ > 0 ], strength of persistent random walk" << endl;
    cerr << " -persduration INT_NUMBER" << endl;
    cerr << " -mutrate FLOAT_NUMBER [0,1) # mutation rate for key and lock" << endl;
    cerr << " -mu FLOAT_NUMBER [0,1) # mutation rate for genome regulation" << endl;
    cerr << " -mustd FLOAT_NUMBER [0,1) # mutation size for genome regulation" << endl;
    cerr << " -casize INT_NUMBER INT_NUMBER # dimensions of the CA" << endl;
    cerr << " -popsize INT_NUMBER [-pop_as_initpop] # population size, optional same as n_init_cells" << endl;
    cerr << " -pop_as_initpop # initial population size = popsize" << endl;
    cerr << " -n_nextgen INT_NUMBER # number of cells that are taken to next generation" << endl;
    cerr
            << " -noevolsim # No evolution at all: sim ends in 1 season, when [howmany_makeit_for_nextgen] cells pass [the_line])"
            << endl;
    cerr << " -nofood # No food distributed in the simulation" << endl;
    cerr << " -noevolreg # No evolution of regulation parameters" << endl;
    cerr << " -scatter # spread cells after a season" << endl;
    cerr << " -nodivisions # do not execute divisions -> number of cells remains the same" << endl;
    cerr << " -latticefile path/to/latticefile # to start simulation from backup" << endl;
    cerr << " -season [INT_NUMBER] # season duration" << endl;
    cerr << " -metabperiod [INT_NUMBER] how often we deduce 1 food of each cell in MCS" << endl;
    cerr << " -gradscale [FLOAT_NUMBER] slope of the gradient (in percent units)" << endl;
    cerr << " -foodpatches [INT_NUMBER] initial number of food patches resources placed in the field" << endl;
    cerr << " -foodpatchperiod [INT_NUMBER] new food patch timer (a new patch will be created every X MCS)" << endl;
    cerr
            << " -foodpatchlength [INT_NUMBER] side of each food patch (total amount of spots per patch will be this squared)"
            << endl;
    cerr << " -foodperspot [INT_NUMBER] how much food each spot contains" << endl;
    cerr << " -maxfoodpatches [INT_NUMBER] maximum number of food patches allowed to coexist in the system" << endl;
    cerr << " -foodstart [INT_NUMBER] initial food for cells" << endl;
    cerr << " -eatperiod [INT_NUMBER] how often a cell can eat" << endl;
    cerr << " -gradnoise [FLOAT_NUMBER] chances that any grid point has gradient, rather than being empty" << endl;
    cerr << " -chemmu [FLOAT_NUMBER] scaling factor for chemotaxis in the Hamiltonian" << endl;
    cerr << " -fitscale [FLOAT_NUMBER] point in field where deathrate is half value" << endl;
    cerr << " -genomefile [string] starting genome with which to seed the field" << endl;
    cerr << " -target_area [INT_NUMBER] that (initial) target area of cells" << endl;
    cerr << " -init_cell_config [0-3] initial configuration of cells when placed in center, see ca.cpp" << endl;
    cerr << " -cell_placement [1-4] field position of cells, (0=center) see ca.cpp" << endl;
    cerr << " -gompertz_alpha [FLOAT_NUMBER] alpha parameter for the hazard function a*e^(b*x)" << endl;
    cerr << " -gompertz_beta [FLOAT_NUMBER] beta parameter for the hazard function a*e^(b*x)" << endl;
    cerr << endl << "Will not execute if celldatafile and moviedir already exist" << endl;
    cerr << "Also, parameter file and Jtable should be in the same directory (unless you used option -keylockfilename)"
         << endl;
    cerr << "Have fun!" << endl;
}

int Parameter::ReadArguments(int argc, char *argv[]) {
    cerr << endl << "Reading arguments from command line" << endl;
    //starts from 2 because 0 is filename, 1 is parameter file path
    for (int i = 2; i < argc; i++) {
        if (0 == strcmp(argv[i], "-celldatafile")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in celldatafile?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            //strcpy(celldatafile, argv[i]); //this can be buggy because it copies over a dynamically allocated char* (celldatafile) that can be a lot shorter
            free(celldatafile);
            celldatafile = (char *) malloc(
                    5 + strlen(argv[i]) * sizeof(char)); //strlen(argv[i]) is ok because argv[i] is null terminated
            celldatafile = strdup(argv[i]);
            cerr << "New value for celldatafile: " << celldatafile << endl;
//       exit(1);
        } else if (0 == strcmp(argv[i], "-fooddatafile")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in fooddatafile?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            //strcpy(fooddatafile, argv[i]); //this can be buggy because it copies over a dynamically allocated char* (celldatafile) that can be a lot shorter
            free(fooddatafile);
            fooddatafile = (char *) malloc(
                5 + strlen(argv[i]) * sizeof(char)); //strlen(argv[i]) is ok because argv[i] is null terminated
            fooddatafile = strdup(argv[i]);
            cerr << "New value for fooddatafile: " << fooddatafile << endl;
        } else if (0 == strcmp(argv[i], "-moviedir")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in moviedir?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            //strcpy(moviedir, argv[i]);
            free(moviedir);
            moviedir = (char *) malloc(
                    5 + strlen(argv[i]) * sizeof(char)); //strlen(argv[i]) is ok because argv[i] is null terminated
            moviedir = strdup(argv[i]);

            cerr << "New value for moviedir: " << moviedir << endl;

        } else if (0 == strcmp(argv[i], "-latticedir")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in latticedir?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            //strcpy(moviedir, argv[i]);
            free(latticedir);
            latticedir = (char *) malloc(
                5 + strlen(argv[i]) * sizeof(char)); //strlen(argv[i]) is ok because argv[i] is null terminated
            latticedir = strdup(argv[i]);

            cerr << "New value for latticedir: " << latticedir << endl;

        } else if (0 == strcmp(argv[i], "-celldatadir")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in celldatadir?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            //strcpy(moviedir, argv[i]);
            free(celldatadir);
            celldatadir = (char *) malloc(
                5 + strlen(argv[i]) * sizeof(char)); //strlen(argv[i]) is ok because argv[i] is null terminated
            celldatadir = strdup(argv[i]);

            cerr << "New value for celldatadir: " << celldatadir << endl;

        } else if (0 == strcmp(argv[i], "-cellgravesdatadir")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in cellgravesdatadir?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            //strcpy(moviedir, argv[i]);
            free(cellgravesdatadir);
            cellgravesdatadir = (char *) malloc(
                5 + strlen(argv[i]) * sizeof(char)); //strlen(argv[i]) is ok because argv[i] is null terminated
            cellgravesdatadir = strdup(argv[i]);

            cerr << "New value for cellgravesdatadir: " << cellgravesdatadir << endl;

        } else if (0 == strcmp(argv[i], "-fooddatadir")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in fooddatadir?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            //strcpy(moviedir, argv[i]);
            free(fooddatadir);
            fooddatadir = (char *) malloc(
                5 + strlen(argv[i]) * sizeof(char)); //strlen(argv[i]) is ok because argv[i] is null terminated
            fooddatadir = strdup(argv[i]);

            cerr << "New value for fooddatadir: " << fooddatadir << endl;
        } else if (0 == strcmp(argv[i], "-networkdir")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in networkdir?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            //strcpy(moviedir, argv[i]);


        } else if (0 == strcmp(argv[i], "-genomefile")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in genomefile?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            //strcpy(moviedir, argv[i]);
            free(genomefile);
            genomefile = (char *) malloc(
                    5 + strlen(argv[i]) * sizeof(char)); //strlen(argv[i]) is ok because argv[i] is null terminated
            genomefile = strdup(argv[i]);

            cerr << "New value for genomefile: " << genomefile << endl;

        } else if (0 == strcmp(argv[i], "-seed")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in seed?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            rseed = atoi(argv[i]);
            cerr << "New value for seed: " << rseed << endl;
        } else if (0 == strcmp(argv[i], "-maxtime")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in maxtime?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            mcs = atoi(argv[i]);
            cerr << "New value for maxtime (mcs in the code): " << mcs << endl;
        } else if (0 == strcmp(argv[i], "-circle_segments")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in circle_segments?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            circle_segments = atoi(argv[i]);
            cerr << "New value for circle_segments (mcs in the code): " << circle_segments << endl;
        } else if (0 == strcmp(argv[i], "-circle_thickness")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in circle_thickness?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            circle_thickness = atoi(argv[i]);
            cerr << "New value for circle_thickness (mcs in the code): " << circle_thickness << endl;
        } else if (0 == strcmp(argv[i], "-circle_dist")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in circle_dist?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            circle_dist = atoi(argv[i]);
            cerr << "New value for circle_dist (mcs in the code): " << circle_dist << endl;
        } else if (0 == strcmp(argv[i], "-mutrate")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in mutrate?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            mut_rate = atof(argv[i]);
            cerr << "New value for mutation rate: " << mut_rate << endl;
        } else if (0 == strcmp(argv[i], "-circle_thickness")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in circle_thickness?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            circle_thickness = atof(argv[i]);
            cerr << "New value for circle_thickness: " << circle_thickness << endl;
        } else if (0 == strcmp(argv[i], "-mu")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in mu?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            mu = atof(argv[i]);
            cerr << "New value for genome mutation rate: " << mu << endl;
        } else if (0 == strcmp(argv[i], "-mustd")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in mustd?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            mustd = atof(argv[i]);
            cerr << "New value for mustd: " << mustd << endl;
        } else if (0 == strcmp(argv[i], "-persduration")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in persduration?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            persduration = atoi(argv[i]);
            cerr << "New value for persistence of movement: " << persduration << endl;
        } else if (0 == strcmp(argv[i], "-latticefile")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in latticefile?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            free(latticefile);
            latticefile = (char *) malloc(
                    5 + strlen(argv[i]) * sizeof(char)); //strlen(argv[i]) is ok because argv[i] is null terminated
            latticefile = strdup(argv[i]);

            cerr << "New value for latticefile: " << latticefile << endl;
        } else if (0 == strcmp(argv[i], "-popsize")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in popsize?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            popsize = atoi(argv[i]);
            cerr << "New value for population size: " << popsize << endl;
            i++;
            if (i == argc) return 0;
            if (0 == strcmp(argv[i], "-pop_as_initpop")) {
                n_init_cells = popsize;
                cerr << "n_init_cells = Pop size" << popsize << endl;
            } else {
                i--;
            }
        } else if (0 == strcmp(argv[i], "-casize")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in casize?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            sizex = atoi(argv[i]);
            i++;
            if (i == argc) {
                cerr << "Something odd in casize?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            sizey = atoi(argv[i]);
            cerr << "New value for CA size x and y: " << sizex << " " << sizey << endl;
        } else if (0 == strcmp(argv[i], "-n_nextgen")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in n_nextgen?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            howmany_makeit_for_nextgen = atoi(argv[i]);
            cerr << "New value for n_nextgen (howmany_makeit_for_nextgen in the code): " << howmany_makeit_for_nextgen
                 << endl;
        } else if (0 == strcmp(argv[i], "-season")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in season?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            season_duration = atoi(argv[i]);
            cerr << "New value for season (season_duration in the code): " << season_duration << endl;
        } else if (0 == strcmp(argv[i], "-nofood")) {
            is_there_food = false;
            cerr << "No food in this simulation" << endl;
        } else if (0 == strcmp(argv[i], "-noevolsim")) {
            evolsim = false;
            cerr << "No evolution in this simulation (sim ends when [howmany_makeit_for_nextgen] cells pass [the_line])"
                 << endl;
        } else if (0 == strcmp(argv[i], "-noscatter_start")) {
            scatter_start = false;
            cerr << "Cells will not be scattered at the start of the first season" << endl;
        } else if (0 == strcmp(argv[i], "-nochemgrad")) {
            chemgrad = false;
            cerr << "Not plotting chemotactic gradient" << endl;
        } else if (0 == strcmp(argv[i], "-groupextinction")) {
            groupextinction = true;
            cerr << "Stopping sim on group extinction" << endl;
        }  else if (0 == strcmp(argv[i], "-noevolvable_adh")) {
            evolvable_adh = false;
            cerr << "Adhesion is not evolvable" << endl;
        } else if (0 == strcmp(argv[i], "-nochemcircles")) {
            chemcircles = false;
            cerr << "Not plotting circles around food patches" << endl;
        } else if (0 == strcmp(argv[i], "-nodivisions")) {
            nodivisions = true;
            cerr << "Cells will not actually divide" << endl;
        } else if (0 == strcmp(argv[i], "-store")) {
            store = true;
            cerr << "pictures will be stored" << endl;
        } else if (0 == strcmp(argv[i], "-noevolreg")) {
            evolreg = false;
            cerr << "No evolution of regulation parameters" << endl;
        } else if (0 == strcmp(argv[i], "-metabperiod")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in metabperiod?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            metabperiod = atoi(argv[i]);
            cerr << "New value for metabperiod: " << metabperiod << endl;
        } else if (0 == strcmp(argv[i], "-gradscale")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in gradscale?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            gradscale = atof(argv[i]);
            cerr << "New value for gradscale: " << gradscale << endl;
        } else if (0 == strcmp(argv[i], "-foodpatches")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in foodpatches?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            foodpatches = atoi(argv[i]);
            cerr << "New value for foodpatches: " << foodpatches << endl;
        } else if (0 == strcmp(argv[i], "-foodpatchperiod")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in foodpatchperiod?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            foodpatchperiod = atoi(argv[i]);
            cerr << "New value for foodpatchperiod: " << foodpatchperiod << endl;
        } else if (0 == strcmp(argv[i], "-foodpatchlength")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in foodpatchlength?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            foodpatchlength = atoi(argv[i]);
            cerr << "New value for foodpatchlength: " << foodpatchlength << endl;
        } else if (0 == strcmp(argv[i], "-foodperspot")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in foodperspot?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            foodperspot = atoi(argv[i]);
            cerr << "New value for foodperspot: " << foodperspot << endl;
        } else if (0 == strcmp(argv[i], "-maxfoodpatches")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in maxfoodpatches?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            maxfoodpatches = atoi(argv[i]);
            cerr << "New value for maxfoodpatches: " << maxfoodpatches << endl;
        } else if (0 == strcmp(argv[i], "-foodstart")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in foodstart?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            foodstart = atoi(argv[i]);
            cerr << "New value for foodstart: " << foodstart << endl;
        } else if (0 == strcmp(argv[i], "-eatperiod")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in eatperiod?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            eatperiod = atoi(argv[i]);
            cerr << "New value for eatperiod: " << eatperiod << endl;
        } else if (0 == strcmp(argv[i], "-chemmu")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in chemmu?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            init_chemmu = atof(argv[i]);
            cerr << "New value for chemmu: " << init_chemmu << endl;
        } else if (0 == strcmp(argv[i], "-target_area")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in target_area?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            target_area = atoi(argv[i]);
            cerr << "New value for target_area: " << target_area << endl;
        } else if (0 == strcmp(argv[i], "-persmu")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in persmu?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            startmu = atof(argv[i]);
            cerr << "New value for persmu: " << startmu << endl;
        } else if (0 == strcmp(argv[i], "-init_cell_config")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in init_cell_config?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            init_cell_config = atoi(argv[i]);
            cerr << "New value for init_cell_config: " << init_cell_config << endl;
        } else if (0 == strcmp(argv[i], "-cell_placement")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in cell_placement?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            cell_placement = atoi(argv[i]);
            cerr << "New value for cell_placement: " << cell_placement << endl;
        } else if (0 == strcmp(argv[i], "-gradnoise")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in gradnoise?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            gradnoise = atof(argv[i]);
            cerr << "New value for gradnoise: " << gradnoise << endl;
        } else if (0 == strcmp(argv[i], "-name")) {
            i++;
            if (i == argc) {
                cerr << "Something odd in name?" << endl;
                return 1;  //check if end of arguments, exit with error in case
            }
            // I'm just going to work in c++ strings - a lot easier
            free(moviedir);
            free(latticedir);
            free(celldatadir);
            free(cellgravesdatadir);
            free(fooddatadir);

            string maybepath_and_name(argv[i]);
            size_t botDirPos = maybepath_and_name.find_last_of("/");
            string dir("");
            string name;
            if (botDirPos != std::string::npos) {
                // then there is a character '/' in name, which means that
                // we are going to save data in some path, hence
                // we have to split where this is happening
                dir = maybepath_and_name.substr(0, botDirPos + 1);
                name = maybepath_and_name.substr(botDirPos + 1, maybepath_and_name.length());
            } else {
                name = maybepath_and_name;
            }

            string name_outfile = dir; //will this contain the last '/''
            name_outfile.append("data_");
            name_outfile.append(name);
            name_outfile.append(".txt");

            string name_moviedir = dir;
            name_moviedir.append("movie_");
            name_moviedir.append(name);

            string name_latticedir = dir;
            name_latticedir.append("lattice_");
            name_latticedir.append(name);

            string name_celldatadir = dir;
            name_celldatadir.append("celldata_");
            name_celldatadir.append(name);

            string name_cellgravedatadir = dir;
            name_cellgravedatadir.append("cellgravedata_");
            name_cellgravedatadir.append(name);

            string name_fooddatadir = dir;
            name_fooddatadir.append("fooddata_");
            name_fooddatadir.append(name);

            string name_networkdir = dir;
            name_networkdir.append("networkdir_");
            name_networkdir.append(name);

            std::cerr << "New value for output filename: " << name_outfile << '\n';
            std::cerr << "New value for name_moviedir: " << name_moviedir << '\n';
            std::cerr << "New value for name_latticedir: " << name_latticedir << '\n';
            std::cerr << "New value for name_celldatadir: " << name_celldatadir << '\n';
            std::cerr << "New value for name_cellgravedatadir: " << name_cellgravedatadir << '\n';
            std::cerr << "New value for name_fooddatadir: " << name_fooddatadir << '\n';
            std::cerr << "New value for name_networkdir: " << name_networkdir << '\n';

            moviedir = (char *) malloc(50 + strlen(argv[i]) * sizeof(char));
            (char *) malloc(50 + strlen(argv[i]) * sizeof(char));
            latticedir = (char *) malloc(50 + strlen(argv[i]) * sizeof(char));
            celldatadir = (char *) malloc(50 + strlen(argv[i]) * sizeof(char));
            cellgravesdatadir = (char *) malloc(50 + strlen(argv[i]) * sizeof(char));
            fooddatadir = (char *) malloc(50 + strlen(argv[i]) * sizeof(char));
            strdup(name_outfile.c_str());
            moviedir = strdup(name_moviedir.c_str());
            latticedir = strdup(name_latticedir.c_str());
            celldatadir = strdup(name_celldatadir.c_str());
            cellgravesdatadir = strdup(name_cellgravedatadir.c_str());
            fooddatadir = strdup(name_fooddatadir.c_str());
            // this took a while to code :P
        } else {
            cerr << "Something went wrong reading the commandline arguments" << endl;
            return 1;
        }
    }
    return 0;
}

void Parameter::Read(const char *filename) {

    static bool ReadP = false;

    if (ReadP) {

        //throw "Run Time Error in parameter.cpp: Please Read parameter file only once!!";
        CleanUp();

    } else
        ReadP = true;

    FILE *fp = OpenReadFile(filename);


    T = fgetpar(fp, "T", 50., true);
    target_area = igetpar(fp, "target_area", 100, true);
    half_div_area = igetpar(fp, "half_div_area", 100, true);
    half_div_area_2 = igetpar(fp, "half_div_area_2", 100, true);
    target_length = igetpar(fp, "target_length", 60, true);
    lambda = fgetpar(fp, "lambda", 50, true);
    lambda2 = fgetpar(fp, "lambda2", 5.0, true);
    //Jtable = sgetpar(fp, "Jtable", "J.dat", true);
    conn_diss = igetpar(fp, "conn_diss", 2000, true);
    vecadherinknockout = bgetpar(fp, "vecadherinknockout", false, true);
    extensiononly = bgetpar(fp, "extensiononly", false, true);
    chemotaxis = igetpar(fp, "chemotaxis", 1000, true);
    border_energy = igetpar(fp, "border_energy", 100, true);
    neighbours = igetpar(fp, "neighbours", 2, true);
    min_area_for_life = igetpar(fp, "min_area_for_life", 5, true);
    Jmed = fgetpar(fp, "Jmed", 14, true);
    Jalpha = fgetpar(fp, "Jalpha", 7, true);
    key_lock_len = igetpar(fp, "key_lock_len", 6, true);
    n_chem = igetpar(fp, "n_chem", 0, true);
    if (n_chem) {
        diff_coeff = dgetparlist(fp, "diff_coeff", n_chem, true);
        decay_rate = dgetparlist(fp, "decay_rate", n_chem, true);
        secr_rate = dgetparlist(fp, "secr_rate", n_chem, true);
        saturation = fgetpar(fp, "saturation", 0, true);
        dt = fgetpar(fp, "dt", 2.0, true);
        dx = fgetpar(fp, "dx", 2.0e-6, true);
        pde_its = igetpar(fp, "pde_its", 15, true);
    }
    n_init_cells = igetpar(fp, "n_init_cells", 100, true);
    size_init_cells = igetpar(fp, "size_init_cells", 10, true);
    sizex = igetpar(fp, "sizex", 200, true);
    sizey = igetpar(fp, "sizey", 200, true);
    mcs = igetpar(fp, "mcs", 10000, true);
    rseed = igetpar(fp, "rseed", -1, true);
    subfield = fgetpar(fp, "subfield", 1.0, true);
    relaxation = igetpar(fp, "relaxation", 0, true);
    storage_stride = igetpar(fp, "storage_stride", 10, true);
    graphics = bgetpar(fp, "graphics", true, true);
    store = bgetpar(fp, "store", false, true);
    divisioncolour = bgetpar(fp, "divisioncolour", false, true);
    genomefile = sgetpar(fp, "genomefile", "", true);
    nodivisions = bgetpar(fp, "nodivisions", false, true);
    nr_regnodes = igetpar(fp, "nr_regnodes", 0, true);
    mu = fgetpar(fp, "mu", 0., true);
    mustd = fgetpar(fp, "mustd", 0., true);
    divtime = igetpar(fp, "divtime", 200, true);
    circle_thickness = igetpar(fp, "circle_thickness", 1, true);
    circle_segments = igetpar(fp, "circle_segments", 15, true);
    circle_dist = igetpar(fp, "circle_dist", 40, true);
    divdur = igetpar(fp, "divdur", 50, true);
    gompertz_alpha = fgetpar(fp, "gompertz_alpha", 0.0000075, true);
    gompertz_beta = fgetpar(fp, "gompertz_beta", 0.0004, true);
    scatter_start = bgetpar(fp, "scatter_start", true, true);
    chemgrad = bgetpar(fp, "chemgrad", true, true);
    groupextinction = bgetpar(fp, "groupextinction", false, true);
    evolvable_adh = bgetpar(fp, "evolvable_adh", true, true);
    chemcircles = bgetpar(fp, "chemcircles", true, true);
    motiledeath = fgetpar(fp, "motiledeath", 1.0, true);
    dividingdeath = fgetpar(fp, "dividingdeath", 0.0, true);
    moviedir = sgetpar(fp, "moviedir", "data_film", true);
    celldatafile = sgetpar(fp, "celldatafile", "", true);
    fooddatafile = sgetpar(fp, "fooddatafile", "", true);
    save_data_period = igetpar(fp, "save_data_period", 100, true);
    food_influx_location = sgetpar(fp, "food_influx_location", "nowhere", true);
    initial_food_amount = igetpar(fp, "initial_food_amount", 0, true);
    metabperiod = igetpar(fp, "metabperiod", 20, true);
    ardecay = fgetpar(fp, "ardecay", 0., true);
    growth = fgetpar(fp, "growth", 0., true);
    gradnoise = fgetpar(fp, "gradnoise", 0.1, true);
    gradscale = fgetpar(fp, "gradscale", 1.0, true);
    foodpatches = igetpar(fp, "foodpatches", 1, true);
    foodpatchperiod = igetpar(fp, "foodpatchperiod", 1000, true);
    foodpatchlength = igetpar(fp, "foodpatchlength", 1, true);
    foodperspot = igetpar(fp, "foodperspot", 1, true);
    maxfoodpatches = igetpar(fp, "maxfoodpatches", 10000, true);
    foodstart = igetpar(fp, "foodstart", 1000, true);
    eatperiod = igetpar(fp, "eatperiod", 1000, true);
    min_contact_duration_for_preying = igetpar(fp, "min_contact_duration_for_preying", 1., true);
    frac_contlen_eaten = fgetpar(fp, "frac_contlen_eaten", 1., true);
    metabolic_conversion = fgetpar(fp, "metabolic_conversion", 0.5, true);
    chancemediumcopied = fgetpar(fp, "chancemediumcopied", 0.0001, true);
    colortable_filename = sgetpar(fp, "colortable_filename", "default.ctb", true);
    plots = sgetpar(fp, "plots", "tau", true);
    key_lock_weights = sgetpar(fp, "key_lock_weights", "1 2 3 4 5 6", true);
    evolsim = igetpar(fp, "evolsim", 0, true);
    mut_rate = fgetpar(fp, "mut_rate", 0.01, true);
    circle_coverage = fgetpar(fp, "circle_coverage", 0.5, true);
    persduration = igetpar(fp, "persduration", 0, true);
    startmu = fgetpar(fp, "startmu", 0.0, true);
    init_chemmu = fgetpar(fp, "init_chemmu", 0.0, true);
    scaling_cell_to_ca_time = igetpar(fp, "scaling_cell_to_ca_time", 1, true);
    latticedir = sgetpar(fp, "latticedir", "lattice", true);
    celldatadir = sgetpar(fp, "celldatadir", "celldata", true);
    cellgravesdatadir = sgetpar(fp, "cellgravesdatadir", "cellgravedata", true);
    fooddatadir = sgetpar(fp, "fooddatadir", "fooddata", true);
    save_lattice_period = igetpar(fp, "save_lattice_period", 0, true);
    howmany_makeit_for_nextgen = igetpar(fp, "howmany_makeit_for_nextgen", 1, true);
    popsize = igetpar(fp, "popsize", 1, true);
    the_line = igetpar(fp, "the_line", 1, true);
    is_there_food = bgetpar(fp, "is_there_food", false, true);
    evolreg = bgetpar(fp, "evolreg", false, true);
    zero_persistence_past_theline = bgetpar(fp, "zero_persistence_past_theline", false, true);
    season_experiment = bgetpar(fp, "season_experiment", false, true);
    season_duration = igetpar(fp, "season_duration", 1, true);
    init_cell_config = igetpar(fp, "init_cell_config", 0, true);
    cell_placement = igetpar(fp, "cell_placement", 0, true);
}

// In the future the parser for the rules for key to J val tau,medium
// will be more developed, maybe even evolvable 8O
// int Parameter::SumLookupTableValue(int *lookup_table){
//   return -1;
// }
// int Parameter::MultiplyLookupTableValue(int *lookup_table){
//   return -1;
// }


const char *sbool(const bool &p) {

    const char *true_str = "true";
    const char *false_str = "false";
    if (p)
        return true_str;
    else
        return false_str;
}

void Parameter::Write(ostream &os) const {
    setlocale(LC_NUMERIC, "C");

    os << " T = " << T << endl;
    os << " target_area = " << target_area << endl;
    os << " div_area = " << half_div_area << endl;
    os << " div_area_2 = " << half_div_area_2 << endl;
    os << " target_length = " << target_length << endl;
    os << " lambda = " << lambda << endl;
    os << " lambda2 = " << lambda2 << endl;
    //if (Jtable)
    //  os << " Jtable = " << Jtable << endl;
    os << " conn_diss = " << conn_diss << endl;
    os << " vecadherinknockout = " << sbool(vecadherinknockout) << endl;
    os << " extensiononly = " << sbool(extensiononly) << endl;
    os << " chemotaxis = " << chemotaxis << endl;
    os << " border_energy = " << border_energy << endl;
    os << " neighbours = " << neighbours << endl;
    os << " min_area_for_life = " << min_area_for_life << endl;
    os << " Jmed = " << Jmed << endl;
    os << " Jalpha = " << Jalpha << endl;
    os << " key_lock_len = " << key_lock_len << endl;
    os << " n_chem = " << n_chem << endl;
    os << " diff_coeff = " << diff_coeff[0] << endl;
    os << " decay_rate = " << decay_rate[0] << endl;
    os << " secr_rate = " << secr_rate[0] << endl;
    os << " saturation = " << saturation << endl;
    os << " dt = " << dt << endl;
    os << " dx = " << dx << endl;
    os << " pde_its = " << pde_its << endl;
    os << " n_init_cells = " << n_init_cells << endl;
    os << " size_init_cells = " << size_init_cells << endl;
    os << " sizex = " << sizex << endl;
    os << " sizey = " << sizey << endl;
    os << " mcs = " << mcs << endl;
    os << " rseed = " << rseed << endl;
    os << " subfield = " << subfield << endl;
    os << " relaxation = " << relaxation << endl;
    os << " storage_stride = " << storage_stride << endl;
    os << " graphics = " << sbool(graphics) << endl;
    os << " store = " << sbool(store) << endl;
    if (genomefile) {
        os << " genomefile = " << genomefile << endl;
    }
    os << " nr_regnodes = " << nr_regnodes << endl;
    os << " mu = " << mu << endl;
    os << " mustd = " << mustd << endl;
    os << " divtime= " << divtime << endl;
    os << " circle_thickness= " << circle_thickness << endl;
    os << " circle_dist= " << circle_dist << endl;
    os << " circle_segments= " << circle_segments << endl;
    os << " divdur = " << divdur << endl;
    os << " gompertz_alpha = " << gompertz_alpha << endl;
    os << " gompertz_beta = " << gompertz_beta << endl;
    os << " initial_food_amount = " << initial_food_amount << endl;
    os << " food_influx_location = " << food_influx_location << endl;
    os << " metabperiod = " << metabperiod << endl;
    os << " ardecay = " << ardecay << endl;
    os << " growth = " << growth << endl;
    os << " gradnoise = " << gradnoise << endl;
    os << " gradscale = " << gradscale << endl;
    os << " foodpatches = " << foodpatches << endl;
    os << " foodpatchperiod = " << foodpatchperiod << endl;
    os << " foodpatchlength = " << foodpatchlength << endl;
    os << " foodperspot = " << foodperspot << endl;
    os << " maxfoodpatches = " << maxfoodpatches << endl;
    os << " foodstart = " << foodstart << endl;
    os << " eatperiod = " << eatperiod << endl;
    os << " divisioncolour = " << divisioncolour << endl;
    os << " nodivisions = " << nodivisions << endl;
    os << " motiledeath = " << motiledeath << endl;
    os << " dividingdeath = " << dividingdeath << endl;
    os << " min_contact_duration_for_preying = " << min_contact_duration_for_preying;
    os << " frac_contlen_eaten = " << frac_contlen_eaten << endl;
    os << " metabolic_conversion = " << metabolic_conversion << endl;
    os << " chancemediumcopied = " << chancemediumcopied << endl;
    os << " celldatafile = " << celldatafile << endl;
    os << " fooddatafile = " << fooddatafile << endl;
    os << " save_data_period = " << save_data_period << endl;
    os << " colortable_filename = " << colortable_filename << endl;
    os << " plots = " << plots << endl;
    os << " key_lock_weights = " << key_lock_weights << endl;
    os << " mut_rate = " << mut_rate << endl;
    os << " circle_coverage = " << circle_coverage << endl;
    os << " evolsim = " << evolsim << endl;
    os << " persduration = " << persduration << endl;
    os << " startmu = " << startmu << endl;
    os << " scaling_cell_to_ca_time = " << scaling_cell_to_ca_time << endl;
    os << " latticedir = " << latticedir << endl;
    os << " celldatadir = " << celldatadir << endl;
    os << " cellgravedatadir = " << cellgravesdatadir << endl;
    os << " fooddatadir = " << fooddatadir << endl;
    os << " save_lattice_period = " << save_lattice_period << endl;
    if (moviedir)
        os << " moviedir = " << moviedir << endl;
    os << " howmany_makeit_for_nextgen = " << howmany_makeit_for_nextgen << endl;
    os << " popsize = " << popsize << endl;
    os << " the_line = " << the_line << endl;
    os << " is_there_food = " << is_there_food << endl;
    os << " evolreg = " << evolreg << endl;
    os << " zero_persistence_past_theline = " << zero_persistence_past_theline << endl;
    os << " season_experiment = " << season_experiment << endl;
    os << " season_duration = " << season_duration << endl;
    os << " init_cell_config = " << init_cell_config << endl;
    os << " cell_placement = " << cell_placement << endl;
}

ostream &operator<<(ostream &os, Parameter &p) {
    p.Write(os);
    return os;
}

Parameter par;
