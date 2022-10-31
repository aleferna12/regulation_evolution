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
#include <stdio.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <math.h>
#include "dish.h"
#include "random.h"
#include "cell.h"
#include "info.h"
#include "parameter.h"
#include "sqr.h"
#include "output.h"

#ifdef QTGRAPHICS
#include "qtgraph.h"
#else
#include "x11graph.h"
#endif
#include <chrono>
using namespace std::chrono;

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

    if (! strlen(par.backupfile) && !strlen(par.competitionfile)) {
      //THIS IS TO USE FOR NORMAL INITIALISATION
      if(par.scatter_start){
        CPM->PlaceCellsRandomly(par.n_init_cells,par.size_init_cells);
      }else{
        CPM->PlaceCellsOrderly(par.n_init_cells,par.size_init_cells);
      }
      CPM->InitializeEdgeList(false);
      cout << "done initialising edge list"<<endl;

      CPM->ConstructInitCells(*this); //within an object, 'this' is the object itself

      // Assign a random type to each of the cells, i.e. PREYS and PREDATORS
      CPM->SetRandomTypes();
      cout << "done setting types"<<endl;
      //Initialise key-lock pairs - we do it after types, because we need the information
      InitKeyLock();
      cout << "done setting keylock"<<endl;
      //Initialise vector of J values for each cell
      InitVectorJ();
      cout << "done setting J values"<<endl;
      //Initialise the contactlength bookkeeping now that the cells are placed
      // at this stage, cells are only surrounded by medium
      InitContactLength();  // see dish.cpp - you don't need dish->InitContactLength because this part IS in dish
      cout << "done setting contact length"<<endl;

      cout << "Going to initialise genome"<<endl;
      for(auto &c: cell) {
        if(c.Sigma()){
          c.setGTiming((int)(RANDOM()*par.scaling_cell_to_ca_time));
          //c.setGTiming(0);
          c.dividecounter=0;
          c.SetTargetArea(par.target_area); //sets target area because in dividecells the new target area = area
          //initialise a cell's timing for gex Updating
          //c.setGTiming((int)(RANDOM()*par.scaling_cell_to_ca_time));
          //creates a cell's genome, either randomly or from file
          if (strlen(par.genomefile)){
            c.ReadGenomeFromFile(par.genomefile);
          }
          else{
            c.CreateRandomGenome(2, par.nr_regnodes, 1); //important! we only have one output node, but used to make 2 (one unused)
          }
          c.ClearGenomeState();
        }
      }

      cout << "Done initialising genome"<<endl;

      //Set function pointer for food update, depending on parameters
      Food->InitIncreaseVal(CPM); //a pointer to CPM is an argument to InitIncreaseVal
                                     // but NOT to IncreaseVal if it points to IncreaseValEverywhere

      // Initialises food plane (now the gradient plane)
      Food->RandomizeResourcePeaks();
      Food->IncreaseVal(*(Food));
      cout<<"done with food"<<endl;
      //run CPM for some time without persistent motion
      for(int init_time=0;init_time<10;init_time++){
        CPM->AmoebaeMove2(PDEfield);  //this changes neighs
      }
      cout<<"done with update"<<endl;
      InitCellMigration();
      UpdateCellParameters3(0);//update cell status //UpdateCellParameters2();
      par.starttime=0;
    }
    else if (strlen(par.competitionfile)){
      ReadCompetitionFile(par.competitionfile);
      CPM->InitializeEdgeList(false);
      cout << "done initialising edge list for competition"<<endl;
      InitContactLength();
      cout << "done initialising contacts for competition"<<endl;
      InitVectorJ();
      cout << "done initialising Jvalues for competition"<<endl;
      Food->InitIncreaseVal(CPM); //a pointer to CPM is an argument to InitIncreaseVal
      // Initialises food plane (now the gradient plane)
      Food->RandomizeResourcePeaks();
      Food->IncreaseVal(*(Food));
      for(int init_time=0;init_time<10;init_time++){
        CPM->AmoebaeMove2(PDEfield);  //this changes neighs
      }
      cout<<"done with update"<<endl;
      InitCellMigration();
      UpdateCellParameters3(0);//update cell status //UpdateCellParameters2();
      par.starttime=0;
    }

    else {
      cout << "Reading backfile" << endl;
      cout<<"backup file is "<<par.backupfile<<endl;
      par.starttime=ReadBackup(par.backupfile);
      int networktime=par.starttime;
      //now get the right network into the cells. first find which time point to use
      while(networktime%par.season_duration){
        networktime+=1000;
      }
      char fname[300];
      for(auto &c: cell) {
        if(c.Sigma() && c.AliveP()){
          //c.setGTiming((int)(RANDOM()*par.scaling_cell_to_ca_time));
          //c.dividecounter=0;
          //c.SetTargetArea(par.target_area); //sets target area because in dividecells the new target area = area
          sprintf(fname,"%s/t%010d_c%04d.txt",par.genomefile,networktime,c.Sigma());
          c.ReadGenomeFromFile(fname);
        }
      }
      CPM->InitializeEdgeList(false);
      InitContactLength();
      InitVectorJ();
      Food->InitIncreaseVal(CPM);
      Food->IncreaseVal(*(Food));
    }
  } catch(const char* error) {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);
  }
}

TIMESTEP {

  try {
    static Dish *dish=new Dish(); //here ca planes and cells are constructed
    static Info *info=new Info(*dish, *this);
    static int i=par.starttime; //starttime is set in Dish. Not the prettiest solution, but let's hope it works.
    static int sum=0, nr=0;
    int extinct=0;

    if( !(i%100000) ) cerr<<"TIME: "<<i<<endl;

    //auto start = high_resolution_clock::now();
    dish->CellsEat4();
    //auto stop = high_resolution_clock::now();
    //auto duration = duration_cast<microseconds>(stop - start);
    //sum+=duration.count();
    //nr++;
    //cout << duration.count() << endl;
    //This function updates the network and deals with the consequences of the output (motility vs division)

    dish->UpdateCellParameters3(i); // for continuous GRN updating and reproduction

    dish->CellMigration();//updates persistence time and targetvectors

    dish->CPM->AmoebaeMove2(dish->PDEfield);  //this changes neighs

    dish->UpdateNeighDuration();

    //Change the season and do evolution
    if( i%25 == 0){

      if(par.evolsim){
        if(i>par.starttime && i%par.season_duration==0 ){
          // cout <<" duration average "<<sum/nr<<endl;
          //reproduce people based on fitness criterion
          //remove random cells until popsize is back to normal
          //reset food and gradient
          std::cerr << "Time = "<<i << '\n';
          std::cerr << "End of season: there are "<< dish->CountCells() <<" cells" << '\n';
          dish->SaveNetworks(i);
          dish->SaveAdheringNeighbours(i);
          dish->GradientBasedCellKill(par.popsize);
          dish->SaveAncestry(i);
          //dish->RemoveMotileCells(par.popsize); //kill all nondividing cells and more if necessary; for noncontinuous reproduction
          std::cerr << "After remove there are "<< dish->CountCells() <<" cells" << '\n';

          if (par.scatter_cells) dish->ScatterEndOfSeason();
        //  dish->ReproduceEndOfSeason(); //replenish the population; for noncontinuous sim
          //cout<<"done reproducing"<<endl;
          //dish->UpdateCellParameters2();//update cell status; for noncontinuous sim
          //cout<<"done updating"<<endl;
          dish->Food->RandomizeResourcePeaks();
          dish->Food->IncreaseVal(*(dish->Food)); //this has to be last thing to do here
          std::cout << "End of season: Gradient switching at time (+/- 25 MCS) = "<< i << '\n';
          if(strlen(par.competitionfile)){
            //check if one of the groups is extinct. if yes, end simulation
            extinct=dish->CountCellGroups();
            
            if(extinct){
              std::cout << "Group extinct after "<< i <<" time steps. ending simulation..."<< endl;
              exit(0);
            }
          }
        }
      }else{
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
        if (i>par.starttime && i%par.season_duration==0){
          exit(0);
        }
      }
    }


    // TO SCREEN
    // UNUSED
    if (par.graphics && !(i%par.storage_stride)) {

      BeginScene();
      ClearImage();
      if(par.readcolortable){
      //  dish->Plot(this,1);
        dish->Plot(this,2);
      }
      else{
        dish->Plot(this,0);
      }
      //dish->Food->Plot(this, dish->CPM);
      //char title[400];
      //snprintf(title,399,"CellularPotts: %d MCS",i);
      //ChangeTitle(title);
      EndScene();
      info->Menu();
    }
    // TO FILE FOR MOVIE
    if (par.store && !(i%par.storage_stride)) {
      if(par.readcolortable){
        char fname[300];
        sprintf(fname,"%s/angle%09d.png",par.datadir,i);
        BeginScene(); //this is an empty function for X11
        ClearImage(); //
        dish->Plot(this,1); //everything contained here
        EndScene();
        Write(fname);
        sprintf(fname,"%s/order%09d.png",par.datadir,i);
        BeginScene(); //this is an empty function for X11
        ClearImage(); //
        dish->Plot(this,2); //everything contained here
      //dish->Food->Plot(this,dish->CPM); //will this work?  YES !!!
        EndScene();
        Write(fname); //FIXED SO THAT CODE AND IMAGE MATCH!
      }else{
        char fname[300];
        sprintf(fname,"%s/tau%09d.png",par.datadir,i);
        // BeginScene(); //this is an empty function for X11
        ClearImage();
        dish->Plot(this,0); // this is g //everything contained here
        EndScene();
        Write(fname);
      }
    }
    // TO FILE FOR TEXT
    if( !(i%par.save_text_file_period) ){
      int popsize=dish->SaveData(i); //saves data to text
      if( 0 == popsize ){
        cerr << "Global extinction after"<<i<<"time steps, simulation terminates now" << endl;
        exit(0);
      }
    }
    // TO FILE FOR BACKUP
    if( !(i%par.save_backup_period) ){
      dish->MakeBackup(i); //saves all permanent data
    }

    i++;
  } catch(const char* error) {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);
  }

}

int PDE::MapColour(double val) {

  return (((int)((val/((val)+1.))*100))%100)+155;
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
    if(argc>2){
      int exit_valarg = par.ReadArguments(argc,argv);
      if( 0 != exit_valarg ){
        par.PrintWelcomeStatement(); //see parameter.h
        exit(1);
      }
    }

    //Creates a rule for J val with medium that is going to be used in following function
    par.CreateRule(par.Jmed_rule_input);
    // Open file where initial key locks are specified, and assigns them to cells
    par.Read_KeyLock_list_fromfile(par.keylock_list_filename);

    cerr<<endl<<"Warning, this version is ***NOT*** suitable for pde field!!!"<<endl;
    //Depends on this: AddSiteToMoments (and Remove), FindCellDirections2, etc...
    cerr<<endl<<"WARNING, use wrapped boundaries if cells are A LOT smaller than sizex and sizey"<<endl<<endl;
    cerr<<endl<<"WARNING: DO NOT EVOLVE CHEMMU, or if you do, change the replication function (where it is always reset to init_chemmu)"<<endl<<endl;

    //check if directory for movies exists, create it if not, exit otherwise
    DoesDirExistsIfNotMakeit(par.datadir);  //see output.cpp
    DoesDirExistsIfNotMakeit(par.backupdir);  //see output.cpp
    DoesDirExistsIfNotMakeit(par.networkdir);  //see output.cpp

    //check if data file exists, if not exit
    if(FileExistsP(par.datafile)){
      cerr<<"File "<< par.datafile << " already exists, simulation not starting" << endl;
      exit(1);
    }

    if(par.periodic_boundaries && par.lambda2>0.){
      cerr<<"main(): Error. Cannot have wrapped periodic boundaries and lambda2>0"<<endl;
      cerr<<"(because I cannot calculate second moment for cells crossing boundaries)"<<endl;
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
    cerr <<"Using X11 graphics (batch mode). sizex and y are "<<par.sizex<<" "<< par.sizey <<endl;
    X11Graphics g(par.sizex*2,par.sizey*2);
    int t;

    for(t=0;t<=par.mcs;t++){
      //cerr<<"Time: "<<t<<endl;
      g.TimeStep();

    }
#endif

  }catch(const char* error){
    std::cerr << error << "\n";
    exit(1);
  }
  catch(...) {
    std::cerr << "An unknown exception was caught\n";
  }
  return 0;
}
