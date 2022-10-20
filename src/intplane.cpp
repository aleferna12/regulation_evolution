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
#include <math.h>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional> // for bind, see InitIncreaseVal()
#include "crash.h"
#include "parameter.h"
#include "ca.h"
#include "intplane.h"
#include "misc.h"

using std::placeholders::_1;

/* STATIC DATA MEMBER INITIALISATION */
const int IntPlane::nx[9] = {0, 1, 1, 1, 0,-1,-1,-1, 0 };
const int IntPlane::ny[9] = {0, 1, 0,-1,-1,-1, 0, 1, 1 };

extern Parameter par;


/** PRIVATE **/

IntPlane::IntPlane(const int sx, const int sy, int grad_srcs) {
  sigma=0;
  thetime=0;
  sizex=sx;
  sizey=sy;

  grad_sources = grad_srcs;
  // Needs more work, grad_sources scale too much with distance and not enough with area
  min_resource_dist = DetermineMinDist();
  cout << "MINDIST" << min_resource_dist << endl;
  diagonal = sqrt(sizex*sizex + sizey*sizey);
  peaksx = new int[grad_sources] {};
  peaksy = new int[grad_sources] {};

  // This should be always 1 if we keep the gradients linear
  dist_coef = 1;
  // This should be always 1 if we keep the gradients independent from each other
  interference = false;

  sigma=AllocateSigma(sx,sy);
}


IntPlane::IntPlane(void) {

  sigma=0;
  sizex=0; sizey=0;
  thetime=0;

}

// destructor (virtual)
IntPlane::~IntPlane(void) {
  if (sigma) {
    free(sigma[0]);
    free(sigma);
    sigma=0;
  }
}


// Alternatively, we could use the most isolated point to know where to put next peak at each iteration
// I think that doing this would be worse, as it is more computationally intensive and probably will tend to accumulate
// peaks in the corners (?), which may be problematic for small grad_sources numbers
double IntPlane::DetermineMinDist() {
  // Subtract two because I think positions 0 and size are forbidden (?) - yes, they seem to be
  double ratio = (sizey - 2) / (sizex - 2);
  // ratio * sepx = sepy
  // sepx = sepy / ratio
  // (sepx + 1) * (sepy + 1) = grad_sources - 1
  // ratio * pow(sepx, 2) + sepx * (1 + ratio) + 2 - grad_sources = 0
  // Do the same for sepy and solve quadradic equations
  double sepx = SolveQuadradic(ratio, 1 + ratio, 2 - grad_sources);
  double sepy = SolveQuadradic(1/ratio, 1 + 1/ratio, 2 - grad_sources);
  double mindistx = (sizex - 2) / (sepx * 2 + 2);
  double mindisty = (sizey - 2) / (sepy * 2 + 2);
  return sqrt(mindistx * mindistx + mindisty * mindisty);
}

int **IntPlane::AllocateSigma(const int sx, const int sy) {

  int **mem;
  sizex=sx; sizey=sy;

  mem=(int **)malloc(sizex*sizeof(int *));

  if (mem==NULL)
    MemoryWarning();

  mem[0]=(int *)malloc(sizex*sizey*sizeof(int));
  if (mem[0]==NULL)
      MemoryWarning();

  {  for (int i=1;i<sizex;i++)
    mem[i]=mem[i-1]+sizey;}

  /* Clear IntPlane plane */
  { for (int i=0;i<sizex*sizey;i++)
    mem[0][i]=0.; }

   return mem;
}

void IntPlane::Plot(Graphics *g2) {
  // l=layer: default layer is 0
  for (int x=1;x<sizex-1;x++)
    for (int y=1;y<sizey-1;y++) {
      // Make the pixel four times as large
      // to fit with the CPM plane
      g2->Point(sigma[x][y],2*x,2*y);
      g2->Point(sigma[x][y],2*x+1,2*y);
      g2->Point(sigma[x][y],2*x,2*y+1);
      g2->Point(sigma[x][y],2*x+1,2*y+1);
    }

}

// Plot the value of the intplane only in the medium of the CPM
void IntPlane::Plot(Graphics *g2, CellularPotts *cpm) {

  // this has to take into account stuff from cpm plane (maybe x,y info should give a tau of the cpm plane)

  // cpm->sigma[x][y] returns sigma, which I can use to indicise the vector of cells... can I?
  // not really, food doesn't know about cell

  // maybe this function should really be in dish...

  // suspend=true suspends calling of DrawScene
  for(int x=1;x<sizex-1;x++)
    for(int y=1;y<sizey-1;y++)
      //if (cpm->Sigma(x,y)==0) {
      if(sigma[x][y]!=0){
        int colorindex;
        if(sigma[x][y]>0) colorindex = 10+sigma[x][y];
        else colorindex = 5;
	      // Make the pixel four times as large
	      // to fit with the CPM plane
	      g2->Point(colorindex,2*x,2*y);
        g2->Point(colorindex,2*x+1,2*y);
        g2->Point(colorindex,2*x,2*y+1);
        g2->Point(colorindex,2*x+1,2*y+1);
      }
}

// private
void IntPlane::NoFluxBoundaries(void) {

  // all gradients at the edges become zero,
  // so nothing flows out
  // Note that four corners points are not defined (0.)
  // but they aren't used in the calculations


    for (int x=0;x<sizex;x++) {
      sigma[x][0]=sigma[x][1];
      sigma[x][sizey-1]=sigma[x][sizey-2];
    }

    for (int y=0;y<sizey;y++) {
      sigma[0][y]=sigma[1][y];
      sigma[sizex-1][y]=sigma[sizex-2][y];
    }

}


// private
void IntPlane::AbsorbingBoundaries(void) {

  // all boundaries are sinks,

    for (int x=0;x<sizex;x++) {
      sigma[x][0]=0.;
      sigma[x][sizey-1]=0.;
    }

    for (int y=0;y<sizey;y++) {
      sigma[0][y]=0.;
      sigma[sizex-1][y]=0.;
    }

}

// private
void IntPlane::PeriodicBoundaries(void) {

  // periodic...

    for (int x=0;x<sizex;x++) {
      sigma[x][0]=sigma[x][sizey-2];
      sigma[x][sizey-1]=sigma[x][1];
    }
    for (int y=0;y<sizey;y++) {
      sigma[0][y]=sigma[sizex-2][y];
      sigma[sizex-1][y]=sigma[1][y];
    }

}

void IntPlane::DiffuseParticles(void)
{
  double diff_const=0.01;
  std::vector<std::vector<int>> diff( sizex , std::vector<int>(sizey, 0));

  diff_const /= double(par.scaling_cell_to_ca_time);

  for(int i=1;i<sizex;i++)for(int j=1;j<sizey;j++){
    if(sigma[i][j]!=0){
      int foodtomove = BinomialDeviate( sigma[i][j] , diff_const );
      diff[i][j] -= foodtomove;
      while(foodtomove>0){
        //where does it move? randomly in the 8 neighbourhood (excl. self)
        int xpos = -1 + (int)( 3.*RANDOM() ); // int number in [-1,1]
        int ypos = -1 + (int)( 3.*RANDOM() ); // int number in [-1,1]
        xpos+=i;
        ypos+=j;
        if(par.periodic_boundaries){
          if(xpos>=sizex-1) xpos -= sizex-2;
          if(xpos<=0) xpos += sizex-2;
          if(ypos>=sizey-1) ypos -= sizey-2;
          if(ypos<=0) ypos += sizey-2;
          //std::cerr << "Hello3" << '\n';
          //std::cerr << "Hello4" << '\n';
        }else{
          // if fixed boundaries diffusion of particles on boundaries does not happen
          if(xpos>=sizex-1 || xpos<=0 || ypos>=sizey-1 || ypos<=0){
            xpos=i;
            ypos=j;
          }
        }
        diff[xpos][ypos]++;
        foodtomove--;
      }
    }
  }
  for(int i=1;i<sizex;i++)for(int j=1;j<sizey;j++){
    sigma[i][j]+=diff[i][j];
  }
}
//copy of this function in ca.cpp
int IntPlane::SetNextVal(int pos){
  //the plane has a 1 px boundary on all size, therefore we place the pixels
  //within that boundary
  static int xcount=1, ycount=1;

  if(xcount>=sizex-1 ||ycount>=sizey-1){
    return 1;
  }

  sigma[xcount][ycount]=pos;
  ycount++;
  if(ycount==sizey-1){
    ycount=1;
    xcount++;
  }
  return 0;
}

// This function initialises a functional (from <function>)
// depending on parameters. The function is responsible for updating the field
// If 'nowhere' option is used, parameter food_influx should be set
void IntPlane::InitIncreaseVal(CellularPotts *cpm) {

  // change values of initial_food_amount and others
  // if food_influx_location== "nowhere"
  if(strcmp(par.food_influx_location,"specified_experiment") == 0){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    IncreaseVal = std::bind(&IntPlane::IncreaseValSpecifiedExp, this, cpm);
    //     exit(1);
    ;
  }
  else{
    cerr<<"INIT: Error. Got unidentified food influx location: "<<par.food_influx_location<<endl;
    exit(1);
  }

  //exit(1);
}


peakinfo IntPlane::ClosestPeak(int x, int y, int upto) {
  if (upto == -1) {
    upto = grad_sources;
  }

  peakinfo res;
  long mindist_sq = (long) INFINITY;

  for (int src = 0; src < upto; src++) {
    // - ensures always finding smaller value
    int dx = peaksx[src] - x;
    int dy = peaksy[src] - y;
    long dist_sq = dx * dx + dy * dy;
    if (dist_sq < mindist_sq) {
      mindist_sq = dist_sq;
      res.x = peaksx[src];
      res.y = peaksy[src];
    }
  }

  res.dist = sqrt(mindist_sq);
  return res;
}


double IntPlane::DistMostIsolatedPoint() {
  double dist = 0;
  for(int i=1;i<sizex-1;i++)for(int j=1;j<sizey-1;j++) {
    double closest_dist = ClosestPeak(i, j).dist;
    if (closest_dist > dist) {
      dist = closest_dist;
    }
  }
  return dist;
}


void IntPlane::RandomizeResourcePeaks() {
  peaksx[0] = (int) RandomNumber(sizex - 1);
  peaksy[0] = (int) RandomNumber(sizey - 1);
  for (int src = 1; src < grad_sources; src++) {
    int x = 0;
    int y = 0;
    double dist = 0;
    while (dist < min_resource_dist) {
      x = (int) RandomNumber(sizex - 1);
      y = (int) RandomNumber(sizey - 1);
      dist = ClosestPeak(x, y, src).dist;
    }
    peaksx[src] = x;
    peaksy[src] = y;
  }
  dist_most_isolated = DistMostIsolatedPoint();
}


double IntPlane::FoodEquation(double dist_from_peak) {
  // Prevents negative value generation from resource sources too far away
  dist_from_peak = min(dist_from_peak, dist_most_isolated);
  return par.gradscale * dist_most_isolated/100 * pow(1 - dist_from_peak/dist_most_isolated, dist_coef);
}


double IntPlane::FoodAtPosition(int x, int y) {
  // double pfood_j = 0.125;

  // makes gradient
  // int maxfood = 3;
  // int maxfood = 1+5.* (1. - dist_from_peak/(double)sizey);

  //This is how it was before, worked for field size of 500
  // double dfood = 1+5.* (1. - dist_from_peak/(double)sizey); //this the usable line
  // so maybe - to standardize gradients across field sizes, I could do:
  // dfood = 1 + sizey/100 * (1. - dist_from_peak/(double)sizey)
  // so that the local slope of the gradient stays the same?
  // also- the 1+ part of the equation could go...
  // or even better counter balanced by a lesser gradient in the variable part
  double dfood = 0;
  // TODO: Solve the interference (can it be safely implemented as parameter? Do we even want that?)
  if (not interference) {
    double dist_from_peak = ClosestPeak(x, y).dist;
    dfood = FoodEquation(dist_from_peak);
  } else {
    for (int src = 0; src < grad_sources; src++) {
      int dx = peaksx[src] - x;
      int dy = peaksy[src] - y;
      double dist_from_peak = sqrt(dx * dx + dy * dy);
      dfood += FoodEquation(dist_from_peak);
    }
  }
  // Make sure grad is never 0 (has to do with builtin color generation)
  dfood++;
  return dfood;
}


[[deprecated]]
int IntPlane::WritePeaksData() {
  int num_rows = 5;
  ofstream file;
  file.open(par.peaksdatafile);
  for (int i = 0; i < num_rows; i++) {
    int row = i * sizex / num_rows + 1;
    cout << row << endl;
    for (int col = 1; col < sizey - 1; col++) {
      file << sigma[row][col] << ",";
    }
    file << endl;
  }
  file.close();
  return 0;
}


// I am going to change the direction of the gradient every so often
void IntPlane::IncreaseValSpecifiedExp(CellularPotts *cpm)
{
  RandomizeResourcePeaks();

  maxfood = 0;
  for(int i=1;i<sizex-1;i++)for(int j=1;j<sizey-1;j++){
    double dfood = FoodAtPosition(i, j);
    int local_maxfood = (int)dfood;
    sigma[i][j] = local_maxfood;
    if(RANDOM() < dfood - local_maxfood) local_maxfood++;
    if(RANDOM() < par.gradnoise)
      sigma[i][j]=local_maxfood;
    if (sigma[i][j] > maxfood)
      maxfood = local_maxfood;
      //sigma[i][j]+=10; //else already set to zero
    // else
    //   sigma[i][j]=0;
    //
    // if(i>sizex/2+25 || i<sizex/2-25) {
    //    sigma[i][j]=0;
    //    continue;
    // }
    // //int foodhere;
    // double dist_from_peak;
    // dist_from_peak= (sizey-j)/(double)(sizey); //linear gradient
    // // dist_from_peak= sqrt( (sizey-j)*(sizey-j) + (sizex/2-i)*(sizex/2-i) );
    // int local_maxfood = 3+7.* dist_from_peak/sizey;
    // double pfood_j = 0.5+ 0.5* (sizey-j)/(double)(sizey);
    // if(RANDOM() < pfood_j) sigma[i][j]=local_maxfood;
    // else sigma[i][j]=0;

    //bool is_there_food = false;
    if(par.is_there_food){
      if(RANDOM()<par.foodinflux) sigma[i][j]=-1; //food
    }
  }
  // If we want to see transversal profiles of the plane
  // WritePeaksData();
}
