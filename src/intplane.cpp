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
#include "conrec.h"

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
  min_resource_dist = sqrt(sizex*sizey) / (sqrt(grad_sources) + 1);
  diagonal = sqrt(sizex*sizex + sizey*sizey);
  peaksx = new int[grad_sources] {};
  peaksy = new int[grad_sources] {};

  // This should be always 1 if we keep the gradients linear
  dist_coef = 1;
  // This should be always 1 if we keep the gradients independent from each other
  interference = false;

  // TODO: Remove
  peakx = 1;
  peaky = sizey/2;

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
  if( strcmp(par.food_influx_location,"everywhere") == 0 ){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    IncreaseVal = std::bind(&IntPlane::IncreaseValEverywhere, this);

    //exit(1);
  }else if(strcmp(par.food_influx_location,"nowhere") == 0){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    par.foodinflux=0.; //this is not necessary, because NotIncreaseVal is an empty function, bu well...
    IncreaseVal = std::bind(&IntPlane::NotIncreaseVal, this);
//     exit(1);
    ;
  }else if(strcmp(par.food_influx_location,"notonprey") == 0){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    IncreaseVal = std::bind(&IntPlane::IncreaseValIfEmpty, this, cpm);
//     exit(1);
    ;
  }else if(strcmp(par.food_influx_location,"patchy") == 0){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    IncreaseVal = std::bind(&IntPlane::IncreaseValPatchy, this, cpm);
    //     exit(1);
    ;
  }else if(strcmp(par.food_influx_location,"somewhere") == 0){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    IncreaseVal = std::bind(&IntPlane::IncreaseValSomewhere, this, cpm);
    //     exit(1);
    ;
  }else if(strcmp(par.food_influx_location,"somewhere_notonprey") == 0){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    IncreaseVal = std::bind(&IntPlane::IncreaseValSomewhereIfEmpty, this, cpm);
    //     exit(1);
    ;
  }else if(strcmp(par.food_influx_location,"patchy_random") == 0){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    IncreaseVal = std::bind(&IntPlane::IncreaseValPatchyRandom, this, cpm);
    //     exit(1);
    ;
  }else if(strcmp(par.food_influx_location,"patchy_random_persistence") == 0){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    IncreaseVal = std::bind(&IntPlane::IncreaseValPatchyRandomPersistence, this, cpm);
    //     exit(1);
    ;
  }else if(strcmp(par.food_influx_location,"food_growth") == 0){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    IncreaseVal = std::bind(&IntPlane::IncreaseValSelfGrowth, this, cpm);
    //     exit(1);
    ;
  }else if(strcmp(par.food_influx_location,"specified_experiment") == 0){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    IncreaseVal = std::bind(&IntPlane::IncreaseValSpecifiedExp, this, cpm);
    //     exit(1);
    ;
  }
  else if(strcmp(par.food_influx_location,"boundarygradient") == 0){
    cerr<<"Hello, got food influx location: "<<par.food_influx_location<<endl;
    IncreaseVal = std::bind(&IntPlane::IncreaseValBoundaryGrad, this, cpm);
    //     exit(1);
    ;
  }
  else{
    cerr<<"INIT: Error. Got unidentified food influx location: "<<par.food_influx_location<<endl;
    exit(1);
  }

  //exit(1);
}

// we check beforeahand how many grid points we update and then radnomly generate coordinates.
// Typically, scaling_cell_to_ca_time>1, so foodinflux/scaling_cell_to_ca_time < 1/2,
// so we update less than half of the pixels, so this method is better
void IntPlane::IncreaseValEverywhere(void)
{
  static double prob_food_influx = par.foodinflux/ double(par.scaling_cell_to_ca_time);
  //if(prob_food_influx)<0.5{
    int totalpixels= (sizex-1)*(sizey-1);
    double howmany_pixels_updated = BinomialDeviate( totalpixels, prob_food_influx );
    for(int i=0; i<howmany_pixels_updated; i++){
      int posi = RandomNumber(sizex-2); //in the interval [1,sizex-1[
      int posj = RandomNumber(sizey-2);
      if(sigma[posi][posj]<10) sigma[posi][posj]++;
    }
//   }else{
//     for (int i=1; i<sizex-1; i++)
//       for (int j=1; j<sizey-1;j++){
//         if(RANDOM()<prob_food_influx && sigma[i][j]<10) //replace random() with proper function
//           sigma[i][j]++;
//       }
//   }


}

void IntPlane::NotIncreaseVal(void)
{
  return;
}

void IntPlane::IncreaseValIfEmpty(CellularPotts *cpm)
{
  static double prob_food_influx = par.foodinflux/ double(par.scaling_cell_to_ca_time);
  int totalpixels= (sizex-1)*(sizey-1);
  double howmany_pixels_updated = BinomialDeviate( totalpixels, prob_food_influx );
  for(int i=0; i<howmany_pixels_updated; i++){
    int posi = RandomNumber(sizex-2); //in the interval [1,sizex-1[
    int posj = RandomNumber(sizey-2);
    if(sigma[posi][posj]<10 && !cpm->Sigma(posi,posj)) sigma[posi][posj]++;
  }

//     for (int i=1; i<sizex-1; i++)
//         for (int j=1; j<sizey-1;j++)
//         {
//             if(RANDOM()<par.foodinflux && sigma[i][j]<10 && !cpm->Sigma(i,j)) //replace random() with proper function
//                 sigma[i][j]++;
// //                 cerr<<"Hello food increased"<<endl;
//         }
}

void IntPlane::IncreaseValPatchy(CellularPotts *cpm)
{
  static double prob_food_influx = par.foodinflux/ double(par.scaling_cell_to_ca_time);
  //if(prob_food_influx)<0.5{
  int totalpixels= (sizex-1)*(sizey-1);
  double howmany_pixels_updated = BinomialDeviate( totalpixels, prob_food_influx );
  for(int i=0; i<howmany_pixels_updated; i++){
    int posi = RandomNumber(sizex-2); //in the interval [1,sizex-1] both included
    int posj = RandomNumber(sizey-2);
    for(int ii=-1;ii<=1;ii++)for(int jj=-1;jj<=1;jj++){
      int posii=posi+ii;
      int posjj=posj+jj;
      if(posii >= sizex-1 || posii < 1) continue;
      if(posjj >= sizey-1 || posjj < 1) continue;
      if(sigma[posi+ii][posj+jj]<10) sigma[posi+ii][posj+jj]=10;
    }
  }
  //   }else{
  //     for (int i=1; i<sizex-1; i++)
  //       for (int j=1; j<sizey-1;j++){
  //         if(RANDOM()<prob_food_influx && sigma[i][j]<10) //replace random() with proper function
  //           sigma[i][j]++;
  //       }
  //   }


}

// food increases everwhere, but somewhere it's a lot, always there
// this is a 3x3 grid
void IntPlane::IncreaseValSomewhere(CellularPotts *cpm){

//   std::vector<int> meanposx = {sizex/8, sizex/4, 3*sizex/8, sizex/2, 5*sizex/8, 3*sizex/4, 7*sizex/8};
//   std::vector<int> meanposy = {sizey/8, sizey/4, 3*sizey/8, sizey/2, 5*sizey/8, 3*sizey/4, 7*sizey/8};
//

  std::vector<int> meanposx = {sizex/4, sizex/2, 3*sizex/4};
  std::vector<int> meanposy = {sizey/4, sizey/2, 3*sizey/4};

//   std::vector<int> meanposx = {sizex/4,  3*sizex/4};
//   std::vector<int> meanposy = {sizey/4,  3*sizey/4};

//   std::vector<int> meanposx = {sizex/2};
//   std::vector<int> meanposy = {sizey/2};

  static double prob_food_influx = par.foodinflux/ double(par.scaling_cell_to_ca_time);
  //if(prob_food_influx)<0.5{
  int totalpixels= (sizex-1)*(sizey-1);
  double howmany_pixels_updated = BinomialDeviate( totalpixels, prob_food_influx );
  for(int i=0; i<howmany_pixels_updated; i++){
    int posi = RandomNumber(sizex-2); //in the interval [1,sizex-1] both included
    int posj = RandomNumber(sizey-2);
    if(sigma[posi][posj]<10) sigma[posi][posj]++;
  }

  for(auto mpx: meanposx){
    for(int k=-10; k<10;k++){
      for(auto mpy: meanposy){
        for(int l=-10; l<10;l++){
          if( mpx+k >= sizex-1 || mpx+k < 1 || mpy+l >= sizey-1 || mpy+l < 1 ) continue;

//           if( sigma[ mpx+k ][ mpy+l ]<10 && RANDOM()<par.foodinflux/double(par.scaling_cell_to_ca_time) ) {
          if( sigma[ mpx+k ][ mpy+l ]<10 && RANDOM()<0.2/double(par.scaling_cell_to_ca_time) ) {
//             cerr<<"Updating: "<< meanposx[i]+k <<" "<< meanposy[i]+l <<endl;
            sigma[ mpx+k ][ mpy+l ]++;
          }

        }
      }
    }

  }

}

void IntPlane::IncreaseValPatchyRandom(CellularPotts *cpm){

  // standard trickling-in of food
  static double prob_food_influx = par.foodinflux/ double(par.scaling_cell_to_ca_time);
  //if(prob_food_influx)<0.5{
  int totalpixels= (sizex-1)*(sizey-1);
  double howmany_pixels_updated = BinomialDeviate( totalpixels, prob_food_influx );
  for(int i=0; i<howmany_pixels_updated; i++){
    int posi = RandomNumber(sizex-2); //in the interval [1,sizex-1] both included
    int posj = RandomNumber(sizey-2);
    if(sigma[posi][posj]<10) sigma[posi][posj]++;
  }

  DiffuseParticles();

  double food_patch_frequency = 0.002;

  if( RANDOM() < food_patch_frequency/double(par.scaling_cell_to_ca_time) ){
    int mpx = RandomNumber(sizex-2); //mean position x of patch
    int mpy = RandomNumber(sizey-2);
    for(int k=-100; k<100;k++){
      for(int l=-100; l<100;l++){
        int distsquare= k*k+l*l;
        if(distsquare > 100*100) continue; //makes a circle
        int px=mpx+k;
        int py=mpy+l;
        //boundary wrapping
        if( mpx+k >= sizex-1) px -= (sizex-1);
        else if(mpx+k < 1) px += (sizex-1);
        if (mpy+l >= sizey-1) py -= (sizey-1);
        else if( mpy+l < 1 ) py += (sizey-1);

        if( distsquare<25*25)
          sigma[ px ][ py ] = 10;
        else if( distsquare< 50*50  && RANDOM()< 0.25 )
          sigma[ px ][ py ] = 10;
        else if( distsquare< 75*75  && RANDOM()< 0.5 )
          sigma[ px ][ py ] = 10;
        else if(RANDOM()< 0.25)
          sigma[ px ][ py ] =  10;

        // if( distsquare>75*75  && RANDOM()< 0.25 )
        //   sigma[ px ][ py ] +=  (10. - sigma[ px ][ py ])*RANDOM();
        // else if( distsquare>50*50  && RANDOM()< 0.5 )
        //   sigma[ px ][ py ] +=  (10. - sigma[ px ][ py ])*RANDOM();
        // else if( distsquare>25*25  && RANDOM()< 0.75 )
        //   sigma[ px ][ py ] +=  (10. - sigma[ px ][ py ])*RANDOM();
        //
        // else if( distsquare>10*10  && RANDOM()< 0.75 )
        // sigma[ px ][ py ] +=  (10. - sigma[ px ][ py ])*RANDOM();
        // else
        //   sigma[ px ][ py ] =  10;

      }
    }
  }

}

//like above, but lower increase rate, and longer time
void IntPlane::IncreaseValPatchyRandomPersistence(CellularPotts *cpm){
  int food_update_persistence=10000; //every 1000 time steps food will update somewhere else
  static int time_since_update=0;
  // double food_patch_frequency = 0.02;
  double food_patch_prob=0.2;
  int food_radius =25;
  int food_radius_square =food_radius*food_radius;

  // static vector<int> v_time_since_update={0,0,0};
  // static vector<int> meanposx = {0,0,0};
  // static vector<int> meanposy = {0,0,0};
  //

  // int mpx,mpy;
  static vector<vector<int>> v_t_since_update_mposxy(3, vector<int>(3,0) );


  //
  // static int mpx=0,mpy=0;

  // standard trickling-in of food
  static double prob_food_influx = par.foodinflux/ double(par.scaling_cell_to_ca_time);
  //if(prob_food_influx)<0.5{
  int totalpixels= (sizex-1)*(sizey-1);
  double howmany_pixels_updated = BinomialDeviate( totalpixels, prob_food_influx );
  for(int i=0; i<howmany_pixels_updated; i++){
    int posi = RandomNumber(sizex-2); //in the interval [1,sizex-1] both included
    int posj = RandomNumber(sizey-2);
    if(sigma[posi][posj]<10) sigma[posi][posj]++;
  }

  for(vector<vector<int>>::iterator it = v_t_since_update_mposxy.begin() ; it != v_t_since_update_mposxy.end(); ++it){
    if( (*it)[0] % food_update_persistence == 0){
      (*it)[0]= food_update_persistence/2 + RANDOM()*food_update_persistence/2; //5000 + random(0,5000) timesteps
      (*it)[1]= RandomNumber(sizex-2);
      (*it)[2]= RandomNumber(sizey-2);

    }
  // }


  // if(time_since_update%food_update_persistence == 0){
  //   time_since_update=0;
  //   mpx = RandomNumber(sizex-2); //mean position x of patch
  //   mpy = RandomNumber(sizey-2);
  //
  // }


  // if( RANDOM() < food_patch_frequency/double(par.scaling_cell_to_ca_time) ){
    int mpx = (*it)[1];
    int mpy = (*it)[2];


    for(int k=-food_radius; k<food_radius;k++){
      for(int l=-food_radius; l<food_radius;l++){
        int distsquare= k*k+l*l;
        if(distsquare > food_radius*food_radius) continue; //makes a circle
        int px=mpx+k;
        int py=mpy+l;
        //boundary wrapping
        if( mpx+k >= sizex-1) px -= (sizex-1);
        else if(mpx+k < 1) px += (sizex-1);
        if (mpy+l >= sizey-1) py -= (sizey-1);
        else if( mpy+l < 1 ) py += (sizey-1);

        if( distsquare<food_radius_square/(4*4) && RANDOM()< food_patch_prob/double(par.scaling_cell_to_ca_time) ){
          if(sigma[ px ][ py ] <10) sigma[ px ][ py ] ++;}
        else if( distsquare< food_radius_square/(2*2)  && RANDOM()< 0.25*food_patch_prob/double(par.scaling_cell_to_ca_time) ) {
          if(sigma[ px ][ py ] <10) sigma[ px ][ py ] ++;}
        else if( distsquare< food_radius_square/((4./3.) * (4./3.) )  && RANDOM()< 0.5*food_patch_prob/double(par.scaling_cell_to_ca_time) ) {
          if(sigma[ px ][ py ] <10) sigma[ px ][ py ] ++;}
        else if(RANDOM()< 0.25*food_patch_prob/double(par.scaling_cell_to_ca_time)){
          if(sigma[ px ][ py ] <10) sigma[ px ][ py ] ++;}

        // if( distsquare>75*75  && RANDOM()< 0.25 )
        //   sigma[ px ][ py ] +=  (10. - sigma[ px ][ py ])*RANDOM();
        // else if( distsquare>50*50  && RANDOM()< 0.5 )
        //   sigma[ px ][ py ] +=  (10. - sigma[ px ][ py ])*RANDOM();
        // else if( distsquare>25*25  && RANDOM()< 0.75 )
        //   sigma[ px ][ py ] +=  (10. - sigma[ px ][ py ])*RANDOM();
        //
        // else if( distsquare>10*10  && RANDOM()< 0.75 )
        // sigma[ px ][ py ] +=  (10. - sigma[ px ][ py ])*RANDOM();
        // else
        //   sigma[ px ][ py ] =  10;

      }
    }

    (*it)[0]++; // increase time in the counter
  }



}

//the idea is a population of bacteria like things that selfreplicate locally
void IntPlane::IncreaseValSelfGrowth(CellularPotts *cpm)
{
  double bacteria_birthrate=0.01;
  for(int i=1;i<sizex;i++)for(int j=1;j<sizey;j++){
    if(sigma[i][j]!=0){
      // there are so many bacteria, each with their birth rate, so:
      int new_borns = BinomialDeviate( sigma[i][j] , bacteria_birthrate/double(par.scaling_cell_to_ca_time) );
      if(new_borns==0) continue;

      while(new_borns>0){
        //where does it replicate? randomly in the 9 neighbourhood (including self)
        int xpos = -1 + (int)( 3.*RANDOM() ); // int number in [-1,1]
        int ypos = -1 + (int)( 3.*RANDOM() ); // int number in [-1,1]
        xpos+=i;
        ypos+=j;
        if(xpos>=sizex-1) xpos -= sizex-2;
        if(xpos<=0) xpos += sizex-2;
        if(ypos>=sizey-1) ypos -= sizey-2;
        if(ypos<=0) ypos += sizey-2;
        //std::cerr << "Hello3" << '\n';
        if(sigma[xpos][ypos]<2) sigma[xpos][ypos]++;
        //std::cerr << "Hello4" << '\n';

        new_borns--;
      }
    }

  }
  // return;
}

//a gradient emanating from one of the boundaries (north)
void IntPlane::IncreaseValBoundaryGrad(CellularPotts *cpm)
{
  int maxfood;
  double pfood_j, dfood;

  peakx=sizex/2;
  peaky=1;

  for(int i=1;i<sizex-1;i++)for(int j=1;j<sizey-1;j++){
    sigma[i][j]=0;

    dfood = par.gradscale*( (double)sizey/100.) * (1. - j/(double)sizey); //this the usable line

    maxfood=(int)dfood;
    if(RANDOM() < dfood - maxfood) maxfood++; //finer gradient made with a little unbiased noise

    //maxfood = 1+5.* (1 - (double)j/par.gradlength);
    // pfood_j =par.gradnoise+ (1.-par.gradnoise)* (1 - (double)j/(double)sizey);
    pfood_j = par.gradnoise;
    if(RANDOM() < pfood_j)  sigma[i][j]=maxfood;
  }
}


double IntPlane::DistClosestResourcePeak(int x, int y, int upto) {
  if (upto == -1) {
    upto = grad_sources;
  }
  long mindist_sq = (long) INFINITY;
  for (int src = 0; src < upto; src++) {
    // - ensures always finding smaller value
    int dx = peaksx[src] - x;
    int dy = peaksy[src] - y;
    long dist_sq = dx * dx + dy * dy;
    mindist_sq = min(dist_sq, mindist_sq);
  }
  return sqrt(mindist_sq);
}


double IntPlane::DistMostIsolatedPoint() {
  double dist = 0;
  for(int i=1;i<sizex-1;i++)for(int j=1;j<sizey-1;j++) {
    double closest_dist = DistClosestResourcePeak(i, j);
    if (closest_dist > dist) {
      dist = closest_dist;
    }
  }
  return dist;
}


void IntPlane::RandomizeResourcePeaks() {
  peaksx[0] = (int) RandomNumber(sizex);
  peaksy[0] = (int) RandomNumber(sizey);
  for (int src = 1; src < grad_sources; src++) {
    int x = 0;
    int y = 0;
    double dist = 0;
    while (dist < min_resource_dist) {
      x = (int) RandomNumber(sizex);
      y = (int) RandomNumber(sizey);
      dist = DistClosestResourcePeak(x, y, src);
    }
    peaksx[src] = x;
    peaksy[src] = y;
  }
  dist_most_isolated = DistMostIsolatedPoint();
}


double IntPlane::FoodEquation(double dist_from_peak) {
  // Prevents negative value generation from resource sources too far away
  dist_from_peak = min(dist_from_peak, dist_most_isolated);
  return par.gradscale * diagonal/100 * pow(1 - dist_from_peak/dist_most_isolated, dist_coef);
}


double IntPlane::FoodAtPostition(int x, int y) {
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
  if (not interference) {
    double dist_from_peak = DistClosestResourcePeak(x, y);
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


int IntPlane::MaxFood() {
  int maxfood = 0;
  for(int i=1;i<sizex-1;i++)for(int j=1;j<sizey-1;j++){
    maxfood = max(maxfood, Sigma(i, j));
  }
  return maxfood;
}


// I am going to change the direction of the gradient every so often
void IntPlane::IncreaseValSpecifiedExp(CellularPotts *cpm)
{
  RandomizeResourcePeaks();
  for(int i=1;i<sizex-1;i++)for(int j=1;j<sizey-1;j++){
    double dfood = FoodAtPostition(i, j);
    sigma[i][j] = dfood;
    int maxfood = (int)dfood;
    if(RANDOM() < dfood - maxfood) maxfood++; //finer gradient made with a little unbiased noise

    // noise
    // double pfood_j = 0.1+ 0.9* (1. - dist_from_peak/(double)sizey);  // this is the usable one
    //double pfood_j = 1.;
    //pfood_j = 0.3; // also like this is works... but why?

    //sigma[i][j]+=maxfood;
    double pfood_j = par.gradnoise;

    if(RANDOM() < pfood_j)
      sigma[i][j]=maxfood;
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
    // int maxfood = 3+7.* dist_from_peak/sizey;
    // double pfood_j = 0.5+ 0.5* (sizey-j)/(double)(sizey);
    // if(RANDOM() < pfood_j) sigma[i][j]=maxfood;
    // else sigma[i][j]=0;

    //bool is_there_food = false;
    if(par.is_there_food){
      if(RANDOM()<par.foodinflux) sigma[i][j]=-1; //food
    }
  }
  WritePeaksData();
}


void IntPlane::IncreaseValSomewhereIfEmpty(CellularPotts *cpm){

  //   std::vector<int> meanposx = {sizex/8, sizex/4, 3*sizex/8, sizex/2, 5*sizex/8, 3*sizex/4, 7*sizex/8};
  //   std::vector<int> meanposy = {sizey/8, sizey/4, 3*sizey/8, sizey/2, 5*sizey/8, 3*sizey/4, 7*sizey/8};
  //

  std::vector<int> meanposx = {sizex/4, sizex/2, 3*sizex/4};
  std::vector<int> meanposy = {sizey/4, sizey/2, 3*sizey/4};

  //   std::vector<int> meanposx = {sizex/4,  3*sizex/4};
  //   std::vector<int> meanposy = {sizey/4,  3*sizey/4};

  //   std::vector<int> meanposx = {sizex/2};
  //   std::vector<int> meanposy = {sizey/2};

  static double prob_food_influx = par.foodinflux/ double(par.scaling_cell_to_ca_time);
  //if(prob_food_influx)<0.5{
  int totalpixels= (sizex-1)*(sizey-1);
  double howmany_pixels_updated = BinomialDeviate( totalpixels, prob_food_influx );
  for(int i=0; i<howmany_pixels_updated; i++){
    int posi = RandomNumber(sizex-2); //in the interval [1,sizex-1] both included
    int posj = RandomNumber(sizey-2);
    if(sigma[posi][posj]<10 && !cpm->Sigma(posi,posj)) sigma[posi][posj]++;
  }

  for(auto mpx: meanposx){
    for(int k=-10; k<10;k++){
      for(auto mpy: meanposy){
        for(int l=-10; l<10;l++){
          if( mpx+k >= sizex-1 || mpx+k < 1 || mpy+l >= sizey-1 || mpy+l < 1 ) continue;

          //           if( sigma[ mpx+k ][ mpy+l ]<10 && RANDOM()<par.foodinflux/double(par.scaling_cell_to_ca_time) ) {
          if( sigma[ mpx+k ][ mpy+l ]<10 && RANDOM()<0.2/double(par.scaling_cell_to_ca_time) && !cpm->Sigma(mpx+k,mpy+l)) {
            //             cerr<<"Updating: "<< meanposx[i]+k <<" "<< meanposy[i]+l <<endl;
            sigma[ mpx+k ][ mpy+l ]++;
          }

        }
      }
    }

  }

}
