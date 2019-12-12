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
#include <list>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <map>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include "cell.h"
#include "sticky.h"
#include "parameter.h"
#include "dish.h"

#define HASHCOLNUM 255

extern Parameter par;

int **Cell::J=0;
int Cell::amount=0;
int Cell::capacity=0;
int Cell::maxsigma=0;
int Cell::maxtau=0;

//Cell::Cell(const Dish &who) : Cytoplasm(who);
// Note: g++ wants to have the body of this constructor in cell.hh
// body is defined in "ConstructorBody" below
class Dish;

using namespace std;

Cell::~Cell(void) {

  amount--;
  if (amount==0) {
    // clear J if last cell has been destructed
    free(J[0]);
    free(J);
    capacity=0;
    maxsigma=0;
    J=0;
  }
  if(par.n_chem){
    delete[] chem;
  }
}


void Cell::CellBirth(Cell &mother_cell) {

  colour = mother_cell.colour;
  alive = mother_cell.alive;
  v[0] = mother_cell.v[0];
  v[1] = mother_cell.v[1];

  // Administrate ancestry
  mother_cell.daughter=this->sigma;
  mother=mother_cell.sigma;
  times_divided=++mother_cell.times_divided;
  owner=mother_cell.owner;

  date_of_birth=owner->Time();
  //cerr<<"sigma:"<<sigma<<" I am born at: "<<date_of_birth<<endl<<endl; //this works fine

  colour_of_birth=mother_cell.colour;
  colour=mother_cell.colour;

  alive=mother_cell.alive;

  tau=mother_cell.tau;
  target_area = mother_cell.target_area;
  target_length = mother_cell.target_length;
  half_div_area = mother_cell.half_div_area;

  genome=src.genome;
  genome.MutateGenome();

  // Do not add moments here, they are going to be calculated from scratch
  // and are initialised to zero in ConstructorBody(), called right before this
  meanx=mother_cell.meanx;
  meany=mother_cell.meany;

  startTarVec();  //why this?
  mother_cell.startTarVec(); //randomises mother cell target vector
                             // this fixes the biased movement bug!
                             // also biologically it makes sense,
                             // because cell polarity is all screwed after cell division
  //tvecx=mother_cell.tvecx;
  //tvecy=mother_cell.tvecy;

  prevx=mother_cell.prevx; //in startTarVec prevx is set to meanx, so we re-set this here
  prevy=mother_cell.prevy;
  persdur=mother_cell.persdur;
  perstime=int(persdur*RANDOM()); //assign random persistence duration upon birth.
  mu=mother_cell.mu;

  for (int ch=0;ch<par.n_chem;ch++)
    chem[ch]=mother_cell.chem[ch];

  n_copies=0;

  chemmu=mother_cell.chemmu;
  chemvecx=mother_cell.chemvecx;
  chemvecy=mother_cell.chemvecy;

  grad[0]=mother_cell.grad[0];
  grad[1]=mother_cell.grad[1];

  particles=mother_cell.particles/2; //daugher gets half particles of mother
  mother_cell.particles = particles; // mother gets particled halved

  growth=mother_cell.growth;
  eatprob=mother_cell.eatprob;

  maintenance_fraction = mother_cell.maintenance_fraction; //must be copied for first time step of life
  k_mf_0 = mother_cell.k_mf_0;
  k_mf_A = mother_cell.k_mf_A;
  k_mf_P = mother_cell.k_mf_P;
  k_mf_C = mother_cell.k_mf_C;

  extprotexpress_fraction = mother_cell.extprotexpress_fraction;
  k_ext_0 = mother_cell.k_ext_0;
  k_ext_A = mother_cell.k_ext_A;
  k_ext_P = mother_cell.k_ext_P;
  k_ext_C = mother_cell.k_ext_C;
  k_ext_0t= mother_cell.k_ext_0t;
  k_ext_Pt= mother_cell.k_ext_Pt;


  weight_for_chemotaxis = mother_cell.weight_for_chemotaxis;
  k_chem_0=mother_cell.k_chem_0;
  k_chem_A=mother_cell.k_chem_A;
  k_chem_P=mother_cell.k_chem_P;
  k_chem_C=mother_cell.k_chem_C;

  clearNeighbours(); //neighbours will be reassigned during the division function

  jlock = mother_cell.jlock;
  jkey = mother_cell.jkey;
  vJ = mother_cell.vJ;

  time_since_birth=0;
  mother_cell.SetTimeSinceBirth(0);
}


void Cell::ConstructorBody(int settau,int setrecycledsigma) {

  //cout<<"Tomato 2, do you come here?"<<endl;

  // Note: Constructor of Cytoplasm will be called first
  alive=true;
  colour=1; // undifferentiated

  colour_of_birth=1;
  date_of_birth=0;
  times_divided=0;
  mother=0;
  daughter=0;

  // add new elements to each of the dimensions of "J"


  // maxsigma keeps track of the last cell identity number given out to a cell
  if(setrecycledsigma==-1){
    // amount gives the total number of Cell instantiations (including copies)
    amount++;

    sigma=maxsigma++;
    //cout<<"check: not recycling, new sigma is "<<sigma<<endl;
  }
  else{
    sigma=setrecycledsigma;
    //cout<<"check: YES recycling, new sigma is "<<sigma<<endl;
  }

  //if (!J) {
  //  ReadStaticJTable(par.Jtable);
  //}

  // This should not be here
  //if(!vJ){
  //  ReadKeyLockFromFile(par.keylock_list_filename)
  //}

  //create the genome: use parameters for size. This is a randomly generated genome.
  //genome.InitGenome(par.innodes, par.regnodes, par.outnodes);
  //let's do this separately.
  tau=settau;
  area=0;
  target_area=0;
  half_div_area=0;
  length=0;
  target_length=par.target_length;
  sum_x=0;
  sum_y=0;
  sum_xx=0;
  sum_yy=0;
  sum_xy=0;
  particles=0;
  eatprob=0.;
  growth=par.growth;

  maintenance_fraction = 1;
  k_mf_0 = par.init_k_mf_0;
  k_mf_A = par.init_k_mf_A;
  k_mf_P = par.init_k_mf_P;
  k_mf_C = par.init_k_mf_C;

  extprotexpress_fraction = 1;
  k_ext_0 = par.init_k_ext_0;
  k_ext_A = par.init_k_ext_A;
  k_ext_P = par.init_k_ext_P;
  k_ext_C = par.init_k_ext_C;
  k_ext_0t = par.init_k_ext_0t;
  k_ext_Pt = par.init_k_ext_Pt;


  weight_for_chemotaxis=0.;
  k_chem_0=par.init_k_chem_0;
  k_chem_A=par.init_k_chem_A;
  k_chem_P=par.init_k_chem_P;
  k_chem_C=par.init_k_chem_C;

  //  growth_threshold=par.dthres;
  growth_threshold=0;
  v[0]=0.; v[1]=0.;
  n_copies=0;
  mu=0.0;
  tvecx=0.;
  tvecy=0.;
  prevx=0.;
  prevy=0.;

  chemvecx=0.;
  chemvecy=0.;

  time_since_birth=0;

  persdur=0;
  perstime=0;
  if(par.n_chem){
    chem = new double[par.n_chem];
  }

  //cerr<<"sigma="<<sigma <<" date of birth: "<<date_of_birth<<endl; //this works fine
}


/*! \brief Read a table of static Js.
 First line: number of types (including medium)
 Next lines: diagonal matrix, starting with 1 element (0 0)
 ending with n elements */
void Cell::ReadStaticJTable(const char *fname) {

  cerr << "Reading J's...\n";
  ifstream jtab(fname);
  if (!jtab)
    perror(fname);

  int n; // number of taus
  jtab >> n;
  cerr << "Number of celltypes:" <<  n << endl;
  maxtau=n-1;

  // Allocate
  if (J) { free(J[0]); free(J); }
  J=(int **)malloc(n*sizeof(int *));
  J[0]=(int *)malloc(n*n*sizeof(int));
  for (int i=1;i<n;i++) {
    J[i]=J[i-1]+n;
  }

  capacity = n;
  {for (int i=0;i<n;i++) {
    for (int j=0;j<=i;j++) {
      jtab >> J[i][j];
      // symmetric...
      J[j][i]=J[i][j];
    }

  }}
}


int Cell::EnergyDifference(const Cell &cell2) const
{

  if (sigma==cell2.sigma){
    //cerr<<"EnergyDifference(): Warning. sigma and sigma2 are the same"<<endl;
    return 0;
  }
  //return J[tau][cell2.tau]; //old version
  //cerr<<"Hello EnergyDifference 0.1: this is vJ of cell with sigma = "<<sigma<<endl;
  //for (auto i: vJ)
  //  cerr << i << " ";
  //cerr<<endl;

  //cerr<<"In contrast, this is what J[tau][cell2.tau] would give: "<<J[tau][cell2.tau]<<endl;
  /*if(vJ[cell2.sigma]<=1){
    cerr<<"tau="<<tau<<" sigma="<<sigma<<", tau2="<<cell2.tau<<" sigma2="<<cell2.sigma<<", J=" <<vJ[cell2.sigma]<<"=?"<<cell2.vJ[sigma]<<endl;
  */
//     for (auto i = jkey.begin(); i != jkey.end(); ++i)
//       std::cout << *i ;
//     std::cout << ' ';
//     for (auto i = jlock.begin(); i != jlock.end(); ++i)
//       std::cout << *i ;
//     std::cout<<endl;
//     for (auto i = cell2.jlock.begin(); i != cell2.jlock.end(); ++i)
//       std::cout << *i ;
//     std::cout << ' ';
//     for (auto i = cell2.jkey.begin(); i != cell2.jkey.end(); ++i)
//       std::cout << *i ;
//     std::cout<<endl;
//
//
//     exit(1);
//   }
  return vJ[cell2.sigma]; //new version
}

void Cell::ClearJ(void) {

  for (int i=0;i<capacity*capacity;i++) {
    J[0][i]=EMPTY;
  }
}


void Cell::setNeighbour(int neighbour, int boundarylength, int contactduration)
{

  if(boundarylength==0)//remove this neighbour
    neighbours.erase(neighbour);
  else
    neighbours[neighbour]=make_pair(boundarylength, contactduration); //if the element is already present, the boundarylength will be modified, otherwise a new element will be created.

}

int Cell::returnBoundaryLength(int cell)
{
  if(neighbours.count(cell))
    return neighbours[cell].first;

  return 0;
}


int Cell::returnDuration(int cell)
{
  if(neighbours.count(cell))
    return neighbours[cell].second;

  return 0;
}

void Cell::clearNeighbours()
{
  neighbours.clear();
}

int Cell::updateNeighbourBoundary(int cell, int boundarymodification)
{
  //cerr<<"Hello updNeiBound begin"<<endl;

  if(!neighbours.count(cell) && boundarymodification<0){
    printf("Cell.updateNeighbourBoundary: error: negatively updating contact of cell %d with nonexisting neighbour %d\n",sigma,cell);
    exit(1);
    //return 1;
  }else if(!neighbours.count(cell)){
    //cerr<<"Hello updNeiBound 0"<<endl;
    neighbours[cell]=make_pair(boundarymodification,0);
    //cerr<<"Hello updNeiBound 1"<<endl;
  }else if(neighbours.count(cell)){
    //cerr<<"Hello updNeiBound 2"<<endl;
    neighbours[cell].first += boundarymodification;
    //cerr<<"Hello updNeiBound 3"<<endl;
    //DO this after one MCS for duration reasons
//     if(neighbours[cell].first==0)//remove this neighbour
//     {
//       neighbours.erase(cell);
//     }
  }

  if(neighbours[cell].first<0){
    neighbours.erase(cell);
    printf("Cell.updateNeighbourBoundary: error: updating contact of cell %d with neighbour %d to negative value\n",sigma,cell);
    return 2;
  }

  return 0;
}

int Cell::SetNeighbourDurationFromMother(int cell, int motherduration)
{
  if( !neighbours.count(cell) ){
    printf("Cell.SetNeighbourDurationFromMother: error: Nonexisting neighbour %d\n",cell);
    return 1;
  }
  else
    neighbours[cell].second = motherduration;

  return 0;
}

int Cell::updateNeighbourDuration(int cell, int durationmodification)
{
  if(!neighbours.count(cell))
  {
    printf("Cell.updateNeighbourDuration: error: Nonexisting neighbour %d\n",cell);
    return 1;
  }

  else if(neighbours.count(cell))
  {
    neighbours[cell].second += durationmodification;
  }


  return 0;
}

int Cell::MutateKeyAndLock(void)
{
//   cerr<<"key: ";
//   for(auto x: jkey) cerr<<x<<" ";
//   cerr<<"lock: ";
//   for(auto x: jlock) cerr<<x<<" ";
//   cerr<<endl;


  //double mutrate=0.2; // probability per bit to be flipped
  int nmut;
  //calculate how many mutations from binomial distribution for both key and lock
  int keysize = jkey.size();
  int locksize = jlock.size();
  int keylocksize= keysize + locksize;

//   cerr<<"mut rate is: "<<par.mut_rate<<endl;

  nmut = BinomialDeviate(keylocksize,par.mut_rate);
  //if zero, return 0
  if(0==nmut) return 0;

  //else randomize positions, and flip bits
  //return number of mutations
  int positions[keylocksize];
  for(int i=0; i<keylocksize; i++) positions[i]=i; //initialise array of positions
  int where=keylocksize;
  for(int i=0; i<nmut; i++){
    int mutpos = (int)(where*RANDOM());
    int tmp = positions[mutpos];
    positions[mutpos]=positions[where-1];
    positions[where-1] = tmp;
    where--;

    //key or lock
    if(positions[keylocksize -i -1 ] < keysize){
      int pos = positions[keylocksize-i-1];
      jkey[ pos ] = 1 - jkey[ pos ]; // 0 if 1, and 1 if 0
    }else{
      int pos=positions[keylocksize -i -1 ] - keysize;
      jlock[ pos ] = 1 - jlock[ pos ]; // 0 if 1, and 1 if 0
    }

  }


//   cerr<<"Ney: ";
//   for(auto x: jkey) cerr<<x<<" ";
//   //cerr<<endl;
//   cerr<<"Nock: ";
//   for(auto x: jlock) cerr<<x<<" ";
//   cerr<<endl;
//
//   cerr<<"nmut = "<<nmut<<endl;
//   for(int i=keylocksize-1; i> keylocksize - nmut-1;i--) cerr<<positions[i]<<" ";
//   cerr<<endl;
//   exit(1);
//
  return nmut;

}

//returns a number between 0 and 1 which is either the maintenance_fraction
// or the fraction of expressed surface proteins (for adhesion)
double Cell::CalculateMaintenance_or_ExtProtExpr_Fraction(double k0, double kA,double kP,double kC)
{
  double fraction=k0;
  // std::cerr << "k0,kA,kP,kC = "<<k0<<" "<<kA<<" "<<kP<<" "<<kC  << '\n';

  // std::cerr << "k0 = "<<fraction << '\n';
  fraction += kA * area / (double)half_div_area;

  // std::cerr << "k0+kA*A = "<<fraction  << '\n';

  fraction += kP * particles / 50.; // <-a reasonable scaling factor :P

  // std::cerr << "k0+kA*P+kP*P = "<<fraction  << '\n';


  //next bit averages over non zero contacts:
  // I should get the Jvalues from the sigma in contact with this cell
  // and the length
  // should get a normalised number in [0,1] when I divided by a reasonble scaling factor: 43 :P
  double contlen_total=0.;
  double sum_J_times_contlen=0.;
  for(auto nei:neighbours){
    int sigma_nei = nei.first; // sigma of cell in contact
    if(sigma_nei==0) continue;

    int contlen_nei = nei.second.first; // contact length
    contlen_total += contlen_nei;
    sum_J_times_contlen += vJ[sigma_nei]*contlen_nei; //
  }

  //no contribution from contacts if this cell is in contact with no one
  double toadd_avrgJ = (contlen_total>0)? sum_J_times_contlen/(double)(contlen_total*43.): 0.;
  fraction += kC * toadd_avrgJ;

  // std::cerr << "k0+kA*P+kP*P+kC*C = "<<fraction  << '\n';


  if(fraction<0.) fraction=0.;
  if(fraction>1.) fraction=1.;

  // std::cerr << "Final = "<<fraction  << '\n';


  return fraction;
}
